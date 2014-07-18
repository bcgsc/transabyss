"""
This module provides methods for finding read-support for events detected
through the Trans-ABySS pipeline.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import sys
import re
import math
import operator
import pysam
from intspan import overlap, subsume
from tools import reverse_complement, seq_overlap, remove_from_list

class BAM:
    """Wrapper using PySam to extract alignment info from BAM file,
    provides generalized methods to extract read support for chimeric events
    """
    def __init__(self, bam_file, min_mapq=1, sample_type=None):
	self.source = bam_file
        self.bam = pysam.Samfile(bam_file, "rb")
        self.min_mapq = min_mapq
	self.sample_type = sample_type

    def confirm(self, target, pos=[], feature=None, allele=None, regions=None, seq=False, breakpoint=None, breakpoint_buffer=0, refseq=None, focal=False, max_depth=None, no_proper=False, ctg_len=None):
	"""delegator for finding read-support"""
        feature = feature.lower()
        reads = []
        
        if feature == 'snv' or feature == 'inv':
            reads = self.confirm_snv(target, pos[0], allele, seq=seq)
        elif feature == 'del':
            reads = self.confirm_del(target, pos, seq=seq, refseq=refseq)
	    
	# 'seq' should always be False, other 'big' insertions won't be validated
        elif feature == 'ins':
            reads = self.confirm_ins(target, pos[0], allele, refseq=refseq, seq=seq)
	    
        elif feature == 'breakpoint':
            if not regions:
                reads = self.confirm_breakpoint_contig(target, pos, seq=seq, max_depth=max_depth, ctg_len=ctg_len)
            else:
                reads = self.confirm_breakpoint_genomic(regions, seq=seq, breakpoint=breakpoint, breakpoint_buffer=breakpoint_buffer, focal=focal, no_proper=no_proper)
		
        return reads

    def confirm_perfect(self, target, pos, allele, strand, seq=False, min_from_end=0, expansion=0):
        """Confirms every base in the target is the same as allele
	This is an old method, probably needs to be re-written
	"""
	# makes sure coordinates are integers
	pos[0] = int(pos[0])
	pos[1] = int(pos[1])
	# makes sure coordinate starts from 1
	pos[0] = max(pos[0], 1)
	
        num_reads = 0
        reads = []
	reads_dict = {}
	
	# re-construct allele if this is a repeat expansion
	if expansion > 1:
	    i = 1
	    seq = allele
	    while i < expansion:
		allele = allele + seq
		i = i + 1
	
	# for alleles longer than read length
	coverage = {}
	for i in range(pos[0], pos[1]+1):
	    coverage[i] = 0
	    	
	shorter_than_read = False
	read_lens = {}
        for read in self.bam.fetch(target, pos[0], pos[1] + 1):
	    if len(allele) <= read.rlen:
		shorter_than_read = True
		if read.pos+1 <= pos[0] and read.pos+read.rlen >= pos[1]:
		    if not read_lens.has_key(read.rlen):
			read_lens[read.rlen] = 0
		    read_lens[read.rlen] += 1
		    
		    from_end = min(pos[0] - read.pos -1, read.pos + read.rlen - pos[1])
		    offset = pos[0]-read.pos-1
		    bases = read.seq[offset:offset+len(allele)]
		
		    # needs to reverse-complement bases if contig is aligned to - strand, 
		    # as allele is given as + strand of reference
		    if strand == '-':
			bases = reverse_complement(bases)
		
		    if from_end < min_from_end:
			continue
		    
		    elif bases.lower() == allele.lower():
			if seq:
			    reads.append(read)
			num_reads += 1
		    		
	read_len = None
	if read_lens:
	    sorted_read_lens = sorted(read_lens.iteritems(), key=operator.itemgetter(1))
	    read_len = sorted_read_lens[0][0]
	    
	# none of the reads is longer than the allele(insertion), or no read within pos
	# not sure condition: len(allele) > read_len / 2?
	if not shorter_than_read or (not reads and read_len is not None and len(allele) > read_len / 2):
	    allele_seq = allele
	    if strand == '-':
		allele_seq = reverse_complement(allele)
		
	    # check if allele is covered by reads
	    for pileupcolumn in self.bam.pileup(target, pos[0]-1, pos[1]):
		if pileupcolumn.pos >= pos[0]-1 and pileupcolumn.pos <= pos[1]-1:
		    cov = pileupcolumn.n
		    
		    for pileupread in pileupcolumn.pileups:
			if read_len is None:
			    if not read_lens.has_key(pileupread.alignment.rlen):
				read_lens[pileupread.alignment.rlen] = 0
			    read_lens[pileupread.alignment.rlen] += 1
			    			
			# pysam bug?
			if pileupread.qpos >= len(pileupread.alignment.seq):
			    continue
						
			if pileupread.alignment.seq[pileupread.qpos].lower() != allele_seq[pileupcolumn.pos-pos[0]+1]:
			    cov = cov - 1
			else:
			    reads_dict[pileupread.alignment.qname] = pileupread.alignment
		    coverage[pileupcolumn.pos + 1] = cov
		    			    	    
	    # determine if allele is entirely covered
	    gaps = False
	    for i in range(pos[0], pos[1]+1):
		if coverage[i] == 0:
		    gaps = True
		    break
		
	    if read_len is None and read_lens:
		sorted_read_lens = sorted(read_lens.iteritems(), key=operator.itemgetter(1))
		read_len = sorted_read_lens[0][0]
		
	    if not gaps and (read_len is not None and len(allele) > read_len / 2):
		if seq:
		    reads = reads_dict.values()
		else:
		    num_reads = len(reads_dict.keys())
	
        if seq:
            return reads
        else:
            return num_reads
                               
    def confirm_snv(self, target, pos, allele, min_from_end=8, seq=False, max_depth=2000, max_match=1000):  
	"""Finds reads that show SNV allele at given position"""
	#invalid allele
        if not allele or len(allele) < 1:
	    if seq:
		return None
	    else:
		return 0
		
	num_reads = 0
	reads = []
	# for counting depth
	count = 0
	for read in self.bam.fetch(target, pos, pos + 1):
	    #no duplicate-read allowed
	    if read.is_duplicate:
		continue
	    
	    #filter out chastity-failed reads
	    if read.flag & 512 != 0:
		continue
	    
	    base = self.get_aligned_base(read, pos)
	    
	    if base and base.lower() == allele.lower():
		# determines distance from edge
		from_end = pos - read.pos - 1
		if read.rlen - from_end < from_end:
		    from_end = read.pos + read.rlen - pos
		    	
		# base must be some distance from edge of read
		if from_end >= min_from_end:
		    if seq:
			reads.append(read)
		    else:
			num_reads += 1
			
	    # stop when max_match is reached
	    if num_reads == max_match or len(reads) == max_match:
		break
	    
	    # stop if max_depth is reached
	    count += 1
	    if count == max_depth:
		break
	    
	if not seq:
	    return num_reads
        else:
	    return reads
	
    def get_aligned_base(self, read, target_pos):
	"""Extracts base from read given position from genomic alignment
	Sychronize target and query postions considering indels, soft-clips
	"""
	# cannot return base if cigar string not available
	if read.cigar is None:
	    return None
	
	bases = {}
	query_count = 0
	target_count = read.pos
	for i in range(len(read.cigar)):
	    op, size = read.cigar[i]
	    
	    #start is soft-clipped
	    if i == 0 and op == 4:
		query_count += size
		continue
		
	    #match
	    if op == 0:
		for j in range(size):
		    bases[target_count + j] = read.seq[query_count + j]		    
		target_count += size
		query_count += size
		
	    #insertion
	    elif op == 1:
		query_count += size
	    
	    #deletion
	    elif op == 2:
		target_count += size
		   
	if bases.has_key(target_pos - 1):
	    return bases[target_pos - 1]
	else:
	    return None	    
	
    def confirm_del(self, target, pos, refseq=None, fragment_buffer=1000, seq=False):
	"""Extract all reads containing the deletion
	If the deletion size is the same but has different coordinate in the bam file,
	will compare the sequence after deletion to see if the read contains the same event
	"""
	# deletion size
	size = pos[1] - pos[0] + 1
			
	good_reads = []
	bigger_than_rlen = False
	for read in self.bam.fetch(target, pos[0], pos[1] + 1):
	    # mapq must exceed minimum
	    if int(read.mapq) < int(self.min_mapq):
                continue
	    	    
	    # if deletion size bigger than read length
	    # don't need to look further
	    if size > read.rlen:
		bigger_than_rlen = True
		break
	    
	    # read is not mapped, skipped
	    if read.alen is None:
		continue
	    
	    # make sure read covers deletion 
	    if read.pos + 1 <= pos[0] and read.pos + read.alen >= pos[1]:
		align = self.parse_cigar(read)
		dels = [a for a in align if a[0] == 'del' and a[2] - a[1] + 1 == size]
		
		#if there is a deletion of the same size identified from the cigar string
		if dels:		    		    
		    expected_seq = refseq.GetSequence(target, read.pos + 1, pos[0] - 1) + refseq.GetSequence(target, pos[1] + 1, read.pos + max(read.rlen, read.alen) + size)
		    
		    for d in dels:
			if read.pos+1 > d[0] or read.pos + read.alen < d[1]:
			    continue
			
			#if the start and end of the deletion match, good - simple case
			if d[1] == pos[0] and d[2] == pos[1]:
			    good_reads.append(read)
			    break
			
			# otherwise look for same effective sequence
			read_seq = refseq.GetSequence(target, read.pos + 1, d[1] - 1) + refseq.GetSequence(target, d[2] + 1, read.pos + max(read.rlen, read.alen) + size)
			if read_seq.lower() == expected_seq.lower():
			    good_reads.append(read)
			    break
			    
	# if deletion bigger than read length, try use breakpoint method
	# although if deletion is same size of average fragment length, this won't work
	if bigger_than_rlen:
	    genome_buffer = 2 * fragment_buffer
	    regions = [[target, pos[0] - genome_buffer, pos[0]], [target, pos[1], pos[1] + genome_buffer]]
	    breakpoint = [[target, pos[0]], [target, pos[1]]]	                  
	    good_reads = self.confirm_breakpoint_genomic(regions, breakpoint=breakpoint, seq=True)
		
	if not seq:
	    return len(good_reads)
	else:
	    return good_reads
		
    def parse_cigar(self, read):
	"""Reports cigar string alignment in genomic coordinate
	For use in confirm_ins() and confirm_del()
	Should not be necessary in future versions
	"""
	# list of ('match|ins|del|softclipped', start, end)
	align = []
	
	end = read.pos
	for i in range(len(read.cigar)):				    
	    if read.cigar[i][0] == 0:
		align_type = 'match'		
		start = end + 1
		end = start + read.cigar[i][1] - 1
		align.append([align_type, start, end])
		
	    elif read.cigar[i][0] == 1:
		align_type = 'ins'		
		align.append([align_type, end, read.cigar[i][1]])
		
	    elif read.cigar[i][0] == 2:
		align_type = 'del'
		start = end + 1
		end = start + read.cigar[i][1] - 1	
		align.append([align_type, start, end])
		
	    elif read.cigar[i][0] == 4:
		align_type = 'softclipped'
		start = end + 1
		end = start + read.cigar[i][1] - 1	
		align.append([align_type, start, end])
		
	return align
        
    def confirm_ins(self, target, pos, allele, refseq=None, fragment_buffer=1000, big_cutoff=50, seq=False):
	"""Extracts all reads containing the given insertion
	If the insertion size is the same but has different coordinate in the bam file,
	will compare the sequence after insertion to see if the read contains the same event
	If insertion is bigger than read length, it will not return read objects
	and hence should not use seq=True
	"""
	# size of insertion
        size = len(allele)
			
	# list of reads capturing insertions (smaller than read length)
	good_reads = []
	# number of reads inspected
	num_reads_extracted = 0
	# flag if allele bigger than 1 read
	bigger_than_rlen = False
	
	for read in self.bam.fetch(target, pos, pos + 1):
	    # skip unmapped read
	    if read.alen is None:
		continue
	    
	    num_reads_extracted += 1
	    
	    # if size of allele > read length of first read
	    # stop immediately
	    if len(allele) > read.rlen:
		bigger_than_rlen = True
		break
	    
	    # skip if mapq of read less than minimum
	    if int(read.mapq) < int(self.min_mapq):
                continue
	    
	    # make sure read spans the insertion point (pos)
            if read.pos + 1 <= pos and read.pos + read.rlen >= pos:
		align = self.parse_cigar(read)
		# find all insertions within read that is the same size of the allele
		ins = [a for a in align if a[0] == 'ins' and a[2] == size]
		
		#if there is an insertion of the same size identified from the cigar string
		if ins:
		    # the effective sequence after insertion
		    expected_seq = refseq.GetSequence(target, read.pos + 1, pos) + allele + refseq.GetSequence(target, pos + 1, read.pos + max(read.rlen, read.alen) - size)
		    
		    #go through all insertions to check if the effective sequence is the same
		    for i in ins:
			# skip insertions lying outside the alignments - necessary?
			if read.pos + 1 > i[0] or read.pos + read.alen < i[1]:
			    continue
			
			read_ins_start = i[1] - read.pos
			read_allele = read.seq[read_ins_start:read_ins_start + i[2]]
			    
			#if position is the same and allele is the same, keep - simple case
			if i[1] == pos and read_allele.lower() == allele.lower():
			    good_reads.append(read)
			    break
			
			# check if effective sequences are the same
			read_seq = refseq.GetSequence(target, read.pos + 1, i[1]) + read_allele + refseq.GetSequence(target, i[1]+1, read.pos + max(read.rlen, read.alen) - size)			    
			if read_seq.lower() == expected_seq.lower():
			    good_reads.append(read)
			    break
	
	# insertion bigger than read length - won't return actual reads
	if bigger_than_rlen or (num_reads_extracted == 0 and len(allele) > big_cutoff):
	    genome_buffer = 2 * fragment_buffer
	             
	    # find read pairs that have 1 mate anchored up/downstream and the other
	    # mate sequence is a subsequence of the insertion sequence
	    num_pairs_novel_upstream = self.confirm_novel((target, max(1, pos - genome_buffer), pos), allele, breakpoint=pos)
	    num_pairs_novel_downstream = self.confirm_novel((target, pos + 1, pos + 1 + genome_buffer), allele, breakpoint=pos + 1)
	    # find reads where clipped sequence is part of insertion sequence
	    num_clipped_reads = self.confirm_clipped_reads((target, pos - genome_buffer, pos + genome_buffer), breakpoint=pos, novel_seq=allele)
	    
	    num_reads = num_clipped_reads + min(num_pairs_novel_upstream, num_pairs_novel_downstream) 
	# insertions small than read length
	else:
	    num_reads = len(good_reads)
	       
	if seq:
	    return good_reads
	else:
	    return num_reads
	
    def confirm_novel(self, region, novel_seq, breakpoint, min_olap=10, max_unmapped=1000):
	"""Given novel sequence, find read-pairs that have one unmapped mate that overlaps (with minimum overlap) with novel sequence
	Note of maximum number of pairs with unmapped mates examined (1000)
	Reports the number of pairs that have their unmapped mates entirely subsumed
	or overlap the novel sequence
	"""
	# capture reads with unmapped mates
	unmapped = {}
	# for not going over maximum
	count = 0
	for read in self.bam.fetch(region[0], region[1], region[2]):
	    if read.is_paired and read.mate_is_unmapped:
		count += 1
		unmapped[read.qname] = read
		
		if count == max_unmapped:
		    break
		
	found_mates = 0
	if unmapped:
	    # this assumes unmapped mates is put under same location of mapped mates
	    # store sequences of unmapped mates
	    mate_seq = self.get_unmapped_mate_seq(region[0], unmapped.values())
	    
	    for read in unmapped.values():
		# make anchored read is pointing towards breakpoint
		if (breakpoint == region[1] and not read.is_reverse) or \
		   (breakpoint == region[2] and read.is_reverse):
		    continue		
		
		if mate_seq.has_key(read.qname):
		    if read.pos > region[1]:
			from_end = read.pos - region[1] + 1
		    else:
			from_end = region[1] - read.pos + 1
					    
		    m = re.search(mate_seq[read.qname], novel_seq, re.IGNORECASE)
		    if m:
			found_mates += 1
			
		    else:
			m = re.search(reverse_complement(mate_seq[read.qname]), novel_seq, re.IGNORECASE)
			if m:
			    found_mates += 1
			    
			# if unmapped mate sequence (with or without reverse complement)
			# is not entirely embedded in novel sequence
			# check if it overlaps with novel sequence at the edges
			# if it overlaps at least 10 bases, keep it
			else:
			    olap = seq_overlap(mate_seq[read.qname], novel_seq)
			    if olap > min_olap:
				found_mates += 1
				
			    else:
				olap = seq_overlap(reverse_complement(mate_seq[read.qname]), novel_seq)
				if olap > min_olap:
				    found_mates += 1
				    			    
	return found_mates
	    
    def confirm_clipped_reads(self, region, breakpoint, novel_seq, min_clipped=5):
	"""Finds reads that have clipped sequences at the breakpoint and 
	that overlap with the novel sequence
	Minimum of 5 clipped sequence match is considered positive
	Returns the number of positive reads
	"""
	num_clipped_reads = 0	
	for read in self.bam.fetch(region[0], region[1], region[2]):
	    clipped_seq = self.get_clipped_seq(read, breakpoint)
	    
	    if clipped_seq and len(clipped_seq[1]) > min_clipped:
		clipped_seq_matched = False
		
		# whether the clipped sequence is at 'start' or 'end'
		# dictates the order of the overlap
		if clipped_seq[0] == 'start':
		    olap = seq_overlap(novel_seq, read.seq)		    
		    if olap > len(clipped_seq):
			clipped_seq_matched = True
		    
		else:
		    olap = seq_overlap(read.seq, novel_seq)		    
		    if olap > len(clipped_seq):
			clipped_seq_matched = True
			
		if clipped_seq_matched:
		    num_clipped_reads += 1
		
	return num_clipped_reads
    
    def get_unmapped_mate_seq(self, target, reads, max_reads=1000):
	"""Given a list of reads that don't have mates mapped, extract the sequences of their unmapped mates
	Must NOT be run within another loop of iteration of read extraction
	Used by confirm_novel()
	"""
	seq = {}
	for unmapped_read in reads:
	    for read in self.bam.fetch(target, max(1, unmapped_read.pos - unmapped_read.rlen), unmapped_read.pos + unmapped_read.rlen):
		if read.qname == unmapped_read.qname and read.is_read1 != unmapped_read.is_read1:
		    seq[read.qname] = read.seq
        
	return seq
        
    def get_clipped_seq(self, read, breakpoint):
	"""Extracts clipped sequence and its position given a read
	It makes sure 'breakpoint' lies within the aligned portion
	of the read
	Used by confirm_clipped_read()
	"""
	if read.cigar and (read.cigar[0][0] == 4 or read.cigar[-1][0] == 4):
	    if read.cigar[0][0] == 4:
		start = read.pos - read.cigar[0][1]
	    else:
		start = read.pos
		
	    end = start + read.alen
	    
	    #make sure soft-clipped reads overlap with breakpoint
	    if breakpoint >= start and breakpoint <= end:
		if read.cigar[0][0] == 4:
		    return ('start', read.seq[:read.cigar[0][1]])
		else:
		    return ('end', read.seq[-1 * read.cigar[-1][1]:])
	    
	return None
    
    def is_perfect_align(self, read, ctg_len=None):
	"""Checks if read-to-contig alignment is 'perfect' - no indels"""
	# if cigar string only has 1 part and it's a match
	if len(read.cigar) == 1 and read.cigar[0][0] == 0:
		return True
	    
	# if read is clipped at end or start(2), or both(3)
	elif len(read.cigar) >= 2 and len(read.cigar) <= 3:
	    clipped_start = clipped_end = False
	    if ctg_len is None:
		ctg_len = self.bam.lengths[read.tid]
	    
	    # clipped at start
	    if (read.cigar[0][0] == 4 or read.cigar[0][0] == 5) and read.pos == 0:
		if len(read.cigar) == 2 and read.cigar[1][0] == 0:
		    return True
		clipped_start = True
		
	    # clipped at end
	    if (read.cigar[-1][0] == 4 or read.cigar[-1][0] == 5) and read.pos + read.alen == ctg_len:
		if len(read.cigar) == 2 and read.cigar[0][0] == 0:
		    return True
		clipped_end = True
		
	    # clipped at both start and end, and middle part is a match
	    if clipped_start and clipped_end and len(read.cigar) == 3 and read.cigar[1][0] == 0:
		return True
	    
	return False
		
    def confirm_breakpoint_contig(self, target, pos, only_uniq_frags=False, seq=False, max_depth=None, ctg_len=None):
	"""Finds spanning reads for fusion events
	Returns read objects
	"""
        num_reads = 0
	reads = []

	# coordinate must start from 1
        pos[0] = max(pos[0], 1)
	
	# check if region is suspiciously too highly covered (low-complexity repeat contig)
	too_deep = False
	if max_depth is not None:
	    for pilecolumn in self.bam.pileup(target, pos[0], pos[1] + 1):
		if pilecolumn.pos >= pos[0] and pilecolumn.pos <= pos[1] and pilecolumn.n > max_depth:
		    too_deep = True
		    sys.stdout.write('%s:%s-%s too deeply covered %d (limit:%d)\n' % (target, pos[0], pos[1], pilecolumn.n, max_depth))
		    break
	if too_deep:
	    return reads
		
	last_read = None
        for read in self.bam.fetch(target, pos[0], pos[1] + 1):
	    # don't need to check if current read is same as last read
	    # in terms of pos and cigar - for performance reason
	    if last_read is not None and read.pos == last_read.pos and read.cigar == last_read.cigar:
		if reads and reads[-1].pos == last_read.pos and reads[-1].cigar == last_read.cigar:
		    reads.append(read)
		continue
	    
	    # skip if mapq lower than minimum
            if int(read.mapq) < int(self.min_mapq):
		continue
	    
	    # skip if read is unmapped
	    if read.is_unmapped:
                continue
	    	    
	    # skip if read is not perfectly aligned
	    if not self.is_perfect_align(read, ctg_len):
		continue
	    
	    # skip if read and mate are mapped to same position
	    if not read.mate_is_unmapped and read.rnext == read.tid and read.pnext == read.pos:
		continue
	    
	    # keep if breakpoint region is subsumed in read alignment region
	    if subsume(pos, [read.pos + 1, read.pos + read.alen]):		
		reads.append(read)
		    
	    last_read = read
				
	return reads    
		
    def check_genome_region(self, region):
	"""Checks validity of genomic region given""" 
	# get rid of 'chr' if bam file's genome doesn't have chr
        if not 'chr' in self.bam.references[0] and 'chr' in region[0]:
            region[0] = region[0].replace('chr', '')
	    
	# make sure start > 0
	region[1] = max(1, region[1])
	
        if not region[0] in self.bam.references:
	    return False
	
        # don't do whole chromosome comparisons
        if int(region[1]) == 0 and int(region[2]) == 0:
	    return False
	
	return True
    	
    def confirm_breakpoint_genomic(self, regions, max_pairs=None, breakpoint=None, seq=False, breakpoint_buffer=0, focal=False, no_proper=False):
	"""Finds read pairs supporting fusion events
	Returns both the number and the actual read pairs
	"""
	# check validity of genomic region
	if not self.check_genome_region(regions[0]) or not self.check_genome_region(regions[1]):
	    return 0, None
			
        num_reads = 0
        region1 = regions[0]
        region2 = regions[1]
	
	# find reads pointing to breakpoint in first region
	reads1 = self.find_mates_in_genome(region1, breakpoint[0], breakpoint_buffer=breakpoint_buffer, outside_breakpoint=True, maximum=max_pairs, no_proper=no_proper)	    
	reads1_dict = dict((r.qname, r) for r in reads1)
	
	# find reads pointing to breakpoint in second region
	reads2 = self.find_mates_in_genome(region2, breakpoint[1], breakpoint_buffer=breakpoint_buffer, outside_breakpoint=True, maximum=max_pairs, no_proper=no_proper)
	reads2_dict = dict((r.qname, r) for r in reads2)
		
	# find pairs
	pairs = {}
	# order of magnitude of distance between breakpoints
	breakpoint_distance_order = round(math.log10(max(1, abs(breakpoint[0][1] - breakpoint[1][1]))))
	# order of magnitude of distance between breakpoints allowed
	max_distance_order = 0
	if focal:
	    max_distance_order = 1
		
	for read in reads1:
	    if reads2_dict.has_key(read.qname):	
		# ensure that a mate does not come between a read and the breakpoint
		if not read.is_reverse:
		    if reads2_dict[read.qname].pos > read.pos and reads2_dict[read.qname].pos < breakpoint[0][1]:
			continue
		else:
		    if reads2_dict[read.qname].pos > breakpoint[0][1] and reads2_dict[read.qname].pos < read.pos:
			continue
		if not reads2_dict[read.qname].is_reverse:
		    if read.pos > reads2_dict[read.qname].pos and read.pos < breakpoint[1][1]:
			continue
		else:
		    if read.pos > breakpoint[1][1] and read.pos < reads2_dict[read.qname].pos:
			continue
		    
		# mates cannot not have same pos
		if read.pos == reads2_dict[read.qname].pos:
		    continue
		
		# ensure pair support are at least the same order of magnitude apart as are the breakpoints
		mates_distance_order = round(math.log10(max(1, abs(read.pos - reads2_dict[read.qname].pos))))
		if abs(mates_distance_order - breakpoint_distance_order) > max_distance_order:
		    continue
		
		pairs[read.qname] = read, reads2_dict[read.qname]
		# no need to look for more if maximum is give and reached 
		if max_pairs is not None and len(pairs.keys()) >= max_pairs:
		    break
		
	# calculate stacking - not used 
	stacks1 = {}
	stacks2 = {}
	for name, reads in pairs.iteritems():
	    if not stacks1.has_key(reads[0].pos):
		stacks1[reads[0].pos] = 0
	    stacks1[reads[0].pos] += 1
	    
	    if not stacks2.has_key(reads[1].pos):
		stacks2[reads[1].pos] = 0
	    stacks2[reads[1].pos] += 1
		
	num_reads = len(pairs.keys())
			    
	if not seq:
	    return num_reads, None
	else:
	    return num_reads, pairs.values()
	
    def find_mates_in_genome(self, region, breakpoint, 
                             breakpoint_buffer=0, missing_mates=None, outside_breakpoint=True, 
                             maximum=None, no_duplicates=True, no_chastity=True, no_proper=True):
	"""Find mates of fusion events in genome alignments"""
	if not self.check_genome_region(region):
	    return []
	
	# make sure read name doesn't have '_' - GSC hack
	if missing_mates:
	    missing_mates = dict((r.replace('_', ':'), missing_mates[r]) for r in missing_mates)
			
	reads = {}
	for read in self.bam.fetch(region[0], region[1], region[2]):  
	    # if read is not mapped, skip
	    if read.alen is None:
		continue
	    
	    # skip PCR duplicate
	    if not missing_mates and no_duplicates and read.is_duplicate:
		continue
	    
	    # skip read below minimum mapq
	    if int(read.mapq) < int(self.min_mapq):
		continue
	    
	    # filter out chastity-failed reads
	    if not missing_mates and no_chastity and read.flag & 512 != 0:
		continue
	    
	    # if want to use 'proper_pair' flag, then make sure read is not 'proper_pair'
	    if no_proper and read.is_proper_pair:
		continue
	    
	    # make sure read is in between interval
	    if not missing_mates and not subsume([read.pos + 1, read.pos + read.rlen], [region[1] - read.alen, region[2] + read.alen]):
		continue
	    	    
	    # if given list of missing mates, make sure read is one of them
	    if missing_mates and not self.is_mate(read, missing_mates):
		continue
	    	
	    # for checking if mate is point towards breakpoint
	    pointing_correctly = False

	    # for checking if mate lies completely on one side of breakpoint
	    completely_on_one_side = False
	    # allows read to overlap breakpoint by one read length if missing mates given
	    if missing_mates:
		breakpoint_buffer = read.rlen
	    
	    # checks if read is pointing towards breakpoint and lies completely on one side
	    # sense
	    if not read.is_reverse:
		if read.pos < breakpoint[1]:
		    pointing_correctly = True		
		    if read.pos + read.alen - 1 <= breakpoint[1] + breakpoint_buffer:
			completely_on_one_side = True						
	    # anti-sense
	    elif read.pos + read.alen - 1 > breakpoint[1]:
		pointing_correctly = True		
		if read.pos >= breakpoint[1] - breakpoint_buffer:
		    completely_on_one_side = True
		    
		    	
	    if pointing_correctly:		
		if outside_breakpoint:
		    if completely_on_one_side:
			# if 2 mates are both candidates, choose the one closer to breakpoint
			if not reads.has_key(read.qname) or abs(read.pos - breakpoint[1]) < abs(reads[read.qname].pos - breakpoint[1]):
			    reads[read.qname] = read
			    
		else:
		    if not reads.has_key(read.qname) or abs(read.pos - breakpoint[1]) < abs(reads[read.qname].pos - breakpoint[1]):
			reads[read.qname] = read
	    
	    # check if maximum is reached if maximum is imposed
	    if maximum is not None and len(reads) == maximum:
		break
	    	    
	return reads.values()
   			
    def is_mate(self, read, mates):
	"""Given a read and a dictionary of mates (key:read.qname, value: read),
	see if read is mate of one of the mates
	Assume mates have the same name
	"""
	result = False
	
	# make sure read name doesn't have '_' - GSC hack
	key = read.qname.replace('_', ':')
	if mates.has_key(key):
	    for mate in mates[key]:
		if (read.is_read1 and not mate.is_read1) or (not read.is_read1 and mate.is_read1):
		    result = True
		elif read.seq != mate.seq:
		    result = True
		if result:
		    break
		
	return result
    	
    def coverage_simple(self, target, pos=[], total=False, seq=False):
        """Returns average coverage or reads spanning region (total=True)"""  
        coverage_profile = {}
        for i in range(pos[0], pos[1] + 1):
            coverage_profile[i] = 0

        reads = []
        num_reads = 0
        for read in self.bam.fetch(target, pos[0], pos[1]):
            read_end = read.pos + read.rlen - 1
            if read.pos + 1 <= pos[1] and read_end + 1 >= pos[0]:
                num_reads += 1
                
                if seq:
                    reads.append(read)
                    
                for i in range(max(read.pos, pos[0]-1), min(read_end, pos[1]-1)+1):
                    coverage_profile[i + 1] += 1

        if not total:
            return "%.1f" % (sum(coverage_profile.values()) / len(coverage_profile.keys()))
        elif seq:
            return reads
        else:
            return num_reads
                    
    @classmethod
    def base_qual_to_int(cls, qual_char, offset):
	"""Converts base quality base to integer value"""
	if offset == 64 or offset == 33:
	    return ord(qual_char) - offset
	
	return None

    @classmethod
    def int_to_base_qual(cls, qual_int, offset):
	"""Converts integer to base quality base"""
	if offset == 64 or offset == 33:
	    return chr(qual_int + offset)
	
	return None
    