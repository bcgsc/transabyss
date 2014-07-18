"""
This module stores SNV/indel event data and methods

Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
import sys
from utilities.bam import BAM
from utilities import tools

class SNV:
    """Individual indel or SNV event identified by snv_caller.py"""
    def __init__(self, method=None, snv_type=None, ref=None, ref_start=None, ref_end=None, ref_seq=None, 
                 query_strand=None, var=None, var_start=None, var_end=None, var_seq=None, var_len=None):
        self.method = method
        self.snv_type = snv_type
        self.var = var
        self.var_start = var_start
        self.var_end = var_end
        self.var_seq = var_seq
        self.ref = ref
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.ref_seq = ref_seq
        self.var_len = var_len
        
        self.snv_len = 0
        if snv_type in ('ins', 'dup', 'ITD', 'PTD', 'inv') and var_seq:
            self.snv_len = len(var_seq)
            
        elif ref_seq:
            self.snv_len = len(ref_seq)

        self.from_end = None
        self.nreads_contig = 'na'
        self.nreads_genome = 'na'
        self.nreads_event = 'na'
        self.gene = 'na'
        self.query_strand = query_strand
        
        self.expansion = 0

        self.confirm_contig_region = []
        if self.var_start and self.var_end:
            self.confirm_contig_region = [int(self.var_start), int(self.var_end)]
            
        # filtering
        self.enough_coverage = False
        self.exon = False
        self.nonsynon = False
        self.too_close_to_end = False
        self.at_least_1_read_opposite = False

        # annotations
        self.within_simple_repeats = '-'
        self.repeatmasker = '-'
        self.within_segdup = '-'
        self.dbsnp = '-'
        
        # for filtering snvs from bad contig region
        self.artefact = False
                
    def gff(self, event_type):
        """Output event in gff (not maintained anymore)"""

        if re.match('^(chr|scaffold)', self.ref, re.IGNORECASE):
            chrom = self.ref
        else:
            chrom = 'chr' + self.ref
        var_shortened = self.var.split(" ")[0]

        if self.snv_type.upper() == "SNV":
            label_last = self.ref_seq.upper() + "->" + self.var_seq.upper()
        else:
            label_last = self.snv_type.upper()
        label = ":".join([var_shortened, chrom, str(self.ref_start), label_last])

        out = "\t".join([chrom, self.method, event_type, str(self.ref_start), str(self.ref_end), ".", '+', ".", label])
        return out + "\n"

    def tab(self, debug=False):
        """Outputs event in tabular format"""      
        chrom = self.ref
        var_shortened = self.var.split(" ")[0]

        ref = alt = 'na'
        if self.snv_type == "snv" or self.snv_type == 'inv':
            ref = self.ref_seq
            alt = self.var_seq
        elif self.snv_type == "del":
            if len(self.ref_seq) > 10:
                ref = self.ref_seq[:5] + '...' + self.ref_seq[-5:]
            else:
                ref = self.ref_seq
        elif self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
            alt = self.var_seq
            
        if not self.gene:
            self.gene = 'na'
            
        if not self.at_least_1_read_opposite:
            at_least_1_read_opposite = 'false'
        else:
            at_least_1_read_opposite = 'true'
                    
        fields = [self.snv_type, chrom, str(self.ref_start), str(self.ref_end), var_shortened, str(self.var_len), str(self.var_start), str(self.var_end), str(self.snv_len), ref, alt, str(self.nreads_event), str(self.nreads_contig), str(self.nreads_genome), self.gene, str(self.expansion), self.query_strand, str(self.from_end), "%s-%s" % (self.confirm_contig_region[0], self.confirm_contig_region[1]), self.within_simple_repeats, self.repeatmasker, self.within_segdup, at_least_1_read_opposite, self.dbsnp]
        
        if debug:
            fields.append("enough_reads:" + str(self.enough_coverage))
            fields.append("too_close_to_end:" + str(self.too_close_to_end))
            fields.append("at_least_1_read_opposite:" + str(self.at_least_1_read_opposite))
            fields.append("exon:" + str(self.exon))
            fields.append("nonsynon:" + str(self.nonsynon))
            
        return "\t".join(fields) + "\n"

    def coord(self):
        """Outputs event coordinate in UCSC format"""
        return "%s:%s-%s" % (self.ref, self.ref_start, self.ref_end)

    def genome_read_support(self, bam, refseq, from_end=8):
        """Finds support reads from reads-to-genome BAM"""
        chrom = self.ref 
	if not chrom in bam.bam.references:
	    if 'chr' in chrom and chrom.replace('chr', '') in bam.bam.references:
		chrom = chrom.replace('chr', '') 
	    elif not 'chr' in chrom and 'chr' + chrom in bam.bam.references:
		chrom = 'chr' + chrom
        
        if self.snv_type in ('ins', 'dup', 'ITD', 'PTD') or self.snv_type == 'del':  
            if self.snv_type == 'del':
                self.nreads_genome = bam.confirm(chrom, [int(self.ref_start),int(self.ref_end)], feature=self.snv_type, refseq=refseq)
            else:
                if self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
                    pos = int(self.ref_start)
                else:
                    pos = int(self.ref_start) - 1
                            
                self.nreads_genome = bam.confirm(chrom, [pos], feature='ins', allele=self.var_seq, refseq=refseq)
                
        elif self.snv_type == 'snv':
            self.nreads_genome = bam.confirm(chrom, [int(self.ref_start)], feature=self.snv_type, allele=self.var_seq)
            
        elif self.snv_type == 'inv':
            self.nreads_genome = bam.confirm_perfect(chrom, [int(self.ref_start),int(self.ref_end)], self.var_seq, '+', min_from_end=from_end)
            
        return self.nreads_genome

    def contig_read_support(self, bam, lib=None, from_end=8, get_reads=False):
        """Finds support reads from reads-to-contig BAM"""
        reads = []
        if self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
            if int(self.var_end)-int(self.var_start)+1 != int(self.snv_len):
                self.nreads_contig = 0
            else:
                if get_reads:
                    reads = bam.confirm_perfect(self.var, self.confirm_contig_region, self.var_seq, self.query_strand, min_from_end=from_end, expansion=self.expansion, seq=True)
                    self.nreads_contig = len(reads)
                else:
                    self.nreads_contig = bam.confirm_perfect(self.var, self.confirm_contig_region, self.var_seq, self.query_strand, min_from_end=from_end, expansion=self.expansion)
                    
        elif self.snv_type == 'snv' or self.snv_type == 'inv':
            if int(self.var_end)-int(self.var_start)+1 != int(self.snv_len):
                self.nreads_contig = 0
            else:        
                if get_reads:
                    reads = bam.confirm_perfect(self.var, self.confirm_contig_region, self.var_seq, self.query_strand, min_from_end=from_end, seq=True)
                    self.nreads_contig = len(reads)
                else:
                    self.nreads_contig = bam.confirm_perfect(self.var, self.confirm_contig_region, self.var_seq, self.query_strand, min_from_end=from_end)
                
        elif self.snv_type == 'del':
            if get_reads:
                reads = bam.confirm(self.var, [int(self.confirm_contig_region[0])-from_end+1,int(self.confirm_contig_region[1])+from_end-1], feature='breakpoint', seq=True, ctg_len=int(self.var_len))
                self.nreads_contig = len(reads)
            else:
                self.nreads_contig = bam.confirm(self.var, [int(self.confirm_contig_region[0])-from_end+1, int(self.confirm_contig_region[1])+from_end-1], feature='breakpoint')

        if get_reads:
            return reads

    def upshift(self, refseq): 
        """Shifts event start coordinate upstream for repeat-involved event"""
        if self.snv_type in ('ins', 'dup', 'ITD', 'PTD') or self.snv_type == 'del':
            if self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
                if tools.is_homopolymer(self.var_seq):
                    size = 1
                    seq = self.var_seq[0]
                else:
                    size = len(self.var_seq)
                    seq = self.var_seq
            else:
                if tools.is_homopolymer(self.ref_seq):
                    size = 1
                    seq = self.ref_seq[0]
                else:
                    size = len(self.ref_seq)
                    seq = self.ref_seq
                    
            start = int(self.ref_start) - size
            # skip if 0-size event
            if size == 0:
                sys.stderr.write("error in upshift size 0 contig:%s %s\n" % (self.var, self.ref_seq))
                return
            
            # continues checking upstream sequence to see if it's repeat of sequence in question
            while start > 1:
                upstream = refseq.GetSequence(self.ref, start + 1, start + size)
                if seq.upper() != upstream.upper():
                    break
                start = start - size
                
            # changed reference start coordinate (and end coordinate too if deletion)
            if start + size < self.ref_start:
                if self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
                    sys.stderr.write("shifted %s %s %s %s to %d\n" % (self.var, self.snv_type, self.ref, self.ref_start, start+size))
                    self.ref_start = self.ref_end = start + size
                else:
                    sys.stderr.write("shifted %s %s %s %s to %d\n" % (self.var, self.snv_type, self.ref, self.ref_start, start+size+1))
                    self.ref_start = start + size + 1
                    self.ref_end = self.ref_start + len(self.ref_seq) - 1

    def expand_contig_region(self, contig_sequence, query_strand):
        """Expand read-support checking region if repeats are involved"""
        if not self.snv_type in ('ins', 'dup', 'ITD', 'PTD', 'del'):
            return None
        
        # skip if deleted/inserted sequence is longer than contig sequence
        if self.snv_type == 'del' and len(self.ref_seq) > len(contig_sequence):
            return None
        if self.snv_type in ('ins', 'dup', 'ITD', 'PTD') and len(self.var_seq) > len(contig_sequence):
            return None

        if self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
            seq = self.var_seq[:]
        else:
            seq = self.ref_seq[:]
        
        if len(seq) == 0:
            return None
            
        if tools.is_homopolymer(seq) or len(seq) == 1:
            homo = True
        else:
            homo = False

        # keep a record of previous value for reporting expansion
        region_before = self.confirm_contig_region[:]
        
        # arbitrary big number
        limit = 100000
        
        # forward
        expand = 0
        for i in range(limit):
            if homo:
                changed_base = seq[0].upper()
            else:
                changed_base = seq[i % len(seq)].upper()
                
            downstream_base = None
            if self.snv_type == 'del':
                if query_strand == '+':
                    if int(self.var_end) + i < len(contig_sequence) and int(self.var_end) + i >= 0:
                        downstream_base = contig_sequence[int(self.var_end) + i].upper()
                else:
                    if int(self.var_end) - 2 - i >= 0 and int(self.var_end) - 2 - i < len(contig_sequence): 
                        downstream_base = tools.reverse_complement(contig_sequence[int(self.var_end) - 2 - i]).upper()

            elif self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
                if query_strand == '+':
                    if int(self.var_end) + i < len(contig_sequence) and int(self.var_end) + i >= 0:
                        downstream_base = contig_sequence[int(self.var_end) + i].upper()
                else:
                    if int(self.var_start) - i - 2 >= 0 and int(self.var_start) - i - 2 < len(contig_sequence):
                        downstream_base = tools.reverse_complement(contig_sequence[int(self.var_start) - i - 2]).upper()

            if changed_base == downstream_base:
                expand += 1   
            else:
                break
        
        multiples = expand/len(seq)
        if multiples > 0:
            if query_strand == '+':
                self.confirm_contig_region[1] += multiples * self.snv_len
            else:
                self.confirm_contig_region[0] -= multiples * self.snv_len

        # reverse
        seq = seq[::-1]
        expand = 0

        for i in range(limit):
            if homo:
                changed_base = seq[0].upper()
            else:
                changed_base = seq[i%len(seq)].upper()
                
            upstream_base = None
            if self.snv_type == 'del':
                if query_strand == '+':
                    if int(self.var_start) - i - 1 >= 0 and int(self.var_start) - i - 1 < len(contig_sequence):
                        upstream_base = contig_sequence[int(self.var_start) - i - 1].upper()
                else:
                    if int(self.var_start) + i - 1 < len(contig_sequence) and int(self.var_start) + i - 1 >= 0:
                        upstream_base = tools.reverse_complement(contig_sequence[int(self.var_start)+i-1]).upper()
                        
            elif self.snv_type in ('ins', 'dup', 'ITD', 'PTD'):
                if query_strand == '+':
                    if int(self.var_start) - i - 2 >= 0 and int(self.var_start) - i - 2 < len(contig_sequence):
                        upstream_base = contig_sequence[int(self.var_start) - i - 2].upper()
                else:
                    if int(self.var_end) + i < len(contig_sequence) and int(self.var_end) + i >= 0:
                        upstream_base = tools.reverse_complement(contig_sequence[int(self.var_end) + i]).upper()

            if changed_base == upstream_base:
                expand += 1   
            else:
                break

        multiples = expand/len(seq)
        if multiples > 0:
            if query_strand == '+':
                self.confirm_contig_region[0] -= multiples * self.snv_len
            else:
                self.confirm_contig_region[1] += multiples * self.snv_len
                
        expanded_sequence = contig_sequence[self.confirm_contig_region[0]-1:self.confirm_contig_region[1]]
        
        # coordinate given in 1-based
        if region_before[0] != self.confirm_contig_region[0] or region_before[1] != self.confirm_contig_region[1]:
            self.expansion = (self.confirm_contig_region[1] - self.confirm_contig_region[0] + 1) / self.snv_len
            sys.stderr.write("expand confirm contig region %s %s -> %s %s %s %sx\n" % (self.var, region_before, self.confirm_contig_region, expanded_sequence, len(expanded_sequence), self.expansion))
            
