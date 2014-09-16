"""
This module finds indels, inversion, duplications, and SNVs through un-split alignments

Author: Readman Chiu rchiu@bcgsc.ca
"""
__version__ = '1.5.1'

import os
import sys
import re
import glob
from optparse import OptionParser, OptionGroup
from snv import SNV
from feature import FeatureFinder
from utilities import tools, intspan
from utilities.assembly import Assembly
from utilities.track import Track
from utilities.bam import BAM
from utilities.align_parsers import psl, sam
from utilities.overlap_coord import OverlapCoord
from annotations import repeat, dbsnp
from bubble import BubbleMapper
import utilities.cfg

PACKAGE_DIR = "/".join(os.path.abspath(sys.argv[0]).split("/")[:-2])

class SNVCaller:
    """Calls indels and SNVs from single, non-split contig alignment""" 
    event_types = ('snv', 'ins', 'del', 'inv', 'dup', 'ITD')

    output_headers = ['id',
                      'type',
                      'chr',
                      'chr_start',
                      'chr_end',
                      'ctg',
                      'ctg_len',
                      'ctg_start',
                      'ctg_end',
                      'len',
                      'ref',
                      'alt',
                      'event_reads',
                      'contig_reads',
                      'genome_reads',
                      'gene',
                      'repeat-length',
                      'ctg_strand',
                      'from_end',
                      'confirm_contig_region',
                      'within_simple_repeats',
                      'repeatmasker',
                      'within_segdup',
                      'at_least_1_read_opposite',
                      'dbsnp'
                      ]

    def __init__(self, genome, get_snv, get_indel, align_file=None, contigs_file=None, min_reads=0, 
                 debug=False, min_from_end=None, bubble_mapping=None, sample_type='transcriptome', gene_model=None, annodir=None, mmcfg=None):
        self.align_file = align_file
        self.contigs_file = contigs_file
        self.genome = genome
        self.min_reads = min_reads
        self.debug = debug
        self.min_from_end = min_from_end
        self.gene_model = gene_model
        self.annodir = annodir
        self.mmcfg = mmcfg

        self.chrom_proper, self.refseq, self.ff, self.repeat_overlaps, self.splice_motif_file = None, None, None, None, None
        if self.genome:
            self.prepare_annotation()

        self.snvs = []
        self.grouped_snvs = {}
        
        self.bubble_mapping = bubble_mapping
        self.sample_type = sample_type
                
        # actions to perform
        self.get_snv = get_snv
        self.get_indel = get_indel
                
    def prepare_annotation(self):
        """Prepares for overlapping annotations"""
        if self.genome:
            # for conversion of chromosome name between FASTA and annotation (hg19)
            self.chrom_proper = tools.ucsc_chroms(self.genome, self.annodir)
            
            splice_motif_file = os.path.join(self.annodir, self.genome, 'splice_motifs.txt')
            if os.path.isfile(splice_motif_file):
                self.splice_motif_file = splice_motif_file
                
            # for constructing cDNA sequences for in/out frame determination
            self.refseq = tools.get_refseq_from_2bit(self.annodir, self.genome)
            
            # for finding genes, exons, etc
            if self.gene_model:
                self.ff = FeatureFinder(self.genome, self.gene_model, self.annodir, self.mmcfg)

            # for getting rid of contigs that mapped entirely with repeats
            self.repeat_overlaps = repeat.prepare_overlap(self.genome, self.annodir)
            
    def extract(self, cutoff=None, no_group=False, match_percent=None, identity=None, no_segdup=False): 
        """Wrapper function for identifying indels and SNVs from non-split alignments"""
        splice_motifs = tools.get_splice_motifs(self.splice_motif_file)
        filters = {'unique':True, 'bestn':1, 'match':match_percent, 'identity':identity}
                
        # extracts alignments
        out_format = os.path.splitext(self.align_file)[1]
        aligns = {
            '.psl': psl.parse,
            '.sam': sam.parse
            }[out_format](self.align_file, filters, splice_motif_file=self.splice_motif_file, refseq=self.refseq, noline=True)
                                            
        # links contig sequence to alignment
        ass = Assembly(None, k=1)
        ass.fasta = self.contigs_file
        contigs = ass.get_contigs(sequence=True)
        contig_dict = dict((contig.num, contig) for contig in contigs)
        for align in aligns:
            if contig_dict.has_key(align.query):
                align.contig = contig_dict[align.query]

        for align in aligns:
            if self.bubble_mapping and not self.bubble_mapping.is_bubble_mapped_to_contig(align.query):
                print "remove bubble", align.query
                continue
                         
            snvs = self.get_snvs(align, splice_motifs=splice_motifs, cutoff=cutoff, no_segdup=no_segdup)

            for snv in snvs:
                snv.var_len = align.query_len
                snv.from_end = min(int(snv.var_start) - int(align.qstart), int(align.qend) - int(snv.var_end))

                # identifies repeat units
                snv.upshift(self.refseq)
                snv.expand_contig_region(align.contig.sequence, align.query_strand)
                target = snv.ref
                
                # re-labels 'ins' as 'dup' if expansion >=2 and length > 3
                if snv.snv_type == 'ins' and snv.expansion >= 2 and snv.snv_len > 3:
                    snv.snv_type = 'dup'
                    snv.ref_start += 1
                    snv.ref_end = snv.ref_start + snv.snv_len - 1
                                                                        
                self.snvs.append(snv)
                                
        # group events
        if not no_group:
            self.grouped_snvs = self.group(self.snvs)
            
    def annotate_genes(self):
        """Annotates events with genes"""
        if self.ff:
            for groups in self.grouped_snvs.values():
                for snvs in groups.values():
                    gene = self.ff.get_feature(' '.join((tools.proper_chrom(snvs[0].ref, 
                                                                            chrom_proper=self.chrom_proper), 
                                                         str(snvs[0].ref_start), 
                                                         str(snvs[0].ref_end))), 
                                               refseq=self.refseq, 
                                               variant=snvs[0].var_seq, 
                                               change=snvs[0].snv_type, 
                                               chrom=snvs[0].ref)
                
                    for snv in snvs:
                        snv.gene = gene
                        
                    # relabel duplication within exon as ITD
                    if snvs[0].snv_type == 'dup' and 'exon' in snvs[0].gene\
                       and not 'utr' in snvs[0].gene and not 'intron' in snvs[0].gene:
                        for snv in snvs:
                            snv.snv_type = 'ITD'
                            
    def is_within_repeats(self, proper_chrom, span):
        """Determines if given coordinate span overlaps segdups or simple_repeats"""
        overlaps = repeat.find_overlaps({'chrom':proper_chrom, 'start':span[0], 'end':span[1]}, self.repeat_overlaps)
        if overlaps:
            for repeat_type in overlaps.keys():
                if overlaps[repeat_type]:
                    if repeat_type == 'segdup' or repeat_type == 'simple_repeats':
                        return True
                
        return False
                                
    def overlap_repeats(self):        
        """Overlaps event coordinates with repeatmasker, simple repeats"""
        event_groups_by_chr = self.group_by_chr()
        
        for chrom, events in event_groups_by_chr.iteritems():          
            proper_chrom = tools.proper_chrom(chrom, chrom_proper=self.chrom_proper)
            
            for snv_type, snv_groups in events.iteritems():
                print 'processing repeat', snv_type                
                for snvs in snv_groups:
                    overlaps = repeat.find_overlaps({'chrom':proper_chrom, 'start':int(snvs[0].ref_start), 
                                                     'end':int(snvs[0].ref_end)}, self.repeat_overlaps)
                    if overlaps:
                        attrs = {}
                        for repeat_type, types in overlaps.iteritems():
                            if repeat_type == 'simple_repeats':
                                attr = 'within_simple_repeats'
                            elif repeat_type == 'segdup':
                                attr = 'within_segdup'
                            elif repeat_type == 'rmsk':
                                attr = 'repeatmasker'
                        
                            if types:
                                # only report one with shortest name
                                types_sorted = types.keys()
                                types_sorted.sort(lambda x,y: len(x)-len(y))
                                attrs[attr] = types_sorted[0]
                        
                            if attrs:
                                for snv in snvs:
                                    tools.set_attrs(snv, attrs)          
        # clears cache                                        
        for repeat_olap in self.repeat_overlaps.values():
            repeat_olap.finish()
                        
    def group_by_chr(self):
        """Groups events by genomic location for overlapping with annotation such as dbSNP"""
        events_by_chr = {}
        for event_type in self.event_types:
            if not self.grouped_snvs.has_key(event_type):
                continue
            
            for key, snvs in self.grouped_snvs[event_type].iteritems():
                if not events_by_chr.has_key(snvs[0].ref):
                    events_by_chr[snvs[0].ref] = {}
                    events_by_chr[snvs[0].ref] = {event_type:[]}
                    
                if not events_by_chr[snvs[0].ref].has_key(event_type):
                    events_by_chr[snvs[0].ref][event_type] = []
                    
                events_by_chr[snvs[0].ref][event_type].append(snvs)
                
        return events_by_chr
     
    def overlap_dbsnp(self):
        """Overlaps events with dbSNP"""
        event_groups_by_chr = self.group_by_chr()
        
        for chrom, events in event_groups_by_chr.iteritems():          
            proper_chrom = tools.proper_chrom(chrom, chrom_proper=self.chrom_proper)
            snp_overlap = dbsnp.prepare_overlap(self.genome, proper_chrom, self.annodir)
            
            for snv_type, snv_groups in events.iteritems():
                event_type_check = snv_type
                if snv_type in ('dup', 'ITD', 'PTD'):
                    event_type_check = 'ins'
                
                for snvs in snv_groups:
                    start, end = int(snvs[0].ref_start), int(snvs[0].ref_end)
                    if snv_type in ('dup', 'ITD', 'PTD'):
                        start, end = int(snvs[0].ref_start) - 1, int(snvs[0].ref_start) - 1
                        
                    known = dbsnp.find_concordants({'type':event_type_check, 'chrom':proper_chrom, 'start':start, 'end':end, 'allele':snvs[0].var_seq.lower(), 'size':int(snvs[0].snv_len)}, snp_overlap, refseq=self.refseq, target=chrom)                    
                    if known:
                        for snv in snvs:
                            snv.dbsnp = ','.join(known)
                        
            snp_overlap.finish()
            
    def compare_contig(self, snv1, snv2):
        """For sorting events by comparing contig number and coordinates"""
        if snv1.var.isdigit() and snv2.var.isdigit():
            if int(snv1.var) < int(snv2.var):
                return -1
            elif int(snv1.var) > int(snv2.var):
                return 1
            elif int(snv1.var_start) < int(snv2.var_start):
                return -1
            elif int(snv1.var_start) > int(snv2.var_start):
                return 1
            else:
                return 0            
        elif snv1.var.isdigit() and not snv2.var.isdigit():
            return -1
        elif not snv1.var.isdigit() and snv2.var.isdigit():
            return 1        
        else:
            if snv1.var < snv2.var:
                return -1
            elif snv1.var > snv2.var:
                return 1
            elif int(snv1.var_start) < int(snv2.var_start):
                return -1
            elif int(snv1.var_start) > int(snv2.var_start):
                return 1
            else:
                return 0
            
    def add_support(self, lib=None, from_end=None, genome_bamfile=None, contigs_bamfile=None):
        """Calculates read support"""
        genome_bam = None
        contigs_bam = None
        if genome_bamfile and os.path.exists(genome_bamfile):
            genome_bam = BAM(genome_bamfile)
                                
            for event_type in sorted(self.grouped_snvs.keys()):
                for coord_allele in self.grouped_snvs[event_type]:
                    snvs = self.grouped_snvs[event_type][coord_allele]
                    snvs[0].genome_read_support(genome_bam, self.refseq, from_end=from_end)                        
                    for i in range(1, len(snvs)):
                        snvs[i].nreads_genome = snvs[0].nreads_genome
            
        if contigs_bamfile and os.path.exists(contigs_bamfile):
            contigs_bam = BAM(contigs_bamfile, min_mapq=0)
            self.snvs.sort(self.compare_contig)            
            for event in self.snvs:
                event.contig_read_support(contigs_bam, lib=lib, get_reads=True, from_end=from_end)
                        
    def group(self, snvs):
        """Groups events by event type, coordinate, and allele sequence"""
        grouped_snvs = {}
        for snv in snvs:
            if snv.artefact:
                continue
            
            seq = ""
            if snv.snv_type in self.event_types:
                seq = snv.var_seq[:]
            coord = snv.coord()
            snv_key = coord + "-" + seq
            
            if not grouped_snvs.has_key(snv.snv_type):
                grouped_snvs[snv.snv_type] = {snv_key:[snv]}
            else:
                if not grouped_snvs[snv.snv_type].has_key(snv_key):
                    grouped_snvs[snv.snv_type][snv_key] = [snv]
                else:
                    grouped_snvs[snv.snv_type][snv_key].append(snv)
                    
        return grouped_snvs
                
    def add_contig_reads(self, grouped_snvs):
        """Sums contig read supports of members within same group"""
        for type in self.event_types:
            if not grouped_snvs.has_key(type):
                continue
        
            for coord in grouped_snvs[type].keys():
                event_reads = 0
                for snv in grouped_snvs[type][coord]:
                    if snv.nreads_contig != 'na':
                        event_reads += int(snv.nreads_contig)
                    
                for snv in grouped_snvs[type][coord]:
                    snv.nreads_event = event_reads

    def output_groups(self, grouped_snvs, output, debug=False):
        """Outputs grouped events with identifiers indicating groups"""
        # group count
        count1 = 1
        for type in self.event_types:
            if not grouped_snvs.has_key(type):
                continue
        
            coords = sorted(grouped_snvs[type].keys(), key=lambda coord_allele: (grouped_snvs[type][coord_allele][0].ref, int(grouped_snvs[type][coord_allele][0].ref_start)))
            for coord in coords:
                # group-member count
                count2 = 1
                for snv in grouped_snvs[type][coord]:
                    if len(grouped_snvs[type][coord]) == 1:
                        count = count1
                    else:
                        count = "%d.%d" % (count1, count2)

                    output.write(str(count) + "\t" + snv.tab(debug=debug))
                
                    count2 += 1

                count1 += 1

    def report(self, outdir, post_filter=False):
        """Produces different output files"""
        txt_file = outdir + "/events.tsv"
        self.output_txt(self.snvs, txt_file)

        if post_filter:
            self.post_filter(min_from_end=self.min_from_end)
        
            if self.sample_type == 'genome':
                filtered_snvs = [snv for snv in self.snvs if not snv.artefact and snv.enough_coverage and not snv.too_close_to_end and snv.at_least_1_read_opposite ]
            else:
                filtered_snvs = [snv for snv in self.snvs if not snv.artefact and snv.enough_coverage and not snv.too_close_to_end ]
            
            txt_file = outdir + "/events_filtered.tsv"
            self.output_txt(filtered_snvs, txt_file)
        
            exon_snvs = [snv for snv in filtered_snvs if snv.exon and snv.nonsynon]
            txt_file = outdir + "/events_exons.tsv"
            self.output_txt(exon_snvs, txt_file)
        
            debug_out = open(outdir + "/filter_debug.tsv", 'w')
            self.output_groups(self.grouped_snvs, debug_out, debug=True)
            debug_out.close()
        
            novel_snvs = [snv for snv in filtered_snvs if snv.dbsnp == '-']
            txt_file = outdir + "/events_filtered_novel.tsv"
            self.output_txt(novel_snvs, txt_file) 
        
            novel_exon_snvs = [snv for snv in filtered_snvs if snv.exon and snv.nonsynon and snv.dbsnp == '-']
            txt_file = outdir + "/events_exons_novel.tsv"
            self.output_txt(novel_exon_snvs, txt_file) 
            
    def output_txt(self, snvs, outfile):
        """Groups and outputs given events in tabular format"""
        grouped_snvs = self.group(snvs)        
        out = open(outfile, 'w')
        out.write("\t".join(self.output_headers) + "\n")
        self.output_groups(grouped_snvs, out)
        out.close()
        
    def post_filter(self, min_from_end=None, no_mito=False):
        """Filters events based on read-support, mitonchondria, distance from contig edge, etc"""
        snvs_filtered = []
        
        # for extracting kmer length from contig name
        kmer = re.compile(r'k(\d+):')
        
        # filters out SNVs that are highly clustered
        self.filter_bad_region()

        for snv in self.snvs:
            if snv.nreads_event != 'na' and int(snv.nreads_event) >= int(self.min_reads):
                snv.enough_coverage = True
                
            if snv.nreads_genome != 'na' and int(snv.nreads_genome) < int(self.min_reads_genome):
                snv.enough_coverage = False
                        
            if snv.gene and ('exon' in snv.gene or 'utr' in snv.gene):
                snv.exon = True
                
            if snv.exon and not ':synon' in snv.gene:
                snv.nonsynon = True
                
            # min_from_end only applies for non-bubbles (contig name ended with \d)
            if re.match('\d', snv.var[-1]):
                if min_from_end == None:
                    m = kmer.search(snv.var)
                    if m and m.group:
                        min_from_end_contig = int(m.group(1))
                   
                        if int(snv.from_end) < min_from_end_contig:
                            snv.too_close_to_end = True
                        
                elif int(snv.from_end) < int(min_from_end):
                    snv.too_close_to_end = True
                    
            # mitochondria
            if no_mito and re.match('m', snv.ref, re.IGNORECASE):
                snv.artefact = True
                                
    def filter_bad_region(self):
        """Filters out SNVs that are highly clustered.
        Cannot have more than 10 in 500bp windown.
        """
        # group events by contig
        snvs_by_contig = {}
        for snv in self.snvs:
            if snv.snv_type != 'snv' or len(snv.var_seq) != 1:
                continue
            
            if not snvs_by_contig.has_key(snv.var):
                snvs_by_contig[snv.var] = []
                
            snvs_by_contig[snv.var].append(snv)
            
        for contig, snvs in snvs_by_contig.iteritems():
            snvs.sort(lambda x,y: int(x.var_start)-int(y.var_start))
            
            for i in range(len(snvs)):
                if i < len(snvs)-1:
                    diff = int(snvs[i+1].var_start) - int(snvs[i].var_start)
                else:
                    diff = 'na'
                
            max_in_window = 10
            window = 500
            for i in range(len(snvs)-max_in_window+1):
                bad = [i]
                for j in range(i+1, len(snvs)):
                    if int(snvs[j].var_start) - int(snvs[i].var_start) < window:
                        bad.append(j)
                    else:
                        break
                    
                if len(bad) > max_in_window:
                    for idx in bad:
                        snvs[idx].artefact = True
                        
    def get_snvs(self, align, splice_motifs=None, cutoff=None, no_segdup=False):
        """Extracts indels and snvs from alignment"""
        all_snvs = []
        sys.stderr.write("processing %s\n" % (align.query))
                                    
        if self.refseq:
            if self.sample_type == 'transcriptome' and self.fix_align:
                align.correct_blocks(splice_motifs, self.refseq, align.contig.sequence)
                
            if self.get_indel:
                all_snvs.extend(self.gap_snv(align, splice_motifs, align.contig.sequence, cutoff=cutoff))

            if self.get_snv and (align.mismatch is None or int(align.mismatch) > 0):
                all_snvs.extend(self.match_blocks(align, self.refseq, align.contig.sequence))
                
        if all_snvs and no_segdup:
            proper_chrom = tools.proper_chrom(align.target, chrom_proper=self.chrom_proper)
            if self.is_within_repeats(proper_chrom, [int(align.tstart), int(align.tend)]):
                print 'skip contig %s (%sbp): %s:%s-%s entirely with segdup/repeat' % (align.query, align.query_len, align.target, align.tstart, align.tend)
                del all_snvs[:]
            
        return all_snvs
    
    def match_intron(self, ss, splice_motifs):
        """Determines splite sites correspond to intron by comparing to splice motifs"""
        if ss and (splice_motifs.has_key(ss.lower()) or splice_motifs.has_key(tools.reverse_complement(ss).lower())):
            return True
        else:
            return False
        
    def gap_snv(self, align, splice_motifs, query_seq, cutoff=None):
        """Identifies insertions, deletions, inversion from gapped alignments"""
        if self.debug:
            print align.target, align.blocks
            print align.query, align.query_blocks
            print align.splice_sites
            
        snvs = []
        # cannot identify indels without splice site information
        if self.sample_type == 'transcriptome' and not align.splice_sites:
            return snvs
            
        for i in range(len(align.blocks)-1):
            if self.sample_type != 'transcriptome' or not self.match_intron(align.splice_sites[i], splice_motifs):
                if align.query_strand == '+':
                    qstart = align.query_blocks[i][1]+1
                    qend = align.query_blocks[i+1][0]-1
                    query = query_seq[qstart-1:qend]
                else:
                    qend = align.query_blocks[i][1]-1
                    qstart = align.query_blocks[i+1][0]+1
                    query = query_seq[qstart-1:qend]
                    query = tools.reverse_complement(query)
                    
                # target strand always + from psl
                tstart = align.blocks[i][1]+1
                tend = align.blocks[i+1][0]-1
                target = ''
                if tstart <= tend:
                    target = self.refseq.GetSequence(align.target, tstart, tend)
                
                #if code cannot extract sequence from reference, there must be a disagreement between alignment and reference - abort analysis
                if tend > tstart-1 and len(target) < 1:
                    sys.stderr.write("cannot extract reference sequence, abort: %s %s %s\n" % (align.target, tstart-1, tend))
                    sys.exit(100)

                snv_type = None               
                if qstart > qend and (tend - tstart) >= 0:
                    size = tend - tstart + 1
                    if align.query_strand == '+':
                        qstart = qend
                    else:
                        qend = qstart
                    snv_type = "del"                    
                elif tstart > tend and (qend - qstart) >= 0:
                    size = qend - qstart + 1
                    tstart = tend
                    snv_type = "ins"                    
                else:
                    size = min(1, tend - tstart + 1)
                    snv_type = "indel"
                    
                # skip if 0 or negative size event detected (or smaller than cutoff)
                if size <= 0 or (cutoff and size > cutoff):
                    continue

                target = target.lower()
                query = query.lower()
                # would not report event with non-AGCT characters
                if not re.search('[^agtcATGC]', target) and not re.search('[^agtcATGC]', query):
                    if snv_type != 'indel':
                        snv = SNV('psl', snv_type, align.target, tstart, tend, target, align.query_strand, align.query, qstart, qend, query)
                        snvs.append(snv)
                    # resolves indels
                    else:
                        if len(query) == len(target) and\
                           (query[::-1].lower() == target.lower() or tools.reverse_complement(query).lower() == target.lower()):
                            # inversion must be longer than 1 base
                            if len(query) > 1:
                                snv = SNV('psl', 'inv', align.target, tstart, tend, target, align.query_strand, align.query, qstart, qend, query)
                                snvs.append(snv)
                            # 1 bp gap in both query and target == snv
                            else:
                                snv = SNV('psl', 'snv', align.target, tstart, tend, target, align.query_strand, align.query, qstart, qend, query)
                                snvs.append(snv)
                        # breaks up indel into ins and del
                        else:
                            if align.query_strand == '+':
                                qcoord = qstart
                            else:
                                qcoord = qend
                            
                            snv = SNV('psl', 'del', align.target, tstart, tend, target, align.query_strand, align.query, qcoord, qcoord, query)
                            snvs.append(snv)
                            
                            tcoord = tstart
                            snv = SNV('psl', 'ins', align.target, tcoord, tcoord, target, align.query_strand, align.query, qstart, qend, query)
                            snvs.append(snv)
                                            
        return snvs

    def parse_results(self, snv_file, select_types=None, chrom=None):
        """Parses results from single file into SNV objects"""        
        names = SNVCaller.output_headers
        # conversion between header name and object attribute
        field_name_conversion = {
            'type': 'snv_type',
            'chr': 'ref',
            'chr_start': 'ref_start',
            'chr_end': 'ref_end',
            'ctg': 'var',
            'ctg_len': 'var_len',
            'ctg_start': 'var_start',
            'ctg_end': 'var_end',
            'len': 'snv_len',
            'ref': 'ref_seq',
            'alt': 'var_seq',
            'event_reads': 'nreads_event',
            'contig_reads': 'nreads_contig',
            'genome_reads': 'nreads_genome',
            'gene': 'gene',
            'from_end': 'from_end',
            'ctg_strand': 'query_strand',
        }

        for line in open(snv_file, 'r'):
            cols = line.rstrip('\n').split('\t')
            
            if cols[0] == 'id':
                continue

            attributes = {}
            for i in range(1, len(cols)):
                name = names[i]
                value = cols[i]
                
                if field_name_conversion.has_key(name):
                    name = field_name_conversion[name]
                    
                if name in ('expansion', 'from_end'):
                    value = int(value)
                elif name == 'confirm_contig_region':
                    value = value.split('-')
                    value[0] = int(value[0])
                    value[1] = int(value[1])
                elif name == 'at_least_1_read_opposite':
                    if value == 'true':
                        value = True
                    else:
                        value = False
                    
                attributes[name] = value
            
            if select_types and not attributes['snv_type'] in select_types:
                continue
            
            if chrom and attributes['ref'] != chrom:
                continue
            
            snv = SNV(method='psl')
            tools.set_attrs(snv, attributes)
            self.snvs.append(snv)
            
    def parse_results_dir(self, path):
        """Parses results into single file given a directory of output directories"""
        output_dirs = sorted(glob.glob(os.path.join(path, '*')))
                
        output_files = []
        num_output_dirs = 0
        missing_dirs = []
        
        ok = True
        
        for job_num in range(1, len(output_dirs)+1):
            cluster_outdir = "%s/%s" % (options.output_file, job_num)
            if os.path.isdir(cluster_outdir):
                num_output_dirs += 1
                    
                snv_file = cluster_outdir + '/events.tsv'
                if os.path.exists(snv_file):
                    output_files.append(snv_file)
                else:
                    missing_dirs.append(str(job_num))
                    print snv_file
                    sys.stdout.write("%s does not have output\n" % (cluster_outdir))
                    ok = False
                        
        sys.stdout.write("output dirs:%s output files:%s\n" % (num_output_dirs, len(output_files)))
                
        if num_output_dirs == len(output_files):
            concat_outfile = options.outdir + "/events_concat.tsv"
            print concat_outfile
            concat_out = open(concat_outfile, 'w')
                    
            # just for check if header is to be written
            count = 0
            for output_file in output_files:
                infile = open(output_file, 'r')
                if count == 0:
                    concat_out.writelines(open(output_file, 'r').readlines())
                else:
                    concat_out.writelines(open(output_file, 'r').readlines()[1:])
                infile.close()
                    
                self.parse_results(output_file)
                count += 1
                        
            concat_out.close()
        else:
            sys.stdout.write("missing (%s):%s\n" % (len(missing_dirs), ','.join(missing_dirs)))
            ok = False
        
        return ok
            
    def match_blocks(self, align, query_seq):
        """Identifies SNVs"""
        snvs = []
        
        for i in range(len(align.blocks)):
            if align.query_strand == '+':
                qseq = query_seq[int(align.query_blocks[i][0])-1:int(align.query_blocks[i][1])]
            else:
                qseq = tools.reverse_complement(query_seq[int(align.query_blocks[i][1]-1):int(align.query_blocks[i][0])])
            tseq = self.refseq.GetSequence(align.target, int(align.blocks[i][0]), int(align.blocks[i][1]))
            
            mismatches = self.find_mismatches(qseq, tseq)
            
            for pos, change in mismatches.iteritems():
                tpos = int(align.blocks[i][0]) + pos
                if int(align.query_blocks[i][0]) < int(align.query_blocks[i][1]):
                    qpos = int(align.query_blocks[i][0]) + pos
                else:
                    qpos = int(align.query_blocks[i][0]) - pos
                snv = SNV('psl', 'snv', align.target, tpos, tpos, change[0], align.query_strand, align.query, qpos, qpos, change[1])
                snvs.append(snv)
                        
        return snvs

    def find_mismatches(self, qseq, tseq):
        """Reports substitutions given query and target sequence of same length"""
        pos = {}
        bases = ['a','g','t','c']
        if len(qseq) == len(tseq):
            for i in range(len(qseq)):
                if qseq[i].lower() != tseq[i].lower() and qseq[i].lower() in bases and tseq[i].lower() in bases:
                    pos[i] = [tseq[i].lower(), qseq[i].lower()]
        return pos
                    
def main(args, options):
    # output log file
    tools.output_log(__version__, sys.argv, [vars(options)], options.outdir)

    snv_caller = SNVCaller(options.genome, 
                           get_snv=options.snv, 
                           get_indel=options.indel, 
                           sample_type=options.sample_type,
                           min_reads=options.min_reads_contigs, 
                           min_from_end=options.min_from_end, 
                           gene_model=options.gene_model, 
                           debug=options.debug,
                           annodir=options.annodir,
                           mmcfg=options.mmcfg)
    snv_caller.min_reads_genome = options.min_reads_genome
    snv_caller.fix_align = options.fix_align
    
    bm = None
    if options.bubble2contig:
        bm = BubbleMapper(options.bubble2contig)
        bm.map_to_contigs()
        snv_caller.bubble_mapping = bm

    # creates output directories
    if options.outdir and not os.path.isdir(options.outdir):
        os.makedirs(options.outdir)
                   
    # extracts alignments
    if options.align_file:        
        if os.path.isdir(options.align_file):
            align_files = sorted(glob.glob(os.path.join(options.align_file, '*.*')))
                
            for i in range(len(align_files)):
                snv_caller.align_file = "%s/seq.%s.psl" % (options.align_file, i+1)
                    
                if os.path.isdir(options.contigs_file):
                    snv_caller.contigs_file = "%s/seq.%s.fa" % (options.contigs_file, i+1)
                else:
                    snv_caller.contigs_file = options.contigs_file
                        
                sys.stdout.write("processing alignment file #%d\n" % (i+1))
                snv_caller.extract(cutoff=options.indel_size, 
                                   no_group=True, 
                                   match_percent=options.match_percent, 
                                   identity=options.identity,
                                   no_segdup=options.no_segdup)                                
        else:
            snv_caller.align_file = options.align_file
            snv_caller.contigs_file = options.contigs_file
            snv_caller.extract(cutoff=options.indel_size, 
                               match_percent=options.match_percent, 
                               identity=options.identity, 
                               no_segdup=options.no_segdup)
                
    # extracts results
    elif options.output_file:        
        if os.path.isdir(options.output_file):
            ok = snv_caller.parse_results_dir(options.output_file)
            if not ok:
                sys.exit(1)
        else:
            snv_caller.parse_results(options.output_file)  
                           
    if snv_caller:
        # group events
        snv_caller.grouped_snvs = snv_caller.group(snv_caller.snvs)
        
        # annotate
        snv_caller.annotate_genes()
                   
        # gather read support whenever bam file(s) given
        if (options.contigs_bam and os.path.exists(options.contigs_bam)) or\
           (options.genome_bam and os.path.exists(options.genome_bam)):
            snv_caller.add_support(from_end=options.read_buffer, genome_bamfile=options.genome_bam, contigs_bamfile=options.contigs_bam)

        # add contig support reads if desired
        if options.add_contig_reads:
            snv_caller.add_contig_reads(snv_caller.grouped_snvs)
            
        # overlaps repeats, dbSNP
        if options.olap_annot:
            snv_caller.overlap_repeats()
            snv_caller.overlap_dbsnp()
        
        # report
        snv_caller.report(options.outdir, post_filter=options.filter)
        
if __name__ == '__main__':
    usage = "Usage: %prog [options]"
    parser = OptionParser(usage=usage, version="%prog " + __version__)
    
    io = OptionGroup(parser, "input and output")
    io.add_option("-a", "--align_file", dest="align_file", help="alignment file")
    io.add_option("-c", "--contigs_file", dest="contigs_file", help="contigs file")
    io.add_option("-g", "--genome", dest="genome", help="genome")
    io.add_option("-o", "--outdir", dest="outdir", help="outdir")
    io.add_option("-f", "--output_file", dest="output_file", help="output file")
    io.add_option("-T", "--type", dest="sample_type", help="sample type: transcriptome(default), genome, exome", default='transcriptome')
    io.add_option("-B", "--bubble2contig", dest="bubble2contig", help="bubble to contig alignment")
    parser.add_option_group(io)
    
    annotations = OptionGroup(parser, "annotations")
    annotations.add_option("--annodir", dest="annodir", help="the Trans-ABySS 'annotation' directory")
    annotations.add_option("-m", "--gene_model", dest="gene_model", help="gene model used for annotation e.g. k where k=known_genes, e=ensembl, r=refseq")
    annotations.add_option("-O", "--olap_annot", dest="olap_annot", help="overlaps repeats and dbsnp when parsing events", action="store_true", default=False)
    annotations.add_option("--mmcfg", dest="mmcfg", help="the path of `model_matcher.cfg'")
    parser.add_option_group(annotations)
    
    read_support = OptionGroup(parser, "read support")
    read_support.add_option("-G", "--genome_bam", dest="genome_bam", help="genomic bam file")
    read_support.add_option("-C", "--contigs_bam", dest="contigs_bam", help="contigs bam file")
    read_support.add_option("-r", "--min_reads_contigs", dest="min_reads_contigs", help="minimum number of read support from contig bam. Default=3", default=3, type='int')
    read_support.add_option("-n", "--min_reads_genome", dest="min_reads_genome", help="minimum number of read support from genome bam. Default=1", default=1, type='int')
    read_support.add_option("-b", "--read_buffer", dest="read_buffer", help="minimum number of bases from read end", default=8, type='int')
    read_support.add_option("-S", "--add_contig_reads", dest="add_contig_reads", action="store_true", default=False)
    parser.add_option_group(read_support)
    
    operation = OptionGroup(parser, "operations")
    operation.add_option("-V", "--snv", dest="snv", help="extract SNV", action="store_true", default=False)
    operation.add_option("-L", "--indel", dest="indel", help="extract indel", action="store_true", default=False)
    operation.add_option("-X", "--filter", dest="filter", help="post-filter", action="store_true", default=False)
    operation.add_option("-d", "--debug", dest="debug", help="debug", action="store_true", default=False)
    operation.add_option("-A", "--fix_align", dest="fix_align", help="fix alignment block coords", action="store_true", default=False)
    parser.add_option_group(operation)
    
    params = OptionGroup(parser, "parameters and filtering")
    params.add_option("-i", "--indel_size", dest="indel_size", help="indel upper size limit for blat", type='int') 
    params.add_option("-E", "--contig_end", dest="min_from_end", help="minimum distance from contig end", type='int')
    params.add_option("-M", "--match_percent", dest="match_percent", help="minimum match percent for screening alignments. Default:95", type='float', default=95.0)
    params.add_option("-I", "--identity", dest="identity", help="minimum identity for screening alignments Default:95", type='float', default=95.0)
    params.add_option("-N", "--no_mito", dest="no_mito", help="skip mitochrondria", action="store_true", default=False)
    params.add_option("-U", "--no_segdup", dest="no_segdup", help="skip contigs lying completely within segdup/simple_repeat", action="store_true", default=False)
    parser.add_option_group(params)
    
    (options, args) = parser.parse_args()
    main(args, options)
    
