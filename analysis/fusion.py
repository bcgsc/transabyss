"""
This module finds structural variants and gene fusions through split alignments

Author: Readman Chiu rchiu@bcgsc.ca
"""
__version__ = '1.5.1'

import os
import sys
import re
import glob
import operator
from optparse import OptionParser, OptionGroup
from math import floor
from utilities.align_parsers import psl, sam
from utilities import intspan, tools, cfg
from utilities.alignment import Alignment
from utilities.assembly import Assembly
from utilities.bam import BAM
from utilities.tools import compare_chr, get_splice_motifs
from feature import FeatureFinder
from transcript import Transcript
from cnv import CNV
from annotations import dbsnp, dgv, repeat
from pep_change import translate
import utilities.cfg

PYTHON = 'python'

class Fusion:
    """Calls fusion events from split contig alignments"""
    def __init__(self, a1, a2):
        self.align1 = a1
        self.align2 = a2

        self.contig = a1.query
        self.contig_size = a1.query_len

        self.gene1 = self.gene2 = 'NA'
        self.sense1 = self.sense2 = 'NA'
        self.exon1 = self.exon2 = 'NA'
        self.txt1 = self.txt2 = None
        self.txt1_name = self.txt2_name = 'NA'
        self.exon_bound1 = self.exon_bound2 = 'NA'
        self.cross_exon1 = self.exon_exon2 = False

        self.overlap_query = None
        self.overlap_query_fraction = None
        self.overlap_target = None
        self.overlap_target_fraction = None
        self.contig_cover = None
        self.contig_cover_fraction = None
        
        #read evidence
        self.num_read_pairs = 0
        self.num_flanking_pairs = 0
        self.num_breakpoint_pairs = [0,0]
        self.num_spanning_reads = 0
        self.num_spanning_reads_forward, self.num_spanning_reads_reverse = 0, 0
        
        self.mutation = 'NA'
        self.breakpoint = 'NA'
        self.lsr_region = None
        self.cnv = None
        self.fusion_type = 'LSR'
        
        #filtering
        self.exclude = False
        self.good_read_pair = False
        self.good_span_read = False
        self.read_support = False
        
        self.contig_local = None
        self.local = None
        self.ff = None
        self.gene_fusion = False
        self.lsr = False
        self.orients = []
        self.reciprocal = None
        self.event_id = None
        self.descriptor = 'NA'
        self.breaks = None
        self.chroms = None
        
        self.dbsnp = 'NA'
        self.dgv = 'NA'
        
        self.gene5prime = 'NA'
        self.exon5prime = 'NA'
        self.gene3prime = 'NA'
        self.exon3prime = 'NA'
        self.align5prime = self.align3prime = None
        self.genes_flipped = False
        
        self.txt5prime = self.txt3prime = None
        self.frame = 'NA'
        
        self.contig_seq = None
        self.contig_reverse = None
        self.probe = 'NA'
        
        self.repeat1 = self.repeat2 = 'NA'
        
        self.artefact = False
        
        self.inserted_seq = ''
        self.micro_homol_seq1 = self.micro_homol_seq2 = ''
        self.micro_homol1 = self.micro_homol2 = []
                        
    def details(self, tsv=True, id=None):
        """Outputs event in tab format"""
        ctg = self.contig
        
        align1, align2 = self.align1, self.align2            
        target = "%s:%s-%s,%s:%s-%s" % (align1.target, align1.tstart, align1.tend, align2.target, align2.tstart, align2.tend)
        query = "%s-%s,%s-%s" % (align1.qstart, align1.qend, align2.qstart, align2.qend)        
        align_strands = "NA,NA"
        if align1.strand:
            align_strands = "%s,%s" % (align1.strand, align2.strand)
        align_len1 = float(align1.qend) - float(align1.qstart) + 1
        align_len2 = float(align2.qend) - float(align2.qstart) + 1
        align_frac1 = align_len1/float(align1.query_len)
        align_frac2 = align_len2/float(align2.query_len)
        params = "TO:%.2f,CO:%.2f,CC:%.2f,I1:%.1f,I2:%.1f,AF1:%.2f,AF2:%.2f" % (self.overlap_target_fraction, self.overlap_query_fraction, self.contig_cover_fraction, float(align1.identity), float(align2.identity), align_frac1, align_frac2)
            
        fusion_type = 'NA'
        if self.gene_fusion:
            fusion_type = 'gene_fusion'
        elif self.lsr:
            fusion_type = 'lsr'
        elif self.local:
            fusion_type = 'local'            
        reciprocal = 'NA'
        if self.reciprocal:
            reciprocal = self.reciprocal
                        
        if self.orients:
            if self.orients[0] and self.orients[1]:
                orient1, orient2 = 'L', 'L'
            elif self.orients[0] and not self.orients[1]:
                orient1, orient2 = 'L', 'R'
            elif not self.orients[0] and self.orients[1]:
                orient1, orient2 = 'R', 'L'
            elif not self.orients[0] and not self.orients[1]:
                orient1, orient2 = 'R', 'R'
        else:
            orient1, orient2 = 'NA', 'NA'
        orients = '%s,%s' % (orient1, orient2)
        
        genes = '%s,%s' % (self.gene1, self.gene2)
        txt1, txt2 = self.txt1_name, self.txt2_name
        if self.txt1:
            txt1 = self.txt1.name
        if self.txt2:
            txt2 = self.txt2.name
        transcripts = '%s,%s' % (txt1, txt2)
        senses = '%s,%s' % (self.sense1, self.sense2)
        exons = '%s,%s' % (self.exon1, self.exon2)
        exon_bounds = '%s,%s' % (self.exon_bound1, self.exon_bound2)
            
        data = [str(self.event_id), str(ctg), str(self.contig_size), target, query, align_strands, str(self.num_flanking_pairs), ','.join([str(self.num_breakpoint_pairs[0]), str(self.num_breakpoint_pairs[1])]), str(self.num_spanning_reads), str(self.num_spanning_reads_forward), str(self.num_spanning_reads_reverse), self.mutation, self.breakpoint, str(self.size()), genes, transcripts, senses, exons, exon_bounds, reciprocal, self.descriptor, orients, self.gene5prime, self.gene3prime, self.exon5prime, self.exon3prime, self.frame, self.probe, self.repeat1, self.repeat2, params, self.fusion_type]                                               
        data.append(self.dbsnp)
        data.append(self.dgv)
        
        if self.cnv:
            data.append(self.cnv)
            
        if tsv:
            return "\t".join(data)
                     
    def report_align_params(self):
        """Reports alignment parameters"""
        align_len1 = float(self.align1.qend) - float(self.align1.qstart) + 1
        align_len2 = float(self.align2.qend) - float(self.align2.qstart) + 1
        align_frac1 = align_len1/float(self.align1.query_len)
        align_frac2 = align_len2/float(self.align2.query_len)
        params = "TO:%.2f,CO:%.2f,CC:%.2f,I1:%.1f,I2:%.1f,AF1:%.2f,AF2:%.2f" % (self.overlap_target_fraction, self.overlap_query_fraction, self.contig_cover_fraction, float(self.align1.identity), float(self.align2.identity), align_frac1, align_frac2)
        
        return params
    
    def has_mito(self):
        """Checks if fusion alignments belongs to mitochondria"""
        for target in (self.align1.target, self.align2.target):
            if target.lstrip('chr').rstrip('T') == 'M':
                return True
            
        return False
    
    def find_read_pairs(self, bam, genome_coord_buffer=0, breakpoint_buffer=0, keep_mito=False, reads=None, only_unique=False, debug=False):
        """Finds read pairs supporting event in reads-to-genome BAM file"""
        if not keep_mito and self.has_mito():
            return
        
        if not self.breaks:
            return
        
        # make sure bam file has chrom
        chrom1 = self.align1.target
        if not chrom1 in bam.bam.references:
            if 'chr' in chrom1 and chrom1.replace('chr', '') in bam.bam.references:
                chrom1 = chrom1.replace('chr', '') 
            elif not 'chr' in chrom1 and 'chr' + chrom1 in bam.bam.references:
                chrom1 = 'chr' + chrom1
        chrom2 = self.align2.target
        if not chrom2 in bam.bam.references:
            if 'chr' in chrom2 and chrom2.replace('chr', '') in bam.bam.references:
                chrom2 = chrom2.replace('chr', '') 
            elif not 'chr' in chrom2 and 'chr' + chrom2 in bam.bam.references:
                chrom2 = 'chr' + chrom2
                        
        breakpoint = [chrom1, int(self.breaks[0])], [chrom2, int(self.breaks[1])]        
        if not self.orients[0]:
            upper_bound = max(int(self.align1.tend), int(self.breaks[0]) + genome_coord_buffer)
            region1 = [chrom1, int(self.breaks[0]), upper_bound]                   
        else:
            lower_bound = min(int(self.align1.tstart), int(self.breaks[0]) - genome_coord_buffer)
            region1 = [chrom1, lower_bound, int(self.breaks[0])]
        if not self.orients[1]:
            upper_bound = max(int(self.align2.tend), int(self.breaks[1]) + genome_coord_buffer)
            region2 = [chrom2, int(self.breaks[1]), upper_bound]                   
        else:
            lower_bound = min(int(self.align2.tstart), int(self.breaks[1]) - genome_coord_buffer)
            region2 = [chrom2, lower_bound, int(self.breaks[1])]      
                                                        
        # finds full flanking pairs
        focal = False
        if self.fusion_type == 'ITD':
            focal = True
            
        no_proper = False
        if self.fusion_type == 'ITD' or self.fusion_type == 'PTD':
            no_proper = True
                        
        self.num_flanking_pairs, pairs = bam.confirm(None, None, feature="breakpoint", regions=[region1, region2], breakpoint=breakpoint, breakpoint_buffer=breakpoint_buffer, focal=focal, seq=debug, no_proper=no_proper)
        
        if debug:
            if pairs:
                for read1, read2 in pairs:
                    print '>%s-1 flanking-pair %s %s %s reverse:%s proper_pair:%s read1:%s is_duplicate:%s\n%s' % (read1.qname, self.align1.query, self.breakpoint, read1.pos, read1.is_reverse, read1.is_proper_pair, read1.is_read1, read1.is_duplicate, read1.seq)
                    print '>%s-2 flanking-pair %s %s %s reverse:%s proper_pair:%s read1:%s is_duplicate:%s\n%s' % (read2.qname, self.align1.query, self.breakpoint, read2.pos, read2.is_reverse, read2.is_proper_pair, read2.is_read1, read2.is_duplicate, read2.seq)
        
        # finds breakpoint pairs - pairs with one mate overlapping junctions
        if reads:
            reads1 = reads2 = []  
            if reads.has_key(self.breaks[0]) and reads[self.breaks[0]]:
                reads1 = bam.find_mates_in_genome(region1, [self.align1.target, self.breaks[0]], missing_mates=reads[self.breaks[0]], breakpoint_buffer=breakpoint_buffer, no_duplicates=False, no_chastity=False)        
            if reads.has_key(self.breaks[1]) and reads[self.breaks[1]]:
                reads2 = bam.find_mates_in_genome(region2, [self.align2.target, self.breaks[1]], missing_mates=reads[self.breaks[1]], breakpoint_buffer=breakpoint_buffer, no_duplicates=False, no_chastity=False)
                            
            reads1_dict = dict((read.qname, read) for read in reads1)
            reads2_dict = dict((read.qname, read) for read in reads2)
            mutual = [read for read in reads2 if read.qname in reads1_dict]
            if mutual:
                for read in mutual:
                    if abs(read.pos - self.breaks[0]) < abs(read.pos - self.breaks[1]):
                        reads2.remove(reads2_dict[read.qname])
                    else:
                        reads1.remove(reads1_dict[read.qname])
                                            
            if len(reads1) > 0:                
                self.num_breakpoint_pairs[0] = len(reads1)
                
                if only_unique:
                    unique = self.find_unique(reads1)
                    self.num_breakpoint_pairs[0] = len(unique.keys())
                    
                if debug:
                    for read in reads1:
                        print '>%s breakpoint-pair1 %s %s:%s %s %s\n%s' % (read.qname, self.align1.query, self.align1.target, self.breaks[0], read.pos, read.is_reverse, read.seq)

            if len(reads2) > 0:
                self.num_breakpoint_pairs[1] = len(reads2)
                                
                if only_unique:
                    unique = self.find_unique(reads2)
                    self.num_breakpoint_pairs[1] = len(unique.keys())
                    
                if debug:
                    for read in reads2:
                        print '>%s breakpoint-pair2 %s %s:%s %s %s\n%s' % (read.qname, self.align1.query, self.align2.target, self.breaks[1], read.pos, read.is_reverse, read.seq)
                                                
        self.num_read_pairs = self.total_read_pairs()
        
    def find_spanning_reads(self, bam, breakpoint_buffer=0, debug=False, max_depth=None, only_unique=False):
        """Finds reads spanning events in reads-to-contigs BAM file"""
        contig_coords = [int(self.align1.qstart), int(self.align1.qend), int(self.align2.qstart), int(self.align2.qend)]
        contig_coords.sort(key = int)
        # makes sure coordinates do not go overboard in after adding buffer
        start = max(1, contig_coords[1] - breakpoint_buffer)
        end = min(int(self.align1.query_len), contig_coords[2] + breakpoint_buffer)
        
        spanning_reads = bam.confirm(target=self.align1.query, pos=[start,end], feature="breakpoint", seq=True, max_depth=max_depth)
        self.num_spanning_reads = len(spanning_reads)
                
        # only count unique spanning reads
        if only_unique:
            unique = self.find_unique(spanning_reads)
            self.num_spanning_reads = len(unique.keys())
            
            if self.num_spanning_reads != len(spanning_reads):
                spanning_reads = []
                for key, reads in unique.iteritems():
                    spanning_reads.append(reads[0])
                          
        # calculates number reads in each orientation
        if self.contig_reverse is not None:
            strands = {'forward':[], 'reverse':[]}
            for read in spanning_reads:
                if not self.contig_reverse:
                    if not read.is_reverse:
                        strand = 'forward'
                    else:
                        strand = 'reverse'
                    
                elif self.contig_reverse:
                    if not read.is_reverse:
                        strand = 'reverse'
                    else:
                        strand = 'forward'
                        
                strands[strand].append(read)
                
                if debug:
                    print '>%s spanning %s %s %s %s\n%s' % (read.qname, self.align1.query, self.breakpoint, read.pos, read.is_reverse, read.seq)
            
            self.num_spanning_reads_forward = len(strands['forward'])
            self.num_spanning_reads_reverse = len(strands['reverse'])
    
        return spanning_reads
                
    def find_unique(self, reads):
        """Groups reads with same alignment"""
        unique = {}
        for read in reads:
            strand = '+'
            if read.is_reverse:
                strand = '-'
            key = ' '.join([str(read.pos), strand, str(read.alen)])
            if not unique.has_key(key):
                unique[key] = []
            unique[key].append(read)
                    
        return unique

    def is_local(self, model=None, annot_file=None, local_gap=5000):
        """Determines if event is local (presumably experimental artefacts)"""
        # can't be local if aligned to different chromosomes
        if self.align1.target != self.align2.target:
            self.local = False
            return self.local
        
        # don't ever exclude ITD and PTD events although they are local
        if self.fusion_type in ('ITD', 'PTD'):
            self.local = False
            return self.local
        
        # if events overlap by a small margin (5k)
        if self.align1.target == self.align2.target:            
            spans = [[int(self.align1.tstart), int(self.align1.tend)], [int(self.align2.tstart), int(self.align2.tend)]]
            if int(self.align1.tstart) > int(self.align2.tstart):
                spans.reverse()
            spans[0][1] += local_gap

            if intspan.overlap(spans[0], spans[1]):
                self.local = True

        if self.gene1 and self.gene2 and self.gene1.upper() != 'NA' and self.gene2.upper() != 'NA':
            # different-gene event won't be local
            if self.gene1 != self.gene2:
                self.local = False
                return self.local
            
            # if breakpoints mapped to exons that would lead to sense fusion, don't call local
            if self.exon5prime and self.exon3prime:
                if self.exon5prime != self.exon3prime:
                    self.local = False
                    return self.local
                
            # if breakpoints mapped to introns, could be sense fusion, don't call local
            if self.exon1 and self.exon2:
                if 'intron' in self.exon1 and 'intron' in self.exon2:
                    if self.exon1 != self.exon2:
                        self.local = False
                        return self.local
            
            if not self.local:
                if self.gene1 == self.gene2:
                    self.local = True
                    self.exclude = True

            if not self.local and model and annot_file:
                if not self.txt1:
                    txts = Transcript.find_by_name(model, annot_file, self.gene1)
                    if txts:
                        txts.sort(lambda x,y:y.length - x.length)
                        self.txt1 = txts[0]

                if not self.txt2:
                    txts = Transcript.find_by_name(model, annot_file, self.gene2)
                    if txts:
                        txts.sort(lambda x,y:y.length - x.length)
                        self.txt2 = txts[0]

                # if largest transcript overlapping gene1 and largest transcript overlapping gene2, call local
                if self.txt1 and self.txt2 and self.txt1.chrom == self.txt2.chrom and intspan.overlap((self.txt1.txStart, self.txt1.txEnd), (self.txt2.txStart, self.txt2.txEnd)):
                    self.local = True

        return self.local
    
    def define_mutation(self):
        """Defines event type
        Algorithm:
        1. sort the target coordinates based on query strands
        2. define mutation based on whether target coordinates ascend or descend in each region
        3. define and keep breakpoints in each target region: self.breaks
        4. formulate breakpoint
        5. define 'orientation' in each segment for finding reciprocal events
        6. give proper genetic name of event if cytobands provided
        """        
        if self.align1.target != self.align2.target:
            self.mutation = 'translocation'
            
            span1 = [int(self.align1.tstart), int(self.align1.tend)]
            span2 = [int(self.align2.tstart), int(self.align2.tend)]
            if self.align1.strand == '-':
                span1.reverse()
            if self.align2.strand == '-':
                span2.reverse()
                
            if int(self.align1.qstart) < int(self.align2.qstart):
                self.breakpoint = "%s:%s|%s:%s" % (self.align1.target, span1[1], self.align2.target, span2[0])
                self.breaks = [span1[1], span2[0]]
                self.chroms = [self.align1.target, self.align2.target]
                    
                if span1[1] == min(span1):
                    self.orients.append(False)
                else:
                    self.orients.append(True)
                if span2[0] == min(span2):
                    self.orients.append(False)
                else:
                    self.orients.append(True)
                    
            else:
                self.breakpoint = "%s:%s|%s:%s" % (self.align1.target, span1[0], self.align2.target, span2[1])
                self.breaks = [span1[0], span2[1]]
                self.chroms = [self.align1.target, self.align2.target]
                    
                if span1[0] == min(span1):
                    self.orients.append(False)
                else:
                    self.orients.append(True)
                if span2[1] == min(span2):
                    self.orients.append(False)
                else:
                    self.orients.append(True)
                    
            span1.sort(lambda x,y: x-y)
            span1.sort(lambda x,y: x-y)
            self.lsr_region = [[self.align1.target, span1[0], span1[1]], [self.align2.target, span2[0], span2[1]]]
                                    
        else:
            self.chroms = (self.align1.target, self.align2.target)
                            
            if int(self.align1.qstart) < int(self.align2.qstart):
                span1 = [int(self.align1.tstart), int(self.align1.tend)]
                span2 = [int(self.align2.tstart), int(self.align2.tend)]
                if self.align1.strand == '-':
                    span1.reverse()
                if self.align2.strand == '-':
                    span2.reverse()
            else:
                span1 = [int(self.align2.tstart), int(self.align2.tend)]
                span2 = [int(self.align1.tstart), int(self.align1.tend)]
                if self.align2.strand == '-':
                    span1.reverse()
                if self.align1.strand == '-':
                    span2.reverse()
            
            ascend1 = None
            if span1[0] < span1[1]:
                ascend1 = True
            else:
                ascend1 = False                    
            ascend2 = None
            if span2[0] < span2[1]:
                ascend2 = True
            else:
                ascend2 = False
                                                  
            if span1[0] < span2[0]:
                if ascend1 and ascend2:
                    if intspan.overlap(span1, span2):
                        self.mutation = 'duplication'
                    else:
                        self.mutation = 'deletion'
                    self.orients = [True, False]
                    
                elif not ascend1 and not ascend2:
                    self.mutation = 'duplication'
                    self.orients = [False, True]
                    
                elif ascend1 and not ascend2:
                    self.mutation = 'inversion'
                    self.orients = [True, True]
                    
                elif not ascend1 and ascend2:
                    self.mutation = 'inversion'
                    self.orients = [False, False]
                    
                self.breakpoint = "%s:%s|%s:%s" % (self.align1.target, span1[1], self.align2.target, span2[0])
                self.lsr_region = [self.align1.target, span1[1], span2[0]]
                self.breaks = [span1[1], span2[0]]
                                
            elif span1[0] > span2[0]:                   
                if ascend1 and ascend2:
                    self.mutation = 'duplication'
                    self.orients = [False, True]
                    
                elif not ascend1 and not ascend2:
                    if intspan.overlap(span1[::-1], span2[::-1]):
                        self.mutation = 'duplication'
                    else:
                        self.mutation = 'deletion'
                    self.orients = [True, False]
                    
                elif ascend2 and not ascend1:
                    self.mutation = 'inversion'
                    self.orients = [False, False]
                    
                elif not ascend2 and ascend1:
                    self.mutation = 'inversion'
                    self.orients = [True, True]
                                    
                self.breakpoint = "%s:%s|%s:%s" % (self.align2.target, span2[0], self.align1.target, span1[1])
                self.lsr_region = [self.align1.target, span2[0], span1[1]]
                self.breaks = [span2[0], span1[1]]
                                                
            #probably not real, mis-assembly
            elif span1[0] == span2[0]:
                self.mutation = 'NA'
                self.breakpoint = 'NA'
                self.lsr_region = []
                
        if self.size() < 1:
            self.artefact = True
            
        # reverse alignments if necessary
        if self.breaks:
            # make sure breakpoints match alignments
            if (self.breaks[0] == int(self.align2.tstart) or self.breaks[0] == int(self.align2.tend)) and\
               (self.breaks[1] == int(self.align1.tstart) or self.breaks[1] == int(self.align1.tend)):
                tmp_align = self.align1
                self.align1 = self.align2
                self.align2 = tmp_align
            
            # order alignments by breakpoints, so that all events can be consistent for grouping
            self.reorder_aligns_by_coordinate()

                
    def reorder_aligns_by_coordinate(self):
        """Reverses alignments if first alignment is after second alignment in chromosome or coordinate"""
        if self.breaks and len(self.breaks) == 2 and \
           ((self.align1.target != self.align2.target and compare_chr(self.align1.target, self.align2.target) > 0) or
            (self.align1.target == self.align2.target and int(self.breaks[0]) > int(self.breaks[1]))):
            self.breaks.reverse()
            self.lsr_region.reverse()
            self.orients.reverse()
            
            tmp = self.breakpoint.split('|')
            tmp.reverse()
            self.breakpoint = '|'.join(tmp)
            tmp_align = self.align1
            self.align1 = self.align2
            self.align2 = tmp_align
            
            return True
        
        return False

    def proper_name(self, cytobands, chrom_proper=None):
        """Identifies overlappgin cytogenetic band and creates proper genetic name for event"""        
        if self.mutation and self.breaks and cytobands:
            chrom1, chrom2 = self.align1.target, self.align2.target
            if chrom_proper is not None:
                chrom1 = tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper)
                chrom2 = tools.proper_chrom(self.align2.target, chrom_proper=chrom_proper)
            
            band1, band2 = tools.overlap_cytobands(cytobands, chrom1, self.breaks[0], chrom2, self.breaks[1])
            
            #inversion across centromere becomes translocation
            if self.mutation == 'inversion' and band1 is not None and band2 is not None and ((band1[0] == 'p' and band2[0] == 'q') or (band1[0] == 'q' and band2[0] == 'p')):
                self.mutation = 'translocation'
            
            if band1 is not None and band2 is not None:
                return tools.formulate_event(self.mutation, self.align1.target, self.align2.target, band1, band2)
            else:
                return 'NA'
        else:
            return 'NA'
        
    def size(self):
        """Reports size of event"""
        if self.align1.target == self.align2.target and self.breaks:
            return abs(self.breaks[1] - self.breaks[0]) - 1    
        else:
            return '-'
            
    def has_good_spanning_support(self, minimum, debug=False):
        """Determines if event has minimum spanning read support"""
        if self.num_spanning_reads >= minimum and min(self.num_spanning_reads_forward, self.num_spanning_reads_reverse) >= 0:
            return True
        else:
            if debug:
                print 'fail_contig_read_support(%d):%s %s total:%d min_opposite:%d' % (minimum, self.align1.query, self.breakpoint, 
                                                                                   self.num_spanning_reads, min(self.num_spanning_reads_forward, self.num_spanning_reads_reverse))
            return False
        
    def has_good_pair_support(self, minimum, maximum, debug=False):
        """Determines if event has good read pair support"""
        if self.total_read_pairs() >= minimum and self.total_read_pairs() <= maximum and (self.num_flanking_pairs >= 1 or min(self.num_breakpoint_pairs) >= 1):
            return True
        else:
            if debug:
                print 'fail pair_support(%d, %d):%s %s total:%d min_breakpoint_pairs:%d' % (minimum, maximum, self.align1.query, self.breakpoint,
                                                                                        self.total_read_pairs(), min(self.num_breakpoint_pairs))
            return False
        
    def total_read_pairs(self):
        """Sums up the total of flanking and breakpoint pair support"""
        return max(0, self.num_flanking_pairs) + max(0, self.num_breakpoint_pairs[0]) + max(0, self.num_breakpoint_pairs[1])
    
    def classify(self, ff, is_genome=False):
        """Classifies event into PTD, ITD, sense or antisense fusions"""
        if self.gene1.upper() != 'NA' and self.gene2.upper() != 'NA' and self.txt1 and self.txt2:
            if self.gene1 == self.gene2:
                #ITD or PTD
                if self.mutation == 'duplication':
                    if ('exon' in self.exon1 or 'intron' in self.exon1) and ('exon' in self.exon2 or 'intron' in self.exon2):                        
                        orient1, orient2 = 'L', 'L'
                        if not self.orients[0]:
                            orient1 = 'R'
                        if not self.orients[1]:
                            orient2 = 'R'
                                                                                
                        exon1 = self.exon1.lstrip('exon')
                        if 'intron' in self.exon1:
                            exon1 = self.exon1.lstrip('intron')
                        exon2 = self.exon2.lstrip('exon')
                        if 'intron' in self.exon2:
                            exon2 = self.exon2.lstrip('intron')
                            
                        if exon1 == exon2 and 'exon' in self.exon1 and 'exon' in self.exon2:
                            self.fusion_type = 'ITD'
                            
                        #PTD must cross exon boundaries
                        if ff.cross_exon_bounds(self.breaks[0], orient1, self.txt1, exon1) and\
                           ff.cross_exon_bounds(self.breaks[1], orient2, self.txt2, exon2) and\
                           intspan.overlap([self.align1.tstart, self.align1.tend], self.txt1.get_exon_bounds(exon1)) and\
                           intspan.overlap([self.align2.tstart, self.align2.tend], self.txt2.get_exon_bounds(exon2)):
                            self.fusion_type = 'PTD'
            
            # if sample is genome, don't require event to cross exon boundaries to be called gene fusions
            elif self.exon1.upper() != 'NA' and self.exon2.upper != 'NA'\
                 and (is_genome or (self.cross_exon1 and self.cross_exon2))\
                 and not (intspan.overlap([self.txt1.txStart, self.txt1.txEnd], [self.txt2.txStart, self.txt2.txEnd]) and not (self.exon_bound1 == 'yes' or self.exon_bound2 == 'yes')):
                if self.sense1 == self.sense2:
                    self.fusion_type = 'sense_fusion'
                else:
                    self.fusion_type = 'antisense_fusion'
                    
    def pick_best_txt(self, ff, target_proper):
        """Picks best transcript that overlaps both alignments of event"""
        region = ' '.join((target_proper, str(self.breaks[0]), str(self.breaks[1])))
        features = ff.get_feature(region, gene_strand=True, txt_obj=True, all_overlaps=True)
                
        # key = num exon boundaries matched
        candidates = {2:[], 1:[], 0:[]}
        if features:
            txts_checked = {}
            
            for feature in features:
                if feature[1].alias != self.gene1:
                    continue
                
                if txts_checked.has_key(feature[1].name):
                    continue
                
                txts_checked[feature[1].name] = True
                olap_align1 = olap_align2 = False
                exon_bound1 = exon_bound2 = False
                exon1 = exon2 = None
                for i in range(len(feature[1].exons)):
                    exon = i + 1
                    exon_span = feature[1].get_exon_bounds(exon)
                    if intspan.overlap(exon_span, [self.align1.tstart, self.align1.tend]):
                        olap_align1 = True
                        exon1 = exon
                        
                        if self.at_exon_bound('exon%s' % exon, feature[1], self.breaks[0]):
                            exon_bound1 = True
                            exon1 = exon
                            
                    if intspan.overlap(exon_span, [self.align2.tstart, self.align2.tend]):
                        olap_align2 = True
                        exon2 = exon
                        
                        if self.at_exon_bound('exon%s' % exon, feature[1], self.breaks[1]):
                            exon_bound2 = True
                            exon2 = exon
                                                                                        
                if olap_align1 and olap_align2:
                    num_exon_bounds = 0
                    if exon_bound1:
                        num_exon_bounds += 1
                    if exon_bound2:
                        num_exon_bounds += 1
                        
                    candidates[num_exon_bounds].append([feature[1], [exon1, exon2]])
                    
        # go through candidates with number exon boundaries matched, 2 to 0
        for i in range(2, -1, -1):
            if not candidates[i]:
                continue
            
            if len(candidates[i]) == 1:
                self.txt1 = self.txt2 = candidates[i][0][0]
                self.exon1 = 'exon%s' % candidates[i][0][1][0]
                self.exon2 = 'exon%s' % candidates[i][0][1][1]
                break
                
            # pick longest transcript if multiple candidates exist at same number of matched exon boundaries
            else:
                self.txt1 = self.txt2 = Transcript.pick_longest([c[0] for c in candidates[i]], coding=False)
                for c in candidates[i]:
                    if c[0].name == self.txt1.name:
                        self.exon1 = 'exon%s' % c[1][0]
                        self.exon2 = 'exon%s' % c[1][1]
                        break
                break
                                                    
    def order_gene(self):
        """Defines 5' and 3' genes, transcipts, exons when possible"""
        result = None
        if self.gene1.upper() != 'NA' and self.gene2.upper() != 'NA' and self.txt1 and self.txt2:            
            self.gene5prime = '-'
            self.exon5prime = '-'
            self.gene3prime = '-'
            self.exon3prime = '-'
            
            prime1 = self.prime_gene(self.orients[0], self.txt1)
            prime2 = self.prime_gene(self.orients[1], self.txt2)

            if prime1 != prime2 and prime1 in (5,3) and prime2 in (5,3):
                if prime1 == 5:
                    # original gene order doesn't need to be flipped
                    self.genes_flipped = False
                    self.gene5prime = self.gene1
                    self.gene3prime = self.gene2
                    self.txt5prime = self.txt1
                    self.txt3prime = self.txt2
                    self.align5prime = self.align1
                    self.align3prime = self.align2
                    self.sense5prime = self.sense1
                    self.sense3prime = self.sense2
                    
                    if 'exon' in self.exon1:
                        self.exon5prime = self.exon1.lstrip('exon')
                    elif 'intron' in self.exon1:
                        self.exon5prime = self.exon1.lstrip('intron')
                    elif 'utr' in self.exon1:
                        self.exon5prime = self.exon1
                        
                    if 'exon' in self.exon2:
                        self.exon3prime = self.exon2.lstrip('exon')
                    elif 'intron' in self.exon2:
                        self.exon3prime = str(int(self.exon2.lstrip('intron')) + 1)
                    elif 'utr' in self.exon2:
                        self.exon3prime = self.exon2
                                                
                else:
                    # original gene order flipped
                    self.genes_flipped = True
                    self.gene5prime = self.gene2
                    self.gene3prime = self.gene1
                    self.txt5prime = self.txt2
                    self.txt3prime = self.txt1
                    self.align5prime = self.align2
                    self.align3prime = self.align1
                    self.sense5prime = self.sense2
                    self.sense3prime = self.sense1
                    
                    if 'exon' in self.exon2:
                        self.exon5prime = self.exon2.lstrip('exon')
                    elif 'intron' in self.exon2:
                        self.exon5prime = self.exon2.lstrip('intron')
                    elif 'utr' in self.exon2:
                        self.exon5prime = self.exon2
                        
                    if 'exon' in self.exon1:
                        self.exon3prime = self.exon1.lstrip('exon')
                    elif 'intron' in self.exon1:
                        self.exon3prime = str(int(self.exon1.lstrip('intron')) + 1)
                    elif 'utr' in self.exon1:
                        self.exon3prime = self.exon1
                    
    def prime_gene(self, orient, txt):
        """Assigns 5' and 3' to fusion genes when possible
        orient=L, txt_strand='+' => 5'
        orient=L, txt_strand='-' => 3'
        orient=R, txt_strand='+' => 3'
        orient=R, txt_strand='-' => 5'
        """
        prime = None
        if orient:
            if txt.strand == '+':
                prime = 5
            elif txt.strand == '-':
                prime = 3
        else:
            if txt.strand == '+':
                prime = 3
            elif txt.strand == '-':
                prime = 5
                
        return prime
    
    def breakpoints_at_exon_bounds(self):
        """Determines for each breakpoint whether they are at exon boundary"""
        self.exon_bound1, self.exon_bound2 = 'NA', 'NA'
        if 'exon' in self.exon1:
            if self.at_exon_bound(self.exon1, self.txt1, self.breaks[0]):
                self.exon_bound1 = 'yes'
            else:
                self.exon_bound1 = 'no'
                            
        if 'exon' in self.exon2:
            if self.at_exon_bound(self.exon2, self.txt2, self.breaks[1]):
                self.exon_bound2 = 'yes'
            else:
                self.exon_bound2 = 'no'
                    
    def at_exon_bound(self, feature, txt, coord):
        """Determines if given coordinate is at exon boundary given transcript and feature"""
        exon_bound = False
        
        # only consider if 'exon' is in feature
        if 'exon' in feature:
            if txt.strand == '+':
                exon = txt.exons[int(feature.lstrip('exon')) - 1]
            elif txt.strand == '-':
                exon = txt.exons[::-1][int(feature.lstrip('exon')) - 1]
                
            if int(coord) == int(exon[0]) or int(coord) == int(exon[1]):
                exon_bound = True
                
        return exon_bound
                
    def check_frame(self, ff, refseq, contig_seq, min_match=3):
        """Determines if gene fusion is in- or out-of-frame"""
        # both breakpoints within introns, can/should not do translation
        if 'intron' in self.exon1 and 'intron' in self.exon2:
            if self.exon5prime.isdigit() and (self.exon3prime.isdigit() or self.exon3prime == '3utr'):
                self.frame = 'in'
                
        # one breakpoint in intron and the other on exon boundary: no need for translation
        elif 'intron' in self.exon1 or 'intron' in self.exon2:
            if 'intron' in self.exon1 and 'exon' in self.exon2 and self.exon_bound2.lower() == 'yes':
                self.frame = 'in'
            elif 'intron' in self.exon2 and 'exon' in self.exon1 and self.exon_bound1.lower() == 'yes':
                self.frame = 'in'
        
        # don't need to check frame if both breakpoints land on exon boundaries
        elif self.exon5prime.isdigit() and self.exon3prime.isdigit() and self.exon_bound1.lower() == 'yes' and self.exon_bound2.lower() == 'yes':
            self.frame = 'in'
            
        if self.frame != 'in' and self.exon5prime == '5utr' and self.txt3prime.coding_type() == 'CODING':
            if (not self.genes_flipped and self.exon_bound2.lower() == 'yes') or \
               (self.genes_flipped and self.exon_bound1.lower() == 'yes'):
                self.frame = 'in'
                        
            # breakpoint not at exon boundary
            elif self.exon3prime.isdigit():
                if not self.genes_flipped:
                    target3prime = self.align2.target
                else:
                    target3prime = self.align1.target
                    
                cdna3prime = ff.construct_cdna([target3prime], self.txt3prime, refseq)
                pep3prime = translate(cdna3prime, orient=self.txt3prime.strand, frame=0)

                # find out starting amino acid of 3' transcript
                # first, find out corresponding cDNA coordinate
                exons = self.txt3prime.exons
                if self.txt3prime.strand == '-':
                    exons.reverse()

                # get transcript coordinate of exon3prime
                txt_coord = 0
                for i in range(len(exons)):
                    exon_size = int(exons[i][1]) - int(exons[i][0]) + 1
                                        
                    if i < int(self.exon3prime) - 1:
                        txt_coord += exon_size
                    elif i == int(self.exon3prime) - 1:
                        txt_coord += 1
                        break
                    
                # convert transcript coordinate to cds coordinate
                if self.txt3prime.strand == '+':
                    cds_coord = txt_coord - (int(self.txt3prime.cdsStart) - int(self.txt3prime.txStart))
                else:
                    cds_coord = txt_coord - (int(self.txt3prime.txEnd) - int(self.txt3prime.cdsEnd))
                    
                # get peptide sequence starting from exon3prime
                aa_num = int(floor(cds_coord/3))
                if pep3prime:
                    pep3prime_after_bp = pep3prime[aa_num:]
                else:
                    pep3prime_after_bp = None
                    
                # construct all the orfs using the contig sequence
                orfs = translate(contig_seq, all=True, full=True)
                                
                orf_checks = {}
                if pep3prime_after_bp:
                    for o in range(len(orfs)):
                        orf = orfs[o]
                        if (self.sense5prime == '-' and orf[2] == 1) or\
                           (self.sense5prime == '+' and orf[2] == -1):
                            continue        
                                        
                        last_match = tools.substr_search_with_consecutive_mismatches(orf[-1].upper(), pep3prime_after_bp.upper())
                        
                        if last_match is not None:
                            orf_checks[o] = (last_match[1], last_match[1] + 1)
                        
                if orfs and orf_checks:
                    # sort orfs by length of match
                    sorted_orfs = sorted(orf_checks.iteritems(), key=operator.itemgetter(1), reverse=True)

                    # require minimum 3 amino acid match
                    if sorted_orfs[0][1][1] >= min_match:
                        self.frame = 'in'
                
                
        if self.frame != 'in' and self.sense1 == self.sense2 and self.txt5prime is not None and self.txt3prime is not None:
            # get chromosome names for constructing cDNAs
            if self.txt5prime.name == self.txt1.name:
                target5prime = self.align1.target
                target3prime = self.align2.target
            else:
                target5prime = self.align2.target
                target3prime = self.align1.target
                
            # not entire contig may be aligned, just extract aligned contig sequence
            qcoords = [int(self.align1.qstart), int(self.align1.qend), int(self.align2.qstart), int(self.align2.qend)]
            qcoords.sort(key=int)
            contig_seq = contig_seq[qcoords[0]-1:qcoords[-1]]
                
            # extract 5' and 3' cDNA sequence
            cdna5prime = ff.construct_cdna([target5prime], self.txt5prime, refseq)            
            pep5prime = translate(cdna5prime, orient=self.txt5prime.strand, frame=0)
            cdna3prime = ff.construct_cdna([target3prime], self.txt3prime, refseq)            
            pep3prime = translate(cdna3prime, orient=self.txt3prime.strand, frame=0)
            
            # whether contig covers transcription start determines how to compare amino acids
            if (self.txt5prime.strand == '+' and int(self.align5prime.tstart) <= int(self.txt5prime.cdsStart)) or\
               (self.txt5prime.strand == '-' and int(self.align5prime.tend) >= int(self.txt5prime.cdsEnd)):
                contig_has_start = True
            else:
                contig_has_start = False
                            
            # construct all the orfs using the contig sequence
            orfs = translate(contig_seq, all=True, full=True)            
            
            # key = index of orf, value = length of match
            orf_checks = {}
            if pep5prime:
                for o in range(len(orfs)):
                    orf = orfs[o]
                    if (self.sense5prime == '-' and orf[2] == 1) or\
                       (self.sense5prime == '+' and orf[2] == -1):
                        continue                                
                
                    last_match = None
                    if contig_has_start:
                        last_match = tools.substr_search_with_consecutive_mismatches(orf[-1].upper(), pep5prime.upper())

                    else:                            
                        last_match = tools.substr_search_with_consecutive_mismatches(pep5prime.upper(), orf[-1].upper())
                        
                    # keep track of last matched index of orf (= size of match), index of next AA
                    if last_match is not None:
                        orf_checks[o] = (last_match[1], last_match[1] + 1)                        
                        
            else:
                sys.stdout.write("can't determine peptide sequence of 5' gene %s %s %s %s\n" % (self.align1.query, self.txt5prime.alias, self.txt5prime.name, self.txt5prime.coding_type()))
                        
            if orfs and orf_checks:
                # sort orfs by length of match
                sorted_orfs = sorted(orf_checks.iteritems(), key=operator.itemgetter(1), reverse=True)
                
                # best orf = orfs[sorted_orfs[0][0]][-1]
                # get 3' remaing of best orf, start at second amino acid to allow 1 mismatch                    
                orf3prime = orfs[sorted_orfs[0][0]][-1][sorted_orfs[0][1][1]:][1:]
                min_len = min(min_match, len(orf3prime))
                
                # if remaining of contig orf is part of 3' gene orf, then in-frame
                if orf3prime and pep3prime and pep3prime != '' and re.search(orf3prime[:min_len], pep3prime):
                    self.frame = 'in'
                    
                # fused to 3'UTR
                elif orfs[sorted_orfs[0][0]] and len(orfs[sorted_orfs[0][0]][-1]) > 1 and self.exon3prime == '3utr':
                    self.frame = 'in'
                    
                else:
                    self.frame = 'out'
                    
    def create_probe(self, half_size):
        """Creates probe sequence for event"""
        contig_coords = [int(self.align1.qstart), int(self.align1.qend), int(self.align2.qstart), int(self.align2.qend)]
        contig_coords.sort(key = int)
        if self.contig_seq:
            # there is overlap between breakpoints, the smaller of the two middle coordinates is taken as the 
            # breakpoint to construct the probe sequence
            left_coord = max(0, contig_coords[1] - half_size)
            right_coord = min(int(self.align1.query_len), contig_coords[2] - 1 + half_size)
            
            probe_seq = '%s%s' % (self.contig_seq[left_coord : contig_coords[2] - 1].lower(), 
                                  self.contig_seq[contig_coords[2] - 1 : right_coord].upper())
            
            if self.contig_reverse is not None:
                if self.contig_reverse:
                    self.probe = tools.reverse_complement(probe_seq)
                else:
                    self.probe = probe_seq
                    
            else:
                self.probe = '-'
        else:
            self.probe = '-'
        
    def get_subseqs(self):
        """Extracts sub-contig sequences aligned to the 2 regions"""
        coords1 = [int(self.align1.qstart), int(self.align1.qend)]
        coords2 = [int(self.align2.qstart), int(self.align2.qend)]
        coords1.sort(key = int)
        coords2.sort(key = int)
                        
        contig_seq1, contig_seq2 = None, None
        if self.contig_seq:
            contig_seq1 = self.contig_seq[coords1[0] - 1:coords1[1]]
            contig_seq2 = self.contig_seq[coords2[0] - 1:coords2[1]]
         
        return contig_seq1, contig_seq2
    
    def annotate(self, ff, txt_olaps, model, chrom_proper, refseq, contig_seq, cytobands=None, is_genome=False):
        """Annotates breakpoints"""
        if self.breakpoint.upper() == 'NA':
            print 'cannot annotate %s breakpoint:%s' % (self.align1.query, self.breakpoint)
            return
        
        # try to map breakpoints onto exon boundaries first
        feature1, txt1, cross_exon1 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper), self.breaks[0], self.align1.orient, [self.align1.tstart, self.align1.tend], exonbound_only=True, unique=True)  
        feature2, txt2, cross_exon2 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align2.target, chrom_proper=chrom_proper), self.breaks[1], self.align2.orient, [self.align2.tstart, self.align2.tend], exonbound_only=True, unique=True)
        
        # if not successful, try these...
        if not (txt1 and txt2 and txt1.alias != txt2.alias):        
            same_gene = None
            if self.align1.target == self.align2.target:
                chrom = tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper)
                coords = int(self.align1.tstart), int(self.align1.tend), int(self.align2.tstart), int(self.align2.tend)
                region = [chrom, min(coords), max(coords)]
                txts = [t for t in Transcript.find_overlaps(txt_olaps, [model], chrom, region[1], region[2]) if intspan.subsume([region[1], region[2]], [t.txStart, t.txEnd])]
                genes = dict((t.alias, True) for t in txts)
                if len(genes) == 1:
                    same_gene = genes.keys()[0]
            
            # try to see if exon_bound transcript or just one transcript can be mapped
            feature1, txt1, cross_exon1 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper), self.breaks[0], self.align1.orient, [self.align1.tstart, self.align1.tend], gene=same_gene, exonbound_only=True, unique=True)  
            feature2, txt2, cross_exon2 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align2.target, chrom_proper=chrom_proper), self.breaks[1], self.align2.orient, [self.align2.tstart, self.align2.tend], gene=same_gene, exonbound_only=True, unique=True)
        
            # if yes, use the above to guide mapping of the other region
            if txt1 and txt1.alias and not txt2:
                feature2, txt2, cross_exon2 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align2.target, chrom_proper=chrom_proper), self.breaks[1], self.align2.orient, [self.align2.tstart, self.align2.tend], gene=txt1.alias)
            
            elif txt2 and txt2.alias and not txt1:
                feature1, txt1, cross_exon1 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper), self.breaks[0], self.align1.orient, [self.align1.tstart, self.align1.tend], gene=txt2.alias)  
            
            elif not txt1 and not txt2:
                feature1, txt1, cross_exon1 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper), self.breaks[0], self.align1.orient, [self.align1.tstart, self.align1.tend], gene=same_gene)  
                feature2, txt2, cross_exon2 = self.annotate_breakpoint(ff, tools.proper_chrom(self.align2.target, chrom_proper=chrom_proper), self.breaks[1], self.align2.orient, [self.align2.tstart, self.align2.tend], gene=same_gene)
            
        # store up genes, transcripts, and exons
        if feature1 and txt1:
            self.txt1 = txt1
            gene, transcript, exon = feature1.rsplit(':', 2)
            self.gene1 = gene
            self.exon1 = exon.split('(')[0]
            self.cross_exon1 = cross_exon1
        
        if feature2 and txt2:
            self.txt2 = txt2
            gene, transcript, exon = feature2.rsplit(':', 2)
            self.gene2 = gene
            self.exon2 = exon.split('(')[0]
            self.cross_exon2 = cross_exon2
            
        # try to use the same transcript for the same gene
        if self.gene1 == self.gene2 and self.gene1 != 'NA' and self.txt1 and self.txt2 and self.txt1.name != self.txt2.name:
            self.pick_best_txt(ff, tools.proper_chrom(self.align1.target, chrom_proper=chrom_proper))
                    
        # determines senses for each gene
        if self.gene1 and self.txt1 and self.gene1.upper() != 'NA':
            if self.txt1.strand == self.align1.strand:
                self.sense1 = '+'
            else:
                self.sense1 = '-'                
        if self.gene2 and self.txt2 and self.gene2.upper() != 'NA':
            if self.txt2.strand == self.align2.strand:
                self.sense2 = '+'
            else:
                self.sense2 = '-'
                
        # label intron artefact
        if self.txt1 and self.txt2 and self.txt1.name == self.txt2.name and self.mutation == 'deletion':
            if ff.del_is_intron(self.breaks[0], self.breaks[1], self.txt1):
                self.artefact = True
                                
        if not self.artefact:
            self.descriptor = self.proper_name(cytobands, chrom_proper=chrom_proper)
            self.order_gene()
            self.breakpoints_at_exon_bounds()
            self.classify(ff, is_genome=is_genome)
            self.check_frame(ff, refseq, contig_seq)
            
    def annotate_breakpoint(self, ff, target, coord, strand, span, gene=None, exonbound_only=False, unique=False):
        """Annotates individual breakpoint
        Transcript preference order:
        1. exon boundary same as breakpoint
        2. alignment crosses intron-exon boundary
        3. alignment within single intron
        4. alignment within single exon
        5. alignment within UTR
        """
        best_txt, best_feature, cross_exon = None, None, False
        txt_buffer = 5000
        
        region = ' '.join((target, str(coord), str(coord)))
        features = ff.get_feature(region, gene_strand=True, txt_obj=True, strand=strand, all_overlaps=True)                
        if features:
            if gene is not None:
                features_gene = [f for f in features if f[1].alias == gene]
                if features_gene:
                    features = features_gene
                
            # make sure align is with transcript start and end bounds (+ buffer)
            features = [f for f in features if intspan.overlap(span, [int(f[1].txStart) - txt_buffer, int(f[1].txEnd) + txt_buffer])]
                
            exon_txts = [f[1] for f in features if 'exon' in f[0]]
            exons = [f[0] for f in features if 'exon' in f[0]]    
            intron_exon_txts = []            
            intron_txts = []
            for f in features:
                if 'intron' in f[0]:
                    intron_num = f[0].rsplit(':', 2)[-1].lstrip('intron').split('(')[0]
                    intron_bounds = f[1].get_intron_bounds(intron_num)
                    
                    if intspan.subsume(span, intron_bounds):
                        intron_txts.append(f[1])
                    else:
                        intron_exon_txts.append(f[1])
            
            utr_txts = [f[1] for f in features if 'utr' in f[0]]
                                        
            exonbound_txts = []
            if exon_txts and exons:
                for i in range(len(exon_txts)):
                    if self.at_exon_bound(exons[i].rsplit('|',1)[-1].split(':')[-1], exon_txts[i], coord):
                        exonbound_txts.append(exon_txts[i])                            
            if exonbound_txts:
                best_txt = Transcript.pick_longest(exonbound_txts, coding=True)
                if not best_txt:
                    best_txt = Transcript.pick_longest(exonbound_txts, coding=False)
                cross_exon = True
                    
            # if only transcripts with exon boundary matching breakpoint is desired,
            # don't need to continue further
            elif exonbound_only:
                return best_feature, best_txt, cross_exon
            
            # if only unique transcript is desired, then no need to continue further
            elif unique:
                if len(features) == 1:
                    best_feature, best_txt = features[0][0], features[0][1]
                if exonbound_txts or exon_txts or intron_exon_txts:
                    cross_exon = True
                    
                return best_feature, best_txt, cross_exon
                    
            if not best_txt and intron_exon_txts:
                best_txt = Transcript.pick_longest(intron_exon_txts, coding=True)
                cross_exon = True
            if not best_txt and intron_txts:
                best_txt = Transcript.pick_longest(intron_txts, coding=True)
            if not best_txt and exon_txts:
                best_txt = Transcript.pick_longest(exon_txts, coding=True)
                cross_exon = True
            if not best_txt and utr_txts:
                best_txt = Transcript.pick_longest(utr_txts, coding=True)
                        
            if best_txt is not None:
                best_feature = [f[0] for f in features if f[1].name == best_txt.name][0]
                                                                        
        return best_feature, best_txt, cross_exon
    
    def get_inserted_seq(self):
        inserted_seq = ''

        if int(self.align1.qstart) < int(self.align2.qstart):
            query_coords = [self.align1.qstart, self.align1.qend, self.align2.qstart, self.align2.qend]
        else:
            query_coords = [self.align2.qstart, self.align2.qend, self.align1.qstart, self.align1.qend]
            
        if int(query_coords[2]) > int(query_coords[1]):
            inserted_coords = [int(query_coords[1]) + 1, int(query_coords[2]) - 1]
            if self.contig_seq and inserted_coords[0] < inserted_coords[1] and min(inserted_coords) > 0 and max(inserted_coords) < len(self.contig_seq):
                self.inserted_seq = self.contig_seq[inserted_coords[0] - 1 : inserted_coords[1]].upper()                        
            else:
                self.inserted_seq = '' 
            
    def find_homology(self, refseq, is_genomic, feature_finder=None, max_size=200):
        """Determines microhomology of event"""
        subseqs = self.get_subseqs()        
        if not subseqs or len(subseqs) != 2 or subseqs[0] is None or subseqs[1] is None:
            return 
        
        if int(self.align1.qstart) < int(self.align2.qstart):
            query_reused = int(self.align1.qend) - int(self.align2.qstart) + 1
        else:
            query_reused = int(self.align2.qend) - int(self.align1.qstart) + 1
        
        if query_reused > 0:
            self.micro_homol1 = [self.breaks[0], self.breaks[0]]
            self.micro_homol2 = [self.breaks[1], self.breaks[1]]
            
            if self.orients[0]:
                self.micro_homol1[0] -= (query_reused - 1)
            else:
                self.micro_homol1[1] += (query_reused - 1)
                    
            if self.orients[1]:
                self.micro_homol2[0] -= (query_reused - 1)
            else:
                self.micro_homol2[1] += (query_reused - 1)
                   
            self.micro_homol_seq1 = refseq.GetSequence(self.align1.target, self.micro_homol1[0], self.micro_homol1[1])
            self.micro_homol_seq2 = refseq.GetSequence(self.align2.target, self.micro_homol2[0], self.micro_homol2[1])
        
        # consecutive query bases used - no novel sequence
        elif query_reused == 0:
            start = end = None        
            # Right
            if not self.orients[0]:
                if is_genomic or self.exon_bound1 == 'NA' or self.exon_bound1 == 'no':
                    start, end = self.breaks[0] - max_size, self.breaks[0] - 1
                                        
            # Left
            else:
                if is_genomic or self.exon_bound1 == 'NA' or self.exon_bound1 == 'no':
                    start, end = self.breaks[0] + 1, self.breaks[0] + 1 + max_size
                                
            ref1 = None
            if start is not None and end is not None:
                ref1 = refseq.GetSequence(self.align1.target, start, end)
                
            self.micro_homol1 = [self.breaks[0], self.breaks[0]]
            self.micro_homol2 = [self.breaks[1], self.breaks[1]]
                
            if ref1 is not None:                
                subseq = subseqs[1]
                if self.align1.strand == '-':
                    subseq = tools.reverse_complement(subseq)                    
                                            
                if self.orients[0]:
                    start, end = True, False
                else:
                    start, end = False, True
                common_edge = self.find_common_edge(ref1, subseq, start=start, end=end)
                
                if len(common_edge) > 0:
                    if self.orients[0]:
                        self.micro_homol1[1] += len(common_edge)
                    else:
                        self.micro_homol1[0] -= len(common_edge)
                        
                    if self.orients[1]:
                        self.micro_homol2[0] -= (len(common_edge) - 1)
                    else:
                        self.micro_homol2[1] += (len(common_edge) - 1)
                    
                    if self.align1.strand == '-':
                        common_edge = tools.reverse_complement(common_edge)
                    self.micro_homol_seq1 = common_edge
            
            
            start = end = None
            # Right
            if not self.orients[1]:
                if is_genomic or self.exon_bound2 == 'NA' or self.exon_bound2 == 'no':
                    start, end = self.breaks[1] - max_size, self.breaks[1] - 1
            # Left
            else:
                if is_genomic or self.exon_bound1 == 'NA' or self.exon_bound2 == 'no':
                    start, end = self.breaks[1] + 1, self.breaks[1] + 1 + max_size
                                
            ref2 = None
            if start is not None and end is not None:
                ref2 = refseq.GetSequence(self.align2.target, start, end)        
            if ref2 is not None:                
                subseq = subseqs[0]
                if self.align2.strand == '-':
                    subseq = tools.reverse_complement(subseq)                    
                                            
                if self.orients[1]:
                    start, end = True, False
                else:
                    start, end = False, True

                common_edge = self.find_common_edge(ref2, subseq, start=start, end=end)
                
                if len(common_edge) > 0:
                    if self.orients[1]:
                        self.micro_homol2[1] += len(common_edge)
                    else:
                        self.micro_homol2[0] -= len(common_edge)
                        
                    if self.orients[0]:
                        self.micro_homol1[0] -= (len(common_edge) - 1)
                    else:
                        self.micro_homol1[1] += (len(common_edge) - 1)        
                    
                    if self.align2.strand == '-':
                        common_edge = tools.reverse_complement(common_edge)
                    if int(self.align2.qstart) > int(self.align1.qstart):
                        self.micro_homol_seq1 = common_edge + self.micro_homol_seq1 
                    else:
                        self.micro_homol_seq1 += common_edge
            
            self.micro_homol_seq2 = self.micro_homol_seq1
            
            if self.align1.strand == '-':
                self.micro_homol_seq1 = tools.reverse_complement(self.micro_homol_seq1)
            if self.align2.strand == '-':
                self.micro_homol_seq2 = tools.reverse_complement(self.micro_homol_seq2)
            
    def find_common_edge(self, seq1, seq2, start=False, end=False):
        """Finds common bases between 2 sequences at start or end"""
        common_edge = ''
        seq_a = seq_b = None
        if start and not end:
            seq_a, seq_b = seq1, seq2
        if not start and end:
            seq_a, seq_b = seq1[::-1], seq2[::-1]
        
        if seq_a is not None and seq_b is not None:
            x = 0
            while x >= 0 and x < min(len(seq_a), len(seq_b)):
                if seq_a[x].upper() != seq_b[x].upper():
                    break
                x += 1
            
            common_edge = seq_a[:x]
            if x > 0 and not start and end:
                common_edge = common_edge[::-1]
            
        return common_edge

            
class FusionFinder:
    
    output_headers = ['id',
                      'contig',
                      'contig_size',
                      'genomic_regions',
                      'contig_regions',
                      'strands',
                      'flanking_pairs',
                      'breakpoint_pairs',
                      'spanning_reads',
                      'spanning_reads_forward',
                      'spanning_reads_reverse',
                      'rearrangement',
                      'breakpoint',
                      'size',
                      'genes',
                      'transcripts',
                      'senses',
                      'exons/introns',
                      'exon_bounds',
                      'reciprocal',
                      'descriptor',
                      'orientations',
                      "5'gene",
                      "3'gene",
                      "5'exon",
                      "3'exon",
                      'frame',
                      'probe',
                      'repeat1',
                      'repeat2',
                      'alignment_params',
                      'type',
                      'dbsnp',
                      'dgv']

    def __init__(self, params, debug=False, unique=False, is_multi_mapped=True, annodir=None, mmcfg=None):
        self.params = params
        self.fusions = []
        self.debug = debug
        self.unique = unique
        self.annodir = annodir
        self.mmcfg = mmcfg

        self.filtering = {'min_id':98, 'overlap':True, 'haplotype':True, 'readpairs':[0,0], 'min_span_reads':0, 'gene_fusions':False}        
        if params and params.has_key('min_identity'):
            self.filtering['min_id'] = params['min_identity']
        
        # debug option        
        if debug:
            print self.params
                        
        self.is_multi_mapped = is_multi_mapped
        
        self.fusion_groups = {}
        
        # FeatureFinder
        self.ff = None
        self.txt_olaps = None
        self.cytobands = None
        self.contigs_file = None        
        self.splice_motif_file = None
        self.refseq = None
        self.genome = None
        
        # for hg19 chromosome name conversion
        self.chrom_proper = None
        
        # for ORF
        self.contig_seq = None
        
    def prepare_annotation(self):
        """Prepares annotation files for overlapping"""
        if self.genome:
            # for conversion of chromosome name between FASTA and annotation (hg19)
            self.chrom_proper = tools.ucsc_chroms(self.genome, self.annodir)

            # for discerning transcript orientation in determining 5' and 3' genes            
            splice_motif_file = os.path.join(self.annodir, self.genome, 'splice_motifs.txt')
            if os.path.isfile(splice_motif_file):
                self.splice_motif_file = splice_motif_file
            else:
                print 'No splice_motifs.txt ...'
            
            # for cytoband descriptor
            cytoband_file = os.path.join(self.annodir, self.genome, 'cytoBand.txt')
            if os.path.isfile(cytoband_file):
                self.cytobands = tools.extract_cytobands(self.annodir, cytoband_file)
            else:
                print 'No cytoBand.txt ...'
                
            # for constructing cDNA sequences for in/out frame determination
            self.refseq = tools.get_refseq_from_2bit(self.annodir, self.genome)
            
            # for finding genes, exons, etc
            if self.model:
                self.ff = FeatureFinder(self.genome, self.model, self.annodir, self.mmcfg)
                self.txt_olaps = Transcript.prepare_overlap(self.genome, self.model, self.annodir, self.mmcfg)
        
    def get_split_aligns(self, aligns):
        """Extracts fusions through split-alignments"""
        fusions = []
        i = 0
        while i <= len(aligns) - 1:
            start = i
            end = i
            j = i + 1
            while j < len(aligns) and aligns[j].query == aligns[i].query:
                end = j 
                j += 1
                    
            # identifies fusions from alignments of same contig
            fusions.extend(self.extract_fusions(aligns[start:end + 1]))
            i = end + 1
                
        return fusions
    
    def extract_fusions(self, aligns):
        """Identifies fusions from alignments of same contig
        This is done by extending from 1 end of contig to the other end
        Maximum number of links: 10
        """
        fusions = []
                
        # filters out alignments below minimum percent identity
        exclude = {}
        for i in range(len(aligns)):
            if aligns[i].identity < self.params['min_identity']:
                exclude[i] = True
                if self.debug:
                    print 'contig:%s screened out because identity %s less than %s' % (aligns[0].query, aligns[i].identity, self.params['min_identity'])
                            
        # identify start of fusion partners by checking query starts
        start = None
        for i in range(len(aligns)):            
            if exclude.has_key(i):
                continue            
            if start is None:
                start = i        
            elif int(aligns[i].qstart)  < int(aligns[start].qstart):
                start = i               
            # tie in qstart
            elif aligns[i].qstart == aligns[start].qstart and aligns[i].score > aligns[start].score:
                start = i
                    
        # identify end of fusion partners by checking query ends
        end = None
        for i in range(len(aligns)):
            if exclude.has_key(i):
                continue    
            if end is None:
                end = i                
            elif int(aligns[i].qend) > int(aligns[end].qend):
                end = i                        
            #tie in qend
            elif aligns[i].qend == aligns[end].qend and aligns[i].score > aligns[end].score:
                end = i
                
        if self.debug:
            print 'contig:%s start:%s-end:%s' % (aligns[0].query, start, end)
            
        failed = False        
        if start != end and start is not None and end is not None:
            # bridge
            current = start
            used = [start]
            params = {}
            
            # check how many links
            count = 0
            while used[-1] != end:
                next, params[next] = self.plank(current, aligns, used)
                if next is not None:
                    used.append(next)
                    current = next            
                
                count += 1
                if count > 10 or (next is None and used[-1] != end):
                    failed = True
                    break
                
            if self.debug:
                for i in used:
                    print '%d %s:%s-%s %s:%s-%s' % (i, aligns[i].query, aligns[i].qstart, aligns[i].qend, aligns[i].target, aligns[i].tstart, aligns[i].tend)
                
            if not failed:
                # check if path is valid
                passed, params2 = self.valid_path(aligns, used)
                
                if not passed:
                    failed = True
                else:
                    if params2 is not None:
                        for i in used[1:]:
                            for key, value in params2.iteritems():
                                params[i][key] = value
            
            if not failed:
                # creates fusion made up of consecutive partners
                for i in range(len(used)-1):
                    align1, align2 = aligns[used[i]], aligns[used[i+1]]
                    fusion = Fusion(align1, align2)

                    # determine if contig should be flipped, for adding contig-reads and outputing probe and contig sequence
                    fusion.contig_reverse = False
                    if int(fusion.align1.qstart) > int(fusion.align2.qstart):
                        fusion.contig_reverse = True
                    
                    tools.set_attrs(fusion, params[used[i+1]])
                    
                    # filter out mitochondria ones 
                    if self.filtering.has_key('keep_mito') and not self.filtering['keep_mito'] and fusion.has_mito():
                        continue

                    fusions.append(fusion)
                                    
        return fusions
        
    def plank(self, current, aligns, used):
        """Extends bridges of fusion partners"""
        next = None
        min_distance = None
        max_score = None
        params_all = {}
        
        for i in range(len(aligns)):                
            if i in used:
                continue
            
            # not extending
            is_valid_pair, params = self.valid_pair(aligns[current], aligns[i])
            if not is_valid_pair:
                continue
            else:
                params_all[i] = params
            
            distance = abs(int(aligns[current].qend) - int(aligns[i].qstart))
                        
            if next is None:
                next = i
                min_distance = distance
                
            elif self.similar_extension(aligns[next], aligns[i]):
                if int(aligns[i].score) > int(aligns[next].score):
                    next = i
                    min_distance = distance
                continue
                                                
            elif distance < min_distance:
                if abs(distance - min_distance) > 5:
                    next = i
                    min_distance = distance
                elif int(aligns[i].score) > int(aligns[next].score):
                    next = i
                    min_distance = distance
                
        if next is None:
            return None, None
        else:
            return next, params_all[next]
    
    def similar_extension(self, align1, align2):
        """Determines if extension is really close, doesn't really extend"""
        close = 10
        if abs(int(align1.qstart) - int(align2.qstart)) <= close or abs(int(align1.qend) - int(align2.qend)) <= close:
            return True
                
    def valid_pair(self, align1, align2):
        """Determines if given alignments are valid pair for extension
        - no single alignment can cover over 95% of contig
        - one alignment cannot be subsumed in the other
        - query cannot overlap too much
        """
        passed = True
        params = {}
        
        # skip if any member of the alignment set has a single alignment over 95% of the total length
        for align in (align1, align2):
            match_len = int(align.qend) - int(align.qstart)
            if (float(match_len) / float(align.query_len) >= self.params['max_match_fraction']):
                passed = False
                if self.debug:
                    self.output_message('Failed', align1, align2, 'max_match_fraction', float(match_len) / float(align.query_len), self.params['max_match_fraction'])
                return passed, params
        
         
        # skip if one of the alignments is a subset of the other
        if intspan.subsume([align1.qstart, align1.qend], [align2.qstart, align2.qend]) or intspan.subsume([align2.qstart, align2.qend], [align1.qstart, align1.qend]):
            passed = False
            if self.debug:
                self.output_message('Failed', align1, align2, 'query_subset')
            return passed, params

        
        # check how much the query alignments between the contigs overlap
        (overlap_query, overlap_query_fraction) = self.overlap([align1.qstart, align1.qend], [align2.qstart, align2.qend], align1.query_len)
        (overlap_q1, overlap_q1_fraction) = self.overlap([align1.qstart, align1.qend], [align2.qstart, align2.qend], int(align1.qend) - int(align1.qstart) + 1)
        (overlap_q2, overlap_q2_fraction) = self.overlap([align1.qstart, align1.qend], [align2.qstart, align2.qend], int(align2.qend) - int(align2.qstart) + 1)
        params['overlap_query_fraction'] = overlap_query_fraction

        # check if the overlap query between the contigs is too high
        if not overlap_query or overlap_query <= self.params['max_query_olap'] and\
           overlap_q1_fraction <= self.params['max_query_olap_single'] and\
           overlap_q2_fraction <= self.params['max_query_olap_single']:

            # get the target query overlap
            (overlap_target, overlap_target_fraction) = self.overlap([align1.tstart, align1.tend], [align2.tstart, align2.tend], align1.query_len)
            params['overlap_target_fraction'] = overlap_target_fraction
                    
            # if targets are not the same, means no overlap
            if align1.target == align2.target:
                (overlap_target, overlap_target_fraction) = self.overlap([align1.tstart, align1.tend], [align2.tstart, align2.tend], align1.query_len)
            else:
                (overlap_target, overlap_target_fraction) = (0.0, 0.0)
                                                                
        else:
            passed = False
            if self.debug:
                if overlap_query > self.params['max_query_olap']:
                    self.output_message('Failed', align1, align2, 'max_query_olap', overlap_query, self.params['max_query_olap'])
                if overlap_q1_fraction > self.params['max_query_olap_single']:
                    self.output_message('Failed', align1, align2, 'max_query_olap_single', overlap_q1_fraction, self.params['max_query_olap_single'])
                if overlap_q2_fraction > self.params['max_query_olap_single']:
                    self.output_message('Failed', align1, align2, 'max_query_olap_single', overlap_q2_fraction, self.params['max_query_olap_single'])
                
        return passed, params
    
    def valid_path(self, aligns, path):
        """Checks if the fusion path is valid
        - query coordinates are always progressing
        - covers contig by minimum fusion fraction
        """
        passed = False
        params = {}
        
        # checking if path goes backwards
        for i in range(len(path)):
            if i > 0 and int(aligns[path[i]].qend) <= int(aligns[path[i - 1]].qend):
                return passed, params
        
        align_paths = []
        for i in path:
            align_paths.append(aligns[i])
            
        contig_cover, contig_cover_fraction = self.combined_aligned(align_paths)        
        params['contig_cover_fraction'] = contig_cover_fraction

        # check if coverage of the contig between the two alignments is too low
        if contig_cover_fraction >= self.params['min_fusion_fraction']:
            passed = True
            if not passed and self.debug:
                self.output_message('Failed', align_paths[0], align_paths[-1], 'min_fusion_fraction', contig_cover_fraction, self.params['min_fusion_fraction'])
        
        return passed, params
    
    def output_message(self, status, align1, align2, metric, value=None, threshold=None):
        """Output reason for why alignments fail to quality for event"""
        message_data = [status,
                        align1.query, 
                        align1.qstart, align1.qend, align1.target, align1.tstart, align1.tend,
                        align2.qstart, align2.qend, align2.target, align2.tstart, align2.tend,
                        metric]
        
        data = {'status': status,
                'contig': align1.query,
                'qstart1': align1.qstart,
                'qend1': align1.qend,
                'target1':align1.target,
                'tstart1': align1.tstart,
                'tend1': align1.tend,
                'qstart2': align2.qstart,
                'qend2': align2.qend,
                'target2': align2.target,
                'tstart2': align2.tstart,
                'tend2': align2.tend,
                'metric': metric
                }
                
                
        if value is None or threshold is None:
            message = "%(status)s contig:%(contig)s aln1:%(qstart1)s-%(qend1)s %(target1)s %(tstart1)s-%(tend1)s aln2:%(qstart2)s-%(qend2)s %(target2)s %(tstart2)s-%(tend2)s metric:%(metric)s" % data
        else:
            data['value'] = value
            data['threshold'] = threshold
            
            message = "%(status)s contig:%(contig)s aln1:%(qstart1)s-%(qend1)s %(target1)s %(tstart1)s-%(tend1)s aln2:%(qstart2)s-%(qend2)s %(target2)s %(tstart2)s-%(tend2)s metric:%(metric)s value:%(value)s threshold:%(threshold)s" % data
                        
        print message
        
    def overlap(self, span1, span2, total_len):
        """Calculates overlap between two ranges and gives a proportion of the total length covered"""
        s1 = [int(span1[0]), int(span1[1])]
        s2 = [int(span2[0]), int(span2[1])]
        s1.sort(lambda x,y:x-y)
        s2.sort(lambda x,y:x-y)

        olap = intspan.overlap(s1, s2)
        olap_fraction = 0.0
        if olap:
            olap_fraction = float(olap) / float(total_len)

        return (olap, olap_fraction)
        
    def non_polyA_len(self, seq):
        """Guesses non-polyA tail sequence length of contig seqeunce"""
        # check beginning of sequence
        match_start = re.search('^[at]{3,}', seq, re.IGNORECASE)
        
        # check end of sequence
        match_end = re.search('[at]{3,}$', seq, re.IGNORECASE)
        
        length = len(seq)
        if match_start:
            length -= len(match_start.group(0))
        if match_end:
            length -= len(match_end.group(0))
        
        return length
    
    def combined_aligned(self, aligns):
        """Calculates total non-redundant size of query sequence aligned"""
        aligned = []
        olaps = []
        for i in range(len(aligns) - 1):
            align1_aligned =  int(aligns[i].qend) - int(aligns[i].qstart) + 1
            align2_aligned =  int(aligns[i + 1].qend) - int(aligns[i + 1].qstart) + 1
            olap = intspan.overlap([aligns[i].qstart, aligns[i].qend], [aligns[i + 1].qstart, aligns[i + 1].qend])
            aligned.append(align1_aligned)
            olaps.append(olap)
            
            if i + 1 == len(aligns) - 1:
                aligned.append(align2_aligned)
                
        total_aligned = sum(aligned) - sum(olaps)
        if self.params['is_genome']:
            aligned_fraction = float(total_aligned) / float(aligns[0].query_len)
        # exclude polyA sequence for caluclation in case of transcriptome
        else:
            aligned_fraction = float(total_aligned) / float(self.non_polyA_len(self.contig_seq[aligns[0].query]))
        
        return total_aligned, aligned_fraction
    
    def find(self, align_file):
        """Main wrapper function to find fusion events and annotate"""     
        # extracts and filter alignments
        sys.stderr.write("parsing alignments...")
        aligns = []
        filters = {'bestn':self.params['num_aligns']}                        
        ext = os.path.splitext(align_file)[1]
        if ext == '.psl':
            aligns = psl.parse(align_file, filters, noblocks=True, noline=True)
        elif ext == '.sam':
            aligns = sam.parse(align_file, filters, noblocks=True, noline=True)
        else:
            return False
                
        # sets strand of alignments
        for a in aligns:
            if a.target_strand:
                a.strand = a.target_strand
            elif a.query_strand:
                a.strand = a.query_strand
            else:
                a.strand = None                
        sys.stderr.write("done\n")

        # identifies splitpairs
        sys.stderr.write("extracting split alignments...")
        fusions = self.get_split_aligns(aligns)
        sys.stderr.write("done\n")
        
        # only keeps unique pairs and different partners if desired
        if self.unique:
            fusion_by_contig = {}
            same_partners = {}
            for fusion in fusions:
                if fusion_by_contig.has_key(fusion.contig):
                    fusion_by_contig[fusion.contig] += 1
                else:
                    fusion_by_contig[fusion.contig] = 1

                if fusion.align1.target == fusion.align2.target:
                    same_partners[fusion.contig] = True

            unique_fusions = [fusion for fusion in fusions if fusion_by_contig[fusion.contig] == 1 and not same_partners.has_key(fusion.contig)]
            self.fusions.extend(unique_fusions)            
        else:
            self.fusions.extend(fusions)
            
        # defines mutation and breakpoint, updates splice sites
        sys.stderr.write("updating split alignments...")
        if self.splice_motif_file:
            splice_motifs = get_splice_motifs(self.splice_motif_file)
        
        for i in range(len(self.fusions)):
            fusion = self.fusions[i]
            fusion.define_mutation()         
            
            if fusion.breaks is None:
                fusion.artefact = True    
            else:
                fusion.annotate(self.ff, self.txt_olaps, self.model, self.chrom_proper, self.refseq, self.contig_seq[fusion.align1.query], cytobands=self.cytobands)
            
                for align in (fusion.align1, fusion.align2):
                    if align.blocks and self.refseq and splice_motifs:
                        align.get_splice_sites(self.refseq)
                        align.set_orient(splice_motifs)  
        sys.stderr.write("done\n")
        
        # creates probe sequence
        if self.contig_seq:
            print 'creating probe with half length %d' % self.params['probe_len_per_side']
            for fusion in self.fusions:
                if self.contig_seq.has_key(fusion.align1.query):
                    fusion.contig_seq = self.contig_seq[fusion.align1.query]
                    fusion.create_probe(half_size=self.params['probe_len_per_side'])
                    
                    # don't allow non-AGTC in breakpoint/probe sequence
                    if re.search('[^agtc]', fusion.probe, re.IGNORECASE):
                        fusion.artefact = True
                        print 'skipping event of %s because of invalid bases in probe sequence: %s' % (fusion.align1.query, fusion.probe)
                        
        # removes artefacts (e.g. intron as deletion)
        remove = []
        for i in range(len(self.fusions)):
            if self.fusions[i].artefact:
                remove.append(i)
        tools.remove_from_list(self.fusions, remove)
             
    def report(self, outdir, append=False, filtering=False, cnv_file=None, unfiltered=False):
        """Outputs results"""
        if append:
            mode = 'a'
        else:
            mode = 'w'

        unfiltered_outfile = outdir + "/fusions.tsv"
        local_outfile = outdir + "/local.tsv"
        seq_outfile = outdir + "/fusions_filtered.fa"
                
        if cnv_file:
            self.output_headers.append("cnv")
            
        if not filtering or unfiltered:
            if mode == 'w' and os.path.exists(unfiltered_outfile):
                os.remove(unfiltered_outfile)
            out = open(unfiltered_outfile, mode)
            out.write('\t'.join(self.output_headers) + '\n')
            self.output(out, self.fusions)
                                
        if filtering:                    
            remove_bp = []            
            good_contigs = {}
            reciprocals = {}
            good_fusions = []
            for bp, fusions in self.fusion_groups.iteritems():        
                for i in range(len(fusions)):        
                    if not fusions[i].exclude and fusions[i].read_support:
                        good_contigs[fusions[i].align1.query] = True
                    
                        if fusions[i].reciprocal is not None:
                            reciprocals[fusions[i].reciprocal] = True
                            
            #include all reciprocal partners of good contigs
            for bp, fusions in self.fusion_groups.iteritems():
                for fusion in fusions:
                    if reciprocals.has_key(fusion.breakpoint):
                        good_contigs[fusion.align1.query] = True
                    
            for bp, fusions in self.fusion_groups.iteritems():
                for fusion in fusions:
                    if good_contigs.has_key(fusion.align1.query):
                        good_fusions.append(fusion)
                    
                    elif not fusion.exclude and fusion.read_support:
                        good_fusions.append(fusion)
                        
            # filtered results into different types
            groups = ('sense_fusion', 'antisense_fusion', 'PTD', 'ITD', 'LSR')
            for group in groups:
                # identify fusions for specific group
                fusions = [f for f in good_fusions if f.fusion_type == group]
                        
                # re-group to give new order
                self.group_by_breakpoint(fusions, in_place=False)
                        
                outfile = '%s/%s.tsv' % (outdir, group)
                if mode == 'w' and os.path.exists(outfile):
                    os.remove(outfile)
                out = open(outfile, mode)
                out.write('\t'.join(self.output_headers) + '\n')
                self.output(out, fusions)
                out.close()
            
            # sequence file
            outfile = seq_outfile
            if mode == 'w' and os.path.exists(outfile):
                os.remove(outfile)
            out = open(outfile, mode)
            self.output_contig_seq(out, good_fusions)
            out.close()
            
            # local
            del good_fusions[:]
            for bp, fusions in self.fusion_groups.iteritems():
                for fusion in fusions:
                    if fusion.local and fusion.read_support:
                        good_fusions.append(fusion)
            
            # re-group to give new order
            self.group_by_breakpoint(good_fusions, in_place=False)
                        
            if mode == 'w' and os.path.exists(local_outfile):
                os.remove(local_outfile)
            out = open(local_outfile, mode)
            out.write('\t'.join(self.output_headers) + '\n')
            self.output(out, good_fusions)
            out.close() 
            
    def overlap_cnvs(self, fusions, cnvs):
        """Overlaps deletion events with CNV"""
        for fusion in fusions:
            if not fusion.lsr_region:
                continue
            
            if fusion.mutation != 'translocation':
                fusion.cnv = "%s:%s" % (fusion.lsr_region[0], cnvs.overlap(fusion.lsr_region[0], fusion.lsr_region[1], fusion.lsr_region[2])) 
            else:
                region1 = fusion.lsr_region[0]
                region2 = fusion.lsr_region[1]                
                cnv1 = cnvs.overlap(region1[0], region1[1], region1[2])
                cnv2 = cnvs.overlap(region2[0], region2[1], region2[2])                
                fusion.cnv = "%s:%s;%s:%s" % (region1[0], cnv1, region2[0], cnv2)
                
    def overlap_dbsnp(self, fusions_grouped_by_chrom, max_del=1000000):
        """Overlaps deletion events with dbSNP"""
        for chrom in fusions_grouped_by_chrom.keys():
            proper_chrom = tools.proper_chrom(chrom, chrom_proper=self.chrom_proper)
            snp_overlap = dbsnp.prepare_overlap(self.genome, proper_chrom, self.annodir)            
            if snp_overlap is None:
                return
            
            for event, fusions in fusions_grouped_by_chrom[chrom].iteritems():
                fusion = fusions[0]        
                if fusion.mutation == 'deletion':
                    if fusion.size() >= max_del:
                        continue
                    
                    breaks = [fusion.breaks[0], fusion.breaks[1]]
                    breaks.sort(key = int)
                    known = dbsnp.find_concordants({'type':'del', 'chrom':proper_chrom, 'start':breaks[0] + 1, 'end':breaks[1] - 1, 'size': breaks[1] - breaks[0] + 1}, snp_overlap, exact=False, target=chrom)
                
                    if known:
                        for f in fusions:
                            f.dbsnp = ','.join(known)
                    else:
                        for f in fusions:
                            f.dbsnp = '-'
                                                    
            snp_overlap.finish()
            
    def overlap_dgv(self, fusions_grouped_by_chrom):
        """Overlaps events with DGV (only for deletion and inversion)"""
        dgv_overlap = dgv.prepare_overlap(self.genome, self.annodir)
        if dgv_overlap is None:
            return
        
        for chrom in fusions_grouped_by_chrom.keys():
            proper_chrom = tools.proper_chrom(chrom, chrom_proper=self.chrom_proper)
            
            for event, fusions in fusions_grouped_by_chrom[chrom].iteritems():
                fusion = fusions[0]
            
                if fusion.mutation in ('deletion', 'inversion'):
                    proper_chrom = tools.proper_chrom(fusion.align1.target, chrom_proper=self.chrom_proper)                
                    if fusion.mutation == 'deletion':
                        var_type = 'del'
                    else:
                        var_type = 'inv'
                    
                    known = dgv.find_concordants({'type':var_type, 'chrom':proper_chrom, 'start':fusion.breaks[0] + 1, 'end':fusion.breaks[1] - 1}, dgv_overlap) 
                    if known:
                        for f in fusions:
                            f.dgv = ','.join(known)
                    else:
                        for f in fusions:
                            f.dgv = '-'
                    
        dgv_overlap.finish()
        
    def overlap_repeats(self):
        """Overlaps event breakpoints with repeatmasker or simple repeats"""
        repeat_overlaps = repeat.prepare_overlap(self.genome, self.annodir)
                
        for i in range(len(self.fusions)):
            fusion = self.fusions[i]
                        
            if not fusion.breaks:
                continue
            
            fusion.repeat1 = fusion.repeat2 = '-'                        
            chrom1 = tools.proper_chrom(fusion.align1.target, chrom_proper=self.chrom_proper)
            chrom2 = tools.proper_chrom(fusion.align2.target, chrom_proper=self.chrom_proper)            
            overlaps1 = repeat.find_overlaps({'chrom':chrom1, 'start':int(fusion.breaks[0]), 'end':int(fusion.breaks[0])}, repeat_overlaps)
            overlaps2 = repeat.find_overlaps({'chrom':chrom2, 'start':int(fusion.breaks[1]), 'end':int(fusion.breaks[1])}, repeat_overlaps)
            
            repeat1 = []
            if overlaps1:
                for repeat_type, types in overlaps1.iteritems():
                    if types:                        
                        for r in types.keys():
                            if repeat_type == 'segdup':
                                repeat1.append('%s_%s' % (repeat_type, r))
                            else:
                                repeat1.append(r)
                if repeat1:
                    fusion.repeat1 = ','.join(repeat1)
                                
            repeat2 = []
            if overlaps2:
                for repeat_type, types in overlaps2.iteritems():
                    if types:
                        for r in types.keys():
                            if repeat_type == 'segdup':
                                repeat2.append('%s_%s' % (repeat_type, r))
                            else:
                                repeat2.append(r)
                if repeat2:
                    fusion.repeat2 = ','.join(repeat2)        
                        
    def output(self, out, fusions, tsv=True):
        """Outputs given fusions in tab foramt according to pre-assigned output order"""
        for fusion in sorted(fusions, cmp=lambda x,y: x.output_order - y.output_order):
            out.write(fusion.details(tsv=tsv) + "\n")
            
    def output_contig_seq(self, out, fusions):
        """Outputs contig sequence of events"""
        for fusion in fusions:
            if not self.contig_seq.has_key(fusion.align1.query):
                continue
            
            if fusion.contig_reverse:
                out.write('>%s %s\n%s\n' % (fusion.align1.query, fusion.breakpoint, tools.reverse_complement(self.contig_seq[fusion.align1.query])))
            else:
                out.write('>%s %s\n%s\n' % (fusion.align1.query, fusion.breakpoint, self.contig_seq[fusion.align1.query]))
            
    def add_spanning_reads(self, groups):
        """Sums spanning read supports for individual event with group"""
        for group, fusions in groups.iteritems():
            total = 0
            total_forward, total_reverse = 0, 0
            total_breakpoint_pairs = [0,0]
            
            for fusion in fusions:
                total += fusion.num_spanning_reads
                total_forward += fusion.num_spanning_reads_forward
                total_reverse += fusion.num_spanning_reads_reverse
                
                total_breakpoint_pairs[0] += fusion.num_breakpoint_pairs[0]
                total_breakpoint_pairs[1] += fusion.num_breakpoint_pairs[1]
                
            for fusion in fusions:
                fusion.num_spanning_reads = total
                fusion.num_spanning_reads_forward = total_forward
                fusion.num_spanning_reads_reverse = total_reverse
                fusion.num_breakpoint_pairs = total_breakpoint_pairs
                
    def group_by_breakpoint(self, fusions=None, in_place=True):
        """Group events by breakpoint by assigning appropriate event ID and output order"""
        if fusions is None:
            fusions = self.fusions
        
        breakpoints = {}
        for fusion in fusions:
            if fusion.breaks is None:
                continue    
            if not breakpoints.has_key(fusion.breakpoint):
                breakpoints[fusion.breakpoint] = [fusion]
            else:
                breakpoints[fusion.breakpoint].append(fusion)
                
        reciprocal_breakpoints = [breakpoint for breakpoint in breakpoints.keys() if breakpoints[breakpoint][0].reciprocal is not None]
        
        event_id = 1
        event_ids = {}
        output_orders = {}
        order = 0        
        for breakpoint in reciprocal_breakpoints:
            if event_ids.has_key(breakpoint):
                if len(breakpoints[breakpoint]) == 1:
                    fusion = breakpoints[breakpoint][0]
                    fusion.event_id = event_ids[breakpoint]
                    fusion.output_order = output_orders[fusion.event_id]
                    
                else:
                    for i in range(len(breakpoints[breakpoint])):
                        fusion = breakpoints[breakpoint][i]
                        output_order = output_orders[event_ids[breakpoint].split('.')[0]]
                        fusion.event_id = "%s.%d" % (event_ids[breakpoint], i+1)
                        fusion.output_order = output_order
                        output_order += 1
            
            else:
                if len(breakpoints[breakpoint]) == 1:
                    fusion = breakpoints[breakpoint][0]
                    fusion.event_id = "%da" % (event_id)
                    fusion.output_order = order
                    event_ids[fusion.reciprocal] = "%db" % (event_id)
                    order += 1
                    output_orders[event_ids[fusion.reciprocal]] = order
                    order += 1
                    event_id += 1
                    
                else:
                    for i in range(len(breakpoints[breakpoint])):
                        fusion = breakpoints[breakpoint][i]
                        fusion.event_id = "%da.%d" % (event_id, i+1)
                        fusion.output_order = order
                        order += 1
                        
                    event_ids[breakpoints[breakpoint][0].reciprocal] = "%db" % (event_id)
                    output_orders[event_ids[breakpoints[breakpoint][0].reciprocal]] = order
                    event_id += 1
                    order += 1
                    
        for breakpoint in sorted(breakpoints.keys()):
            if breakpoint in reciprocal_breakpoints:
                continue
            if len(breakpoints[breakpoint]) == 1:
                fusion = breakpoints[breakpoint][0]
                fusion.event_id = "%d" % (event_id)
                fusion.output_order = order
                order += 1
                event_id += 1
                    
            else:
                for i in range(len(breakpoints[breakpoint])):
                    fusion = breakpoints[breakpoint][i]
                    fusion.event_id = "%d.%d" % (event_id, i+1)
                    fusion.output_order = order
                    order += 1
                    
                event_id += 1
                                
        # replace objects's fusion groups with grouped results  
        if in_place:
            self.fusion_groups = breakpoints

        return breakpoints
        
    def group_by_chrom(self, mutations=None):
        """Group events by chromosome for overlapping with annotations"""
        fusions_by_chrom = {}
        for breakpoint, fusions in self.fusion_groups.iteritems():
            if fusions[0].align1.target != fusions[0].align2.target:
                continue
            
            chrom = fusions[0].align1.target
            key = fusions[0].mutation + fusions[0].breakpoint
            
            if not fusions_by_chrom.has_key(chrom):
                fusions_by_chrom[chrom] = {}
                fusions_by_chrom[chrom][key] = []
            elif not fusions_by_chrom[chrom].has_key(key):
                fusions_by_chrom[chrom][key] = []
            fusions_by_chrom[chrom][key].extend(fusions)
            
        return fusions_by_chrom
                    
    def parse_fusions_dir(self, path, out=True, outdir=None):
        """Parse directory of fusions output directories
        and checks if there is no partial files
        """
        fusion_dirs = sorted(glob.glob(os.path.join(path, '*')))        
        fusion_files = []
        num_fusion_dirs = 0
        missing_dirs = []
        
        ok = True
        
        for job_num in range(1, len(fusion_dirs)+1):
            cluster_outdir = "%s/%s" % (path, job_num)
            if os.path.isdir(cluster_outdir):
                num_fusion_dirs += 1
                    
                fusion_file = cluster_outdir + '/fusions.tsv'                    
                if os.path.exists(fusion_file):
                    fusion_files.append(fusion_file)
                else:
                    missing_dirs.append(str(job_num))
                    sys.stdout.write("%s does not have fusion file\n" % (cluster_outdir))
                    ok = False
                        
        sys.stdout.write("fusion dirs:%s fusion files:%s\n" % (num_fusion_dirs, len(fusion_files)))
                
        # checks if every directory has fusion file
        if num_fusion_dirs == len(fusion_files):
            if out and outdir is not None:
                concat_outfile = outdir + "/fusions_concat.tsv"
                concat_out = open(concat_outfile, 'w')
                    
            #just for check if header is to be written
            count = 0
            for fusion_file in fusion_files:
                infile = open(fusion_file, 'r')
                if count == 0:
                    if out:
                        concat_out.writelines(open(fusion_file, 'r').readlines())
                else:
                    if out:
                        concat_out.writelines(open(fusion_file, 'r').readlines()[1:])
                infile.close()
                        
                print fusion_file
                self.parse_fusions(fusion_file)
                count += 1
                
            if out:
                concat_out.close()
        else:
            sys.stdout.write("missing (%s):%s\n" % (len(missing_dirs), ','.join(missing_dirs)))
            ok = False
        
        return ok
                        
    def parse_fusions(self, fusions_file):
        """Parses individual fusions file into objects"""    
        # preserve original output order
        output_order = 0
        for line in open(fusions_file, 'r'):
            cols = line.rstrip('\n').split('\t')
            if cols[0] == 'id':
                continue
            
            #empty line
            if len(cols) < 2:
                continue
            
            align1 = Alignment(None)
            align2 = Alignment(None)
            
            align1.query = align2.query = cols[1]
            align1.query_len = align2.query_len = cols[2]

            targets = cols[3].split(',')
            align1.target, align1.tstart, align1.tend = re.split('[:-]', targets[0])
            align2.target, align2.tstart, align2.tend = re.split('[:-]', targets[1])

            queries = cols[4].split(',')
            align1.qstart, align1.qend = re.split('-', queries[0])
            align2.qstart, align2.qend = re.split('-', queries[1])

            align1.strand, align2.strand = cols[5].split(',')

            fusion = Fusion(align1, align2)
            fusion.event_id = cols[0]
            fusion.output_order = output_order

            fusion.num_flanking_pairs = int(cols[6])
            fusion.num_breakpoint_pairs = [int(n) for n in cols[7].split(',')]
            fusion.num_read_pairs = fusion.total_read_pairs()
            fusion.num_spanning_reads = int(cols[8])
            fusion.num_spanning_reads_forward = int(cols[9])
            fusion.num_spanning_reads_reverse = int(cols[10])

            fusion.mutation = cols[11]
            fusion.breakpoint = cols[12]
            if fusion.breakpoint != 'NA':
                fusion.breaks = [int(fusion.breakpoint.split('|')[0].split(':')[1]), int(fusion.breakpoint.split('|')[1].split(':')[1])]
            
            fusion.gene1, fusion.gene2 = cols[14].split(',')
            fusion.txt1_name, fusion.txt2_name = cols[15].split(',')
            fusion.sense1, fusion.sense2 = cols[16].split(',')
            fusion.exon1, fusion.exon2 = cols[17].split(',')
            fusion.exon_bound1, fusion.exon_bound2 = cols[18].split(',')            
            
            if cols[19].upper() == 'NA':
                fusion.reciprocal = None
            else:
                fusion.reciprocal = cols[19]
            fusion.descriptor = cols[20]
                
            fusion.orients = [None, None]
            orient1, orient2 = cols[21].split(',')
            if orient1 == 'L':
                fusion.orients[0] = True
            else:
                fusion.orients[0] = False
            if orient2 == 'L':
                fusion.orients[1] = True
            else:
                fusion.orients[1] = False
            
            fusion.gene5prime = cols[22]
            fusion.gene3prime = cols[23]
            fusion.exon5prime = cols[24]
            fusion.exon3prime = cols[25]
            fusion.frame = cols[26]
            fusion.probe = cols[27]
            fusion.repeat1 = cols[28]
            fusion.repeat2 = cols[29]
            
            params = cols[30].split(',')
            fusion.overlap_target_fraction = float(params[0].split(':')[1])
            fusion.overlap_query_fraction = float(params[1].split(':')[1])
            fusion.contig_cover_fraction = float(params[2].split(':')[1])
            align1.identity = float(params[3].split(':')[1])
            align2.identity = float(params[4].split(':')[1])
                        
            fusion.fusion_type = cols[31]
            fusion.dbsnp = cols[32]
            fusion.dgv = cols[33]
            
            self.fusions.append(fusion)
            
            output_order += 1
            
    def extract_contig_seq(self, contigs_file):
        """Extracts contig sequence file into dictionary"""
        assembly = Assembly(None, fasta=contigs_file)
        contigs = dict((fusion.align1.query, True) for fusion in self.fusions)
        self.contig_seq = dict((contig.num, contig.sequence) for contig in assembly.get_contigs(ids=contigs.keys(), sequence=True))
    
    def annotate(self):
        """Wrapper function for annotating all fusion events"""
        if self.ff:
            for fusion in self.fusions: 
                fusion.annotate(self.ff, self.txt_olaps, self.model, self.chrom_proper, self.refseq, self.contig_seq[fusion.align1.query], 
                                cytobands=self.cytobands, 
                                is_genome=self.params['is_genome'])
                    
    def set_local(self):
        """Labels events as local
        All events derived from a contig will be deemed local
        where one them is
        """
        # group events by contig
        fusions_by_contig = {}
        for i in range(len(self.fusions)):
            if not fusions_by_contig.has_key(self.fusions[i].contig):
                fusions_by_contig[self.fusions[i].contig] = [i]
            else:
                fusions_by_contig[self.fusions[i].contig].append(i)

        for contig, fusion_idx in fusions_by_contig.iteritems():
            contig_local = False
            for idx in fusion_idx:
                fusion = self.fusions[idx]
                
                if not fusion.local:
                    if self.ff and fusion.is_local(self.ff.model, self.ff.annot_file):
                        fusion.local = True
                        contig_local = True
                        break
                    elif fusion.is_local():
                        fusion.local = True
                        contig_local = True
                        break
                                        
            if contig_local:
                for idx in fusion_idx:
                    self.fusions[idx].contig_local = True
                    self.fusions[idx].exclude = True
            
    def post_filter(self):
        """Filters events by read support, overlapping results with mitochondria, repeats"""
        exclude_indices = {}
        mappings = {}
        contig2fusions = {}
                    
        for i in range(len(self.fusions)):
            fusion = self.fusions[i]

            if fusion.contig_local or fusion.local:
                fusion.exclude = True

            if not fusion.exclude and self.filtering.has_key('keep_mito') and not self.filtering['keep_mito']:
                if fusion.has_mito():
                    fusion.exclude = True
                    
            if not fusion.exclude and self.filtering.has_key('no_repeats') and self.filtering['no_repeats']:
                if (fusion.repeat1 != '-' and fusion.repeat1.upper() != 'NA') or (fusion.repeat2 != '-' and fusion.repeat2.upper() != 'NA'):
                    fusion.exclude = True
                    
            if fusion.has_good_pair_support(self.filtering['readpairs'][0], self.filtering['readpairs'][1], debug=self.debug):
                fusion.good_read_pair = True                        

            if fusion.has_good_spanning_support(self.filtering['min_span_reads'], debug=self.debug):
                fusion.good_span_read = True

            if fusion.good_read_pair and fusion.good_span_read:
                fusion.read_support = True
                        
    def post_process(self, options):
        """Links reciprocal events, add spanning read support, overlaps repeats, cnvs, dbSNP, DGV"""
        self.find_reciprocal() 
            
        if options.add_spanning_reads:
            self.add_spanning_reads(self.fusion_groups)
                                                                           
        if options.cnv_file:
            cnvs = CNV(options.cnv_file)
            cnvs.extract()
            self.overlap_cnvs(self.fusions, cnvs)
            
        if options.olap_annot:
            self.overlap_repeats()
            chrom_groups = self.group_by_chrom()
            self.overlap_dbsnp(chrom_groups)
            self.overlap_dgv(chrom_groups)
            
    def read_support(self, genome_bamfile=None, contigs_bamfile=None, 
                     min_mapq=None, keep_mito=False, genome_coord_buffer=None, 
                     only_unique=True):
        """Finding read support for events"""
        sys.stderr.write("finding read support...")
        
        genome_bam = None
        contigs_bam = None
        if genome_bamfile and os.path.exists(genome_bamfile):
            genome_bam = BAM(genome_bamfile, min_mapq=min_mapq)
            
        if contigs_bamfile and os.path.exists(contigs_bamfile):
            contigs_bam = BAM(contigs_bamfile, min_mapq=options.mapq_contigs)
        
        self.fusions.sort(self.compare_breakpoints)
        for fusion in self.fusions:
            if not fusion.breaks:
                continue
            
            spanning_reads = []
            unique = None
            if contigs_bam:
                spanning_reads = fusion.find_spanning_reads(contigs_bam, 
                                                            breakpoint_buffer=self.params['min_ctg_bp_overlap'], 
                                                            max_depth=self.params['max_span_reads'], 
                                                            only_unique=only_unique, 
                                                            debug=self.debug)
                                    
            if genome_bam:
                missing_mates = None
                if spanning_reads:
                    missing_mates = {fusion.breaks[0]:{}, fusion.breaks[1]:{}}
                    for read in spanning_reads:
                        if int(fusion.align2.qend) > int(fusion.align1.qend):
                            if not read.is_reverse:
                                breakpoint = fusion.breaks[1]
                            else:
                                breakpoint = fusion.breaks[0]
                        
                        else:
                            if not read.is_reverse:
                                breakpoint = fusion.breaks[0]
                            else:
                                breakpoint = fusion.breaks[1] 
                                            
                        # don't allow reads ending in /1 or /2
                        read_name = read.qname
                        if re.search('/[12]$', read.qname):
                            read_name = read_name[:-2]
                        
                        if not missing_mates[fusion.breaks[0]].has_key(read_name):
                            missing_mates[fusion.breaks[0]][read_name] = []
                        missing_mates[fusion.breaks[0]][read_name].append(read)
                        if not missing_mates[fusion.breaks[1]].has_key(read_name):
                            missing_mates[fusion.breaks[1]][read_name] = []
                        missing_mates[fusion.breaks[1]][read_name].append(read)
                                            
                fusion.find_read_pairs(genome_bam, genome_coord_buffer=genome_coord_buffer, breakpoint_buffer=self.params['min_genome_bp_overlap'], keep_mito=keep_mito, reads=missing_mates, only_unique=only_unique, debug=self.debug)
                                                                
        sys.stderr.write("done\n")  
        
    def compare_breakpoints(self, f1, f2):
        """For ordering events by breakpoint"""
        if f1.breakpoint == 'NA' or f2.breakpoint == 'NA':
            if f1.breakpoint == 'NA' and f2.breakpoint == 'NA':
                return 0        
            elif f1.breakpoint == 'NA' and f2.breakpoint != 'NA':
                return 1
            else:
                return -1
        elif f1.align1.target < f2.align1.target:
            return -1
        elif f1.align1.target > f2.align1.target:
            return 1
        elif f1.align2.target < f2.align2.target:
            return -1
        elif f1.align2.target > f2.align2.target:
            return 1
        elif f1.breaks[0] < f2.breaks[0]:
            return -1
        elif f1.breaks[0] > f2.breaks[0]:
            return 1
        elif f1.breaks[1] < f2.breaks[1]:
            return -1
        elif f1.breaks[1] > f2.breaks[1]:
            return 1
        else:
            return 0
                
    def realign(self, cfg_file, aligner, outdir):
        """Realigns contig sequence using BLAT to filter out nonspecific mapping
        
        3 types of sequences will be generated:
        1. subsequence that aligns to first region of fusion
        2. subsequence that aligns to second region of fusion
        3. probe sequence around breakpoint
        
        Sub-sequences are used for checking if individually the mapping of the contig
        sequence is specific.  Probe is used for checking if actually the breakpoint
        can be aligned to a single location.
        
        Probe sequence will have a size limit on both sides of the breakpoint (default: 50).
        This is to reduce time for the realignment process.  The capped-off subsequences
        should be extracted from the breakpoint (middle) of the contig sequence.
        
        Sequences of all events are written to a single output FASTA file (realign.fa),
        and BLAT is run locally using parameters described in the config file.  Output's
        filename is always 'realign.psl'.
        
        Args:
            cfg_file: path of TA config file for extracting BLAT parameters
            aligner: aligner ('should be' BLAT)
            outdir: path of fusion output file where FASTA file and output are dumped
        """
        # prepares sequences for realignment
        seq_file = '%s/realign.fa' % outdir
        tmp_seq = open(seq_file, 'w')
                
        for fusion in self.fusions:
            if fusion.artefact:
                continue
                
            # subsequences, maximum length = max_len
            if fusion.contig_seq is None and self.contig_seq.has_key(fusion.align1.query):
                fusion.contig_seq = self.contig_seq[fusion.align1.query]

            subseqs = fusion.get_subseqs()
            if subseqs[0] is not None and subseqs[1] is not None:
                if self.params.has_key('max_subseq_len') and self.params['max_subseq_len'] is not None:
                    if int(fusion.align1.qstart) < int(fusion.align2.qstart):
                        subseq1 = subseqs[0][-1 * min(self.params['max_subseq_len'], len(subseqs[0]))::]
                        subseq2 = subseqs[1][:min(self.params['max_subseq_len'], len(subseqs[1]))]
                    else:
                        subseq1 = subseqs[0][:min(self.params['max_subseq_len'], len(subseqs[0]))]
                        subseq2 = subseqs[1][-1 * min(self.params['max_subseq_len'], len(subseqs[1]))::]                        
                    tmp_seq.write('>%s-1 %d\n%s\n' % (fusion.align1.query, len(subseq1), subseq1))
                    tmp_seq.write('>%s-2 %d\n%s\n' % (fusion.align2.query, len(subseq2), subseq2))
                else:
                    tmp_seq.write('>%s-1 %d\n%s\n' % (fusion.align1.query, len(subseqs[0]), subseqs[0]))
                    tmp_seq.write('>%s-2 %d\n%s\n' % (fusion.align1.query, len(subseqs[1]), subseqs[1]))
                    
            else:
                print 'failed to get subseq', fusion.align1.query
                    
            # probe sequence
            if fusion.probe != '-' and fusion.probe != 'NA':
                tmp_seq.write('>%s %d\n%s\n' % (fusion.align1.query, len(fusion.probe), fusion.probe))
                
        tmp_seq.close()
        
        # realign running Blat locally
        settings = cfg.initialize_settings(cfg_file)
        target = settings['genomes'][self.genome]        
        
        args = {'TARGET': target,
                'QUERY': seq_file,
                'OUTPUT': '%s/realign.psl' % outdir
                }
        if settings is not None:
            cmd = ' '.join([aligner, cfg.create_args(settings, 'commands', aligner, args)])
            print cmd
            os.system(cmd)
            
    def filter_nonspecific(self, outdir):
        """Using BLAT realignments, filters out non-specific fusion calls.
        
        Given the fusion output directory, expects the BLAT realignment file to be there,
        extracts all the alignments without filtering.
        For any given contig, the sub-sequence realignments will have queries
        <CTG>-1, <CTG>-2, and the probe will be <CTG>.
        
        To pass:
        1. all 3 sequences (2 subsequences + probe) have realignment results
        2. the 2 subsequences will realign to the intended regions
        3. the 2 subsequences cannot be mapped to a 'better' location
        4. the probe is not mapped to a single location
        5. the probe is aligned to 2 distinct locations
        
        It will remove the fusion events that don't fulfill the above criteria
        
        Args:
            outdir: path output fusion directory, where 'realign.pl' is expected
                    to be present
        """
        aligns = psl.parse('%s/realign.psl' % outdir, {}, noline=True)
        
        # associate every alignments to generic contig name
        contig2aligns = {}
        for align in aligns:
            contig = align.query
            if re.search('-1$', align.query) or re.search('-2$', align.query):
                contig = align.query[:-2]
                
            if not contig2aligns.has_key(contig):
                contig2aligns[contig] = []
            contig2aligns[contig].append(align)
          
        bad = {}
        for i in range(len(self.fusions)):
            if self.fusions[i].artefact:
                continue
            
            fusion = self.fusions[i]
            probe = fusion.align1.query            
            subseq1 = '%s-1' % probe
            subseq2 = '%s-2' % probe
            if not contig2aligns.has_key(probe):
                print 'failed realign absent all re-alignments %s' % probe
                bad[i] = True
                continue
            
            probe_aligns = [a for a in contig2aligns[probe] if a.query == probe]
            subseq1_aligns = [a for a in contig2aligns[probe] if a.query == subseq1]
            subseq2_aligns = [a for a in contig2aligns[probe] if a.query == subseq2]
                            
            # check for subseq1
            matched = []
            for align in subseq1_aligns:
                if self.is_mapped_to_same_region(align, fusion.align1, fusion.breaks[0]):
                    matched.append(align)
            matched.sort(key=lambda align: int(align.score), reverse=True)
            if not matched:
                print 'failed realign different-mapping subseq1 %s %s:%s-%s' %\
                      (fusion.align1.query, fusion.align1.target, fusion.align1.tstart, fusion.align1.tend)
                bad[i] = True
            elif len(subseq1_aligns) > 1:
                better = self.find_better_align(subseq1_aligns, matched[0])
                if better is not None:
                    print 'failed realign finds better mapping subseq1 %s %s:%s-%s %s:%s-%s %d' % (fusion.align1.query, fusion.align1.target, fusion.align1.tstart, fusion.align1.tend, better.target, better.tstart, better.tend, better.identity)
                    bad[i] = True
            if bad.has_key(i):
                continue
            
            # check for subseq2
            matched = []
            for align in subseq2_aligns:
                if self.is_mapped_to_same_region(align, fusion.align2, fusion.breaks[1]):
                    matched.append(align)
            matched.sort(key=lambda align: int(align.score), reverse=True)
            if not matched:
                print 'failed realign different-mapping subseq2 %s %s:%s-%s' %\
                      (fusion.align2.query, fusion.align2.target, fusion.align2.tstart, fusion.align2.tend)
                bad[i] = True
            elif len(subseq2_aligns) > 1:
                better = self.find_better_align(subseq2_aligns, matched[0])
                if better is not None:
                    print 'failed realign finds better mapping subseq2 %s %s:%s-%s %s:%s-%s %d' % (fusion.align2.query, fusion.align2.target, fusion.align2.tstart, fusion.align2.tend, better.target, better.tstart, better.tend, better.identity)
                    bad[i] = True
            if bad.has_key(i):
                continue
            
            # check if probe can be mapped to single position
            for align in probe_aligns:
                if self.check_single_align_possible(align, fusion):
                    print 'failed realign probe aligned to single position %s %s:%s-%s %s %s-%s %s' %\
                          (align.query, align.target, align.tstart, align.tend, 
                           align.query_len, align.qstart, align.qend, align.identity)
                    bad[i] = True
                    break
            if bad.has_key(i):
                continue
            
            # check if probe alignments are specific
            # assign probe alignments to each region
            probe1_aligns = []
            probe2_aligns = []
            for align in probe_aligns:
                if abs(int(align.qstart) - 1) < abs(int(align.query_len) - int(align.qend)):
                    if int(fusion.align1.qstart) < int(fusion.align2.qstart):
                        probe1_aligns.append(align)
                    else:
                        probe2_aligns.append(align)
                else:
                    if int(fusion.align1.qstart) < int(fusion.align2.qstart):
                        probe2_aligns.append(align)
                    else:
                        probe1_aligns.append(align)
                        
            # cannot align both parts of the probe
            if not probe1_aligns:
                print 'failed realign probe1 %s %s:%s-%s' % (fusion.align1.query, fusion.align1.target, fusion.align1.tstart, fusion.align1.tend)
                bad[i] = True
                continue
            if not probe2_aligns:
                print 'failed realign probe2 %s %s:%s-%s' % (fusion.align2.query, fusion.align2.target, fusion.align2.tstart, fusion.align2.tend)
                bad[i] = True
                continue
                                
        # remove events with non-specific alignments
        tools.remove_from_list(self.fusions, bad.keys())
                        
    def find_better_align(self, aligns, olap_align):
        """Checks if there is a better mapped location given by Blat realignment than 
        given Blat realignment that overlaps the original alignment.
        
        Because it's oftern hard to compare BLAT against GMAP alignments, the best is, 
        to first find out the BLAT alignment that's "the same" as the original GMAP alignment, 
        then uses that BLAT alignment to see if a better BLAT alignment exists among the 
        realignments.  Therefore, the BLAT alignment that overlaps the original GMAP alignment 
        (olap_align) has to be identified first (by is_mapped_to_same_region())
        
        A 'better' alignment is:
        1. match length of query is at least 80% of 'olap_align'
        2. number of blocks cannot be greater than that of 'olap_align'
        3a. if identity of 'olap_align' is 100%, then identity of 'align' must be also
            100% and match length is the same as 'olap_align'
        or
        3b. if identity of 'olap_align' < 100%, number of mismatches of 'align' cannot be
            more than number of mismatches of 'olap_align' by 2
            
        Args:
            aligns: all BLAT realignment objects of same sequence of 'olap_align'
            olap_align: one of the BLAT realignments that is determined to represent
                        the original fusion alignment
                        
        Returns:
            alignment object that is 'better' than 'olap_align', 
            or None (if such alignment cannot be found)
        """
        aligns_sorted = sorted(aligns, key=lambda align: int(align.score), reverse=True)
        olap_align_match_len = int(olap_align.qend) - int(olap_align.qstart) + 1
        for align in aligns_sorted:
            if len(align.blocks) > len(olap_align.blocks):
                continue
            # same alignment
            if align.target == olap_align.target and\
               align.tstart == olap_align.tstart and\
               align.tend == olap_align.tend:
                continue
            match_len = int(align.qend) - int(align.qstart) + 1
            if match_len >= olap_align_match_len * 0.8 and\
               ((olap_align.identity == 100 and match_len == olap_align_match_len and align.identity == 100) or\
                (olap_align.identity < 100 and int(align.mismatch) - int(olap_align.mismatch) <= 2)
                ):
                return align

        return None
                
    def is_mapped_to_same_region(self, align, fusion_align, breakpoint, min_match=0.8):
        """Checks if a re-alignment is mapped to the same neighborhood of the original fusion alignment

        Criteria:
        1. mapped to the same chromosome
        2. breakpoint lies in the span of the alignment
        3. overlap between realignment is at least 80% of the original fusion target span, or
           of the realignment target span
           
        Args:
            align: BLAT realignment object for comparing
            fusion_align: original fusion alignment object for comparing
            breakpoint: coordinate of breakpoint in 'fusion_align'
            min_match(optional): minimum percentage overlap between alignment regions (Default: 80%)
            
        Returns:
            True: same
            False: different
        """
        target_olap = intspan.overlap([align.tstart, align.tend], [fusion_align.tstart, fusion_align.tend])
        target_span = int(align.tend) - int(align.tstart) + 1
        fusion_target_span = int(fusion_align.tend) - int(fusion_align.tstart) + 1
        if align.target == fusion_align.target and\
           int(breakpoint) >= int(align.tstart) and int(breakpoint) <= align.tend and\
           (float(target_olap)/float(target_span) >= min_match or\
            float(target_olap)/float(fusion_target_span) >= min_match
            ):
            same = True
        else:
            same = False
    
        return same
    
    def check_single_align_possible(self, align, fusion):
        """Checks if a probe can be aligned to a single location instead of 2

        Criteria:
        1. qstart <= 2
        2. qend >= query_len - 1
        3. identify >= 90%
        
        Exceptions:
        1. event is a deletion, and BLAT will likely report a single alignment
        2. alignment overlaps with both regions of the fusion alignments
        
        Args:
            align: alignment object (of the BLAT realignment)
            fusion: fusion event object
            
        Returns:
            True: possible
            False: impossible
        """
        possible = False
        if int(align.qstart) <= 2 and int(align.qend) >= int(align.query_len) - 1 and align.identity >= 90.0:
            # unless event is deletion, alignment should not be single alignment
            if fusion.mutation != 'deletion':
                possible = True
                    
            # align overlaps with regions
            elif not (align.target == fusion.align1.target and\
                      intspan.overlap([align.tstart, align.tend], [fusion.align1.tstart, fusion.align1.tend])) or\
                 not (align.target == fusion.align2.target and\
                      intspan.overlap([align.tstart, align.tend], [fusion.align2.tstart, fusion.align2.tend])):
                possible = True
                
        return possible
            
    def find_reciprocal(self, close=100):
        """Links potential reciprocal translocation events"""
        breakpoints = {}        
        for i in range(len(self.fusions)):
            fusion = self.fusions[i]
            # skip if not translocation or breakpoint not available
            if fusion.mutation != 'translocation' or fusion.breakpoint == 'NA' or not fusion.breaks:
                continue
                     
            # group by breakpoints
            if not breakpoints.has_key(fusion.breakpoint):
                breakpoints[fusion.breakpoint] = {'target1': fusion.align1.target, 'target2': fusion.align2.target,
                                                  'break1': fusion.breaks[0], 'break2': fusion.breaks[1],
                                                  'gene1': fusion.gene1, 'gene2': fusion.gene2,
                                                  'orient1': fusion.orients[0], 'orient2': fusion.orients[1],
                                                  'indices': []}
            breakpoints[fusion.breakpoint]['indices'].append(i)
                        
        pairs = {}
        bps = breakpoints.keys()
        for i in range(len(bps) - 1):
            for j in range(i + 1, len(bps)):
                # skip if chromosomes not the same between pairs
                if breakpoints[bps[i]]['target1'] != breakpoints[bps[j]]['target1'] or \
                   breakpoints[bps[i]]['target2'] != breakpoints[bps[j]]['target2']:
                    continue
                
                # Criteria for calling reciprocal:
                # 1. opposite orientations (L,R) 
                # 2. reasonably close or within same genes
                if breakpoints[bps[i]]['orient1'] != breakpoints[bps[j]]['orient1'] and breakpoints[bps[i]]['orient2'] != breakpoints[bps[j]]['orient2'] and \
                   ((abs(breakpoints[bps[i]]['break1'] - breakpoints[bps[j]]['break1']) <= close and abs(breakpoints[bps[i]]['break2'] - breakpoints[bps[j]]['break2']) <= close) or\
                    (breakpoints[bps[i]]['gene1'] != 'NA' and breakpoints[bps[i]]['gene2'] != 'NA' and \
                     breakpoints[bps[i]]['gene1'] == breakpoints[bps[j]]['gene1'] and breakpoints[bps[i]]['gene2'] == breakpoints[bps[j]]['gene2'])):
                    
                    if not pairs.has_key(bps[i]):
                        pairs[bps[i]] = {}
                    pairs[bps[i]][bps[j]] = True
        
        # Establish reciprocal links
        for bp1 in pairs.keys():
            if len(pairs[bp1].keys()) == 1:
                bp2 = pairs[bp1].keys()[0]
                
                for i in breakpoints[bp1]['indices']:
                    self.fusions[i].reciprocal = bp2
                for i in breakpoints[bp2]['indices']:
                    self.fusions[i].reciprocal = bp1

def main(args, options):
    # set up parameters
    params = {}
    params['max_match_fraction'] = options.match_fraction
    params['max_query_olap'] = options.query_olap
    params['max_query_olap_single'] = options.query_olap_single
    params['min_fusion_fraction'] = options.fusion_fraction
    params['min_identity'] = options.min_identity
    params['num_aligns'] = options.num_aligns
    params['num_best'] = options.num_best
    params['min_ctg_bp_overlap'] = options.min_ctg_bp_overlap
    params['min_genome_bp_overlap'] = options.min_genome_bp_overlap
    params['is_genome'] = options.is_genome
    params['max_span_reads'] = options.max_span_reads
    params['max_subseq_len'] = options.max_subseq_len
    params['probe_len_per_side'] = options.probe_len_per_side
    params['is_genome'] = options.is_genome
        
    # output log file
    tools.output_log(__version__, sys.argv, [vars(options)], args[1])
    
    if len(args) == 2:            
        ff = FusionFinder(params, debug=options.debug, unique=options.unique, annodir=options.annodir, mmcfg=options.mmcfg)
        # set up filtering parameters
        ff.filtering['readpairs'] = [options.min_read_pairs, options.max_read_pairs]
        ff.filtering['min_span_reads'] = options.min_span_reads
        ff.filtering['gene_fusions'] = True
        ff.filtering['no_repeats'] = options.no_repeats
        ff.filtering['keep_mito'] = options.keep_mito

        # annotation
        if options.gene_model:
            ff.genome, ff.model = options.gene_model
            ff.prepare_annotation()
            
        # contig sequence
        if options.contigs_file:
            ff.extract_contig_seq(contigs_file=options.contigs_file)
                
        outdir = args[1]
        want_unfiltered = False        
        # parse fusions file
        if options.fusions:
            if (os.path.isdir(args[0])):
                want_unfiltered = True
                ok = ff.parse_fusions_dir(args[0], outdir=outdir)
                if not ok:
                    sys.exit(1)
            else:
                ff.parse_fusions(args[0])
                
        # parse alignment file(s) and find fusions                        
        else:
            if os.path.isdir(args[0]):                
                align_files = sorted(glob.glob(os.path.join(args[0], '*.*')))                
                for i in range(len(align_files)):
                    sys.stdout.write("processing alignment file #%d\n" % (i+1))
                    ff.find(align_files[i])                                                
            else:
                ff.find(args[0])
        
        # group by breakpoint
        ff.group_by_breakpoint()
        
        # realign filtering
        if options.realign and ff.genome:
            ff.realign(options.config_file, options.realigner, outdir=args[1])
            ff.filter_nonspecific(args[1])
            
        # collect read support
        if options.genome_bam or options.contigs_bam:
            ff.read_support(genome_bamfile=options.genome_bam, contigs_bamfile=options.contigs_bam, min_mapq=options.mapq_genome, keep_mito=options.keep_mito, genome_coord_buffer=options.genome_coord_buffer, only_unique=options.unique_spanning_reads)
                        
        ff.report(outdir, filtering=False)
            
        # filtering
        if options.filter:
            # link reciprocals, add spanning reads, overlap dbSNP,DGV, etc
            ff.post_process(options)
            
            # Set events as local
            ff.set_local()
            
            # filter based on read support
            ff.post_filter()

            # output filtered results
            ff.report(outdir, filtering=options.filter, cnv_file=options.cnv_file, unfiltered=want_unfiltered)    

if __name__ == '__main__':
    usage = "Usage: %prog alignment_[file|dir] output_dir"

    parser = OptionParser(usage=usage, version="%prog " + __version__)
    
    io = OptionGroup(parser, "input and output")
    io.add_option("-C", "--contigs_file", dest="contigs_file", help="contigs file")        
    io.add_option("-F", "--fusions", dest="fusions", help="fusions input file given", action="store_true", default=False)
    io.add_option("-g", "--is_genome", dest="is_genome", help="sample is genome", action="store_true", default=False) 
    io.add_option("-V", "--cnv_file", dest="cnv_file")
    parser.add_option_group(io)
    
    operation = OptionGroup(parser, "operations")
    operation.add_option("-X", "--filter", dest="filter", help="post-filter", action="store_true", default=False)
    operation.add_option("-d", "--debug", dest="debug", help="debug", action="store_true", default=False)
    parser.add_option_group(operation)
        
    annotations = OptionGroup(parser, "annotations")
    annotations.add_option("--annodir", dest="annodir", help="the Trans-ABySS 'annotation' directory")
    annotations.add_option("-G", "--gene_model", dest="gene_model", help="gene model used for annotation e.g. hg18 k where k=known_genes, e=ensembl, r=refseq", nargs=2)
    annotations.add_option("-O", "--olap_annot", dest="olap_annot", help="overlaps dbsnp and dgv", action="store_true", default=False)
    annotations.add_option("--mmcfg", dest="mmcfg", help="the path of `model_matcher.cfg'")
    parser.add_option_group(annotations)
    
    alignment_params = OptionGroup(parser, "alignment params")
    alignment_params.add_option("-m", "--max_match_fraction", dest="match_fraction", help="maximum match fraction for individual alignment. Default:0.999", default=0.999, type='float')
    alignment_params.add_option("-q", "--max_query_olap", dest="query_olap", help="maximum query overlap in bases. Default:40", default=40, type='int')
    alignment_params.add_option("-Q", "--max_query_olap_single", dest="query_olap_single", help="maximum single region query overlap. Default:0.3", default=0.3, type='float')
    alignment_params.add_option("-f", "--min_fusion_fraction", dest="fusion_fraction", help="minimum fusion fraction. Default:0.90", default=0.90, type='float')
    alignment_params.add_option("-a", "--num_aligns", dest="num_aligns", help="number alignments. Default:20", default=20, type='int')
    alignment_params.add_option("-A", "--num_best", dest="num_best", help="number of best alignments to use. Default:5", default=5, type='int')
    alignment_params.add_option("-i", "--min_identity", dest="min_identity", help="minimum identity. Default:98", type='float', default=98.0)
    parser.add_option_group(alignment_params)
    
    read_support = OptionGroup(parser, "read support")
    read_support.add_option("-B", "--genome_bam", dest="genome_bam", help="genomic bam file")
    read_support.add_option("-b", "--contigs_bam", dest="contigs_bam", help="contigs bam file")
    read_support.add_option("-s", "--min_span_reads", dest="min_span_reads", help="minimum spanning reads. Default:2", type='int', default=2)
    read_support.add_option("-D", "--max_span_reads", dest="max_span_reads", help="maximum spanning reads. Default:None", type='int', default=None)
    read_support.add_option("-p", "--min_read_pairs", dest="min_read_pairs", help="minimum read pairs. Default:2", type='int', default=2)
    read_support.add_option("-P", "--max_read_pairs", dest="max_read_pairs", help="maximum read pairs. Default:2000", type='int', default=2000)   
    read_support.add_option("-x", "--min_ctg_bp_overlap", dest="min_ctg_bp_overlap", help="minimum contig breakpoint overlap. Default:1", type='int', default=1)
    read_support.add_option("-e", "--min_genome_bp_overlap", dest="min_genome_bp_overlap", help="minimum genome breakpoint overlap. Default:1", type='int', default=1)
    read_support.add_option("-Y", "--mapq_genome", dest="mapq_genome", help="minimum mapping quality for reads to genome. Default:1", type='int', default=0)
    read_support.add_option("-y", "--mapq_contigs", dest="mapq_contigs", help="minimum mapping quality for reads to contigs. Default:0", type='int', default=0)
    read_support.add_option("-c", "--genome_coord_buffer", dest="genome_coord_buffer", help="genome coordinate buffer for extracting reads. Default:2000", default=2000, type='int')
    read_support.add_option("-U", "--unique_spanning_reads", dest="unique_spanning_reads", help="only count unique spanning reads and breakpoint pairs", action="store_true", default=False)
    read_support.add_option("-S", "--add_spanning_reads", dest="add_spanning_reads", action="store_true", default=False)
    parser.add_option_group(read_support)
     
    realign = OptionGroup(parser, "realign")
    realign.add_option("--cfg", dest="config_file", help="the Trans-ABySS configuration file")
    realign.add_option("-R", "--realigner", dest="realigner", help="re-alignment program. Default:blat", default='blat')
    realign.add_option("-r", "--realign", dest="realign", action="store_true", default=False)
    realign.add_option("-l", "--max_subseq_len", dest="max_subseq_len", help="maximum subsequence length", type='int')
    realign.add_option("-o", "--probe_len_per_side", dest="probe_len_per_side", help="probe length on each side of breakpont. Default:50", type='int', default=50)
    parser.add_option_group(realign)
    
    filtering = OptionGroup(parser, "other filtering")
    filtering.add_option("-I", "--filter_repeats", dest="no_repeats", help="no repeat at junction", action="store_true", default=False)
    filtering.add_option("-M", "--keep_mito", dest="keep_mito", help="keep mitochondria", action="store_true", default=False)
    filtering.add_option("-u", "--unique", dest="unique", help="only keep unique pairs (not same target)", action="store_true", default=False)
    parser.add_option_group(filtering)
    
    (options, args) = parser.parse_args()
    main(args, options)
                
