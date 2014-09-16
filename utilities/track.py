"""
This module contains methods for manipulating UCSC tracks (psl, gff, bed).

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
import sys
import re
import alignment
from align_parsers import psl

class Track:
    """Provides methods for generating and parsing UCSC tracks such as PSL, BED, and GFF"""
    def __init__(self, track=None):
        self.track = track

    def parse_gff(self, just_blocks=False):
        """Extracts alignment or alignment blocks from gff file"""
        record = []
        aligns = []
        all_blocks = {}
        
        start = False
        for line in open(self.track, 'r'):
        
            if line[:3].lower() == "chr":
                fields = line.split("\t")

                #hack for merged contigs!
                if len(fields) == 5:
                    fields = [fields[0]+fields[1],fields[2],fields[3],fields[4]]

                if not record or record[-1][-1] == fields[-1]:
                    record.append(fields)
                elif record:
                    record.sort(lambda x,y: int(x[3])-int(y[3]))
                    blocks = []
                    splice_sites = []
                    query_blocks = []
                    for r in record:
                        blocks.append([r[3], r[4]])

                        if ":" in r[1] and not just_blocks:
                            qstart, qend, splice_site = r[1].split(":")
                            query_blocks.append([qstart, qend])
                            splice_sites.append(splice_site)

                    if not just_blocks:
                        align = self.create_align(record[0][1], record[0][-1], record[0][0], record[0][6], blocks)
                        align.query_blocks = query_blocks
                        align.splice_sites = splice_sites
                        aligns.append(align)
                    else:
                        all_blocks[record[0][-1]] = blocks
                    
                    record = [fields]
        if record:
            record.sort(lambda x,y: int(x[3])-int(y[3]))
            blocks = []
            splice_sites = []
            query_blocks = []
            for r in record:
                blocks.append([r[3], r[4]])
                
                if ":" in r[1] and not just_blocks:
                    qstart, qend, splice_site = r[1].split(":")
                    query_blocks.append([qstart, qend])
                    splice_sites.append(splice_site)

            if not just_blocks:
                align = self.create_align(record[0][1], record[0][-1], record[0][0], record[0][6], blocks)
                align.query_blocks = query_blocks
                align.splice_sites = splice_sites
                aligns.append(align)
            else:
                all_blocks[record[0][-1]] = blocks

        if just_blocks:
            return all_blocks
        else:
            return aligns

    def parse_psl(self, minimum=False, refseq_2bit=None, refseq=None, splice_motif_file=None, just_blocks=False):
        """Extracts alignment or alignment blocks from psl file"""
        all_blocks = {}

        refseq = refseq
        if refseq_2bit:
            from tools import get_refseq_from_2bit
            refseq = get_refseq_from_2bit(refseq_2bit)
            
        aligns = psl.parse(self.track, noline=True, minimum=minimum, refseq=refseq, splice_motif_file=splice_motif_file)

        for align in aligns:
            align.query, align.coverage, orient, strand = self.get_contig_name(align.query)

            if just_blocks:
                all_blocks[align.query] = align.blocks

        if not just_blocks:
            return aligns
        else:
            return all_blocks
            
    def create_align(self, method, query=None, target=None, strand=None, blocks=None):
        """Creates alignment from track (Deprecated)"""
        align = alignment.Alignment(method)

        align.query, align.coverage, align.orient, strand = self.get_contig_name(query)
        align.target = target

        if strand == '+':
            align.target_strand = True
        else:
            align.target_strand = False

        align.blocks = blocks

        align.tstart = blocks[0][0]
        align.tend = blocks[-1][-1]

        return align

    def get_contig_name(self, label):
        """Extracts contig name from track field (Deprecated)"""
        contig = label.rstrip()
        coverage = None
        orient = None
        strand = None
        
        if re.match("^\d+", label) and not "k" in label:
            fields = re.split(':', label.rstrip())
            contig = fields[0]
            
            if len(fields) > 1:
                coverage = fields[1]
                
        elif re.match("k\d+:\d+$", label):
            pass
        else:
            fields = re.split(':', label.rstrip())
            #print "fields",fields
            if len(fields) >= 3 and fields[-3][0] == 'k':
                contig = ':'.join(fields[-3:-1])
                coverage = fields[-1]

                #get transcript orientation if it's there in gff file
                if len(fields) > 3:
                    orient = fields[-4][-1]
                    strand = fields[-4][0]
                    
            elif len(fields) >= 2 and fields[-2][0] == 'k':
                contig = ':'.join(fields[-2:])

        return contig, coverage, orient, strand
                
    def output(self, aligns, outfile, genome=None, refseq=None):
        """Outputs track(.psl) file"""
        ext = os.path.splitext(outfile)[1]
        out = open(outfile, 'w')

        out.write("track name=\"%s\" description=\"%s\" visibility=%d itemRgb=\"On\"\n" % (self.name, self.desc, 3))
        for align in aligns:
            if ext == ".gff":
                out.write(align.gff("exon"))
            elif ext == ".psl":
                if align.psl_str:
                    out.write(align.psl_str + '\n')
                else:
                    out.write(align.psl(genome=genome, refseq=refseq) + '\n')
                
        out.close()
                   
    @classmethod
    def ucsc_targets(cls, genome, aligns, annodir):
        """Converts targets in alignment objects to UCSC chromosomes"""
        conversion_file = os.path.join(annodir, genome, "ucsc_chr.txt")

        if os.path.exists(conversion_file):
            conversions = {}
            for line in open(conversion_file, 'r'):
                chr_from, chr_to = line.rstrip('\n').split()
                conversions[chr_from] = chr_to

            if conversions:
                for align in aligns:
                    if conversions.has_key(align.target):
                        if align.psl_str:
                            align.psl_str = align.psl_str.replace(align.target, conversions[align.target])
                        align.target = conversions[align.target]

    @classmethod              
    def create_bed(cls, aligns, outfile=None, out=None, name="", desc=""):
        """Creates bed file of given alignments"""
        if outfile:
            out = open(outfile, 'w')        
        if not out:
            return out

        if name != '' or desc != '':
            out.write("track name=\"%s\" description=\"%s\" visibility=2 itemRgb=\"On\"" % (name, desc))
            
        for align in aligns:
            chrom = str(align.target)
            if align.target[:3] != 'chr':
                chrom = 'chr' + chrom
                
            #zero-based
            chromStart = int(align.tstart) - 1
            chromEnd = int(align.tend)
            
            name = align.query
            
            if align.query_strand:
                strand = align.query_strand
            else:
                strand = '+'
                
            score = 0
            
            if strand == '+':
                itemRgb = '0,0,255'
            else:
                itemRgb = '255,0,0'
            
            out.write('\t'.join([chrom, str(chromStart), str(chromEnd), name, str(score), strand, str(chromStart), str(chromEnd), itemRgb]) + '\n')
        
        if outfile:
            out.close()

     
