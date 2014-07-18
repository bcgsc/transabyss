"""
This module provides methods for parsing a TCGA transcript-genome GAF file.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
from optparse import OptionParser
from analysis import transcript
from ensembl import output

def parse(gaf):
    """Parses GAF file"""
    cs = re.compile('CDSstart=(\d+)')
    ce = re.compile('CDSstop=(\d+)')
    info = {}
    
    txts = []
    for line in open(gaf, 'r'):
        cols = line.rstrip('\n').split('\t')
        if cols[2] == 'transcript' and cols[8] == 'genome':
            chrom, genome_coords, strand = cols[14].split(':')
            transcript_coords = cols[13]
            transcript_id = cols[1]
            gene = cols[15]
            exons_transcript = string_to_spans(transcript_coords)
            exons_genome = string_to_spans(genome_coords)
            
            cds_start = None
            cds_end = None
            if 'CDSstart' in cols[18]:
                m = re.search("CDSstart=(\d+)", cols[18])
                if m:
                    cds_start = int(m.group(1))
                m = re.search("CDSstop=(\d+)", cols[18])
                if m:
                    cds_end = int(m.group(1))
                    
            cds_start_genome = None
            cds_end_genome = None
            if cds_start is not None and cds_end is not None:
                cds_start_genome, cds_end_genome = calc_cds(exons_transcript, exons_genome, cds_start, cds_end, strand)
                            
            if gene is not None and len(gene) > 0:
                txt = transcript.Transcript(transcript_id)
                txt.chrom = chrom
                txt.strand = strand
                txt.txStart, txt.txEnd = exons_genome[0][0] - 1, exons_genome[-1][1]
                if cds_start_genome is not None and cds_end_genome is not None:
                    if cds_start_genome < cds_end_genome:
                        txt.cdsStart, txt.cdsEnd = cds_start_genome - 1, cds_end_genome
                    else:
                        txt.cdsStart, txt.cdsEnd = cds_end_genome - 1, cds_start_genome
                else:
                    txt.cdsStart, txt.cdsEnd = txt.txStart, txt.txStart
                txt.gene = gene
                txt.alias = gene
                txt.exons = exons_genome
                txts.append(txt)
                       
    txts.sort(transcript.Transcript.compare_start)
    return txts
                
def string_to_spans(string):
    """Converts a string('a-b,c-d') to spans"""
    spans = []
    for span in string.split(','):
        start, end = span.split('-')
        spans.append([int(start), int(end)])
        
    return spans

def calc_cds(exons_transcript, exons_genome, cds_start, cds_end, strand):
    """Determine genomic location of CDS start and end for a given transcript"""
    start_genome = None
    end_genome = None
    block = None
    offset = None
    for i in range(len(exons_transcript)):
        if cds_start >= exons_transcript[i][0] and cds_start <= exons_transcript[i][1]:
            offset = cds_start - exons_transcript[i][0]
            block = i
            break
    if block is not None and offset is not None:
        if strand == '+':
            start_genome = exons_genome[block][0] + offset
        else:
            block = len(exons_transcript) - block - 1
            start_genome = exons_genome[block][1] - offset
        
    block = None
    offset = None
    for i in range(len(exons_transcript)):
        if cds_end >= exons_transcript[i][0] and cds_end <= exons_transcript[i][1]:
            offset = cds_end - exons_transcript[i][0]
            block = i
            break
    if block is not None and offset is not None:
        if strand == '+':
            end_genome = exons_genome[block][0] + offset
        else:
            block = len(exons_transcript) - block - 1
            end_genome = exons_genome[block][1] - offset
              
    return start_genome, end_genome

if __name__ == '__main__':
    usage = "Usage: %prog gaf-file outfile"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()
    
    if len(args) == 2:
        output(parse(args[0]), args[1])