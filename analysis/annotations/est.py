"""
This module provides methods for parsing and indexing a UCSC est file.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
import re
from optparse import OptionParser
from analysis import transcript
from utilities.align_parsers import psl
from utilities.overlap_coord import OverlapCoord

def parse_line(line):
    """Parses an individual record of a UCSC est file"""
    align = psl.create_align(re.sub('^\d+\t', '', line))
    
    if align:
        txt = transcript.Transcript(align.query)
        txt.chrom = align.target
        txt.strand = align.query_strand
        txt.txStart = align.tstart
        txt.txEnd = align.tend
        txt.exons = align.blocks
        txt.exonCount = len(align.blocks)
        return txt
    else:
        return None

def index(infile, output):
    """Generates an index file"""
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split("\t")
        start = int(int(cols[16]) / 1000)
        end = int(int(cols[17]) / 1000)
        target = cols[14]
        if not 'chr' in target:
            target = 'chr' + target

        for n in range(start,end+1):
            index = ':'.join((target,str(n)))
            value = str(line_num)

            if not indices.has_key(index):
                indices[index] = [value]
            else:
                indices[index].append(value)

        line_num += 1

    index_file = open(output, 'w')
    for index in sorted(indices.keys()):
        index_file.write(' '.join((index, ','.join(indices[index]))) + "\n")
        
def prepare_overlap(genome, annodir):
    """Converts indexing info into a dictionary"""
    
    genome_dir = os.path.join(annodir, genome)
    
    if os.path.isdir(genome_dir):
        est_file = genome_dir + '/' + 'all_est.txt'
        index_file = est_file + '.idx'
        est_overlap = OverlapCoord(est_file, index_file)
        est_overlap.extract_index()
        
        return est_overlap
    
    return None

if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()
    if options.index:
        index(args[0], options.index)
