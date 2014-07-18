"""
This module provides methods for parsing and indexing a UCSC ensGene file.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser
from analysis import transcript

fields = {0:"name", 1:"chrom", 2:"strand", 3:"txStart", 4:"txEnd", 5:"exonCount", 6:"exonStarts", 7:"exonEnds"}

def parse(infile):
    """Parses old ensGene file"""
    txts = []
    ff = open(infile, 'r')  
    for line in ff.readlines():
        txt = parse_line(line)
        if txt is not None:
            txts.append(txt)
    ff.close()
    
    return txts

def parse_line(line):
    """Parse individual lines of an ensGene file"""
    cols = line.split("\t")
    if cols[0]:
        txt = transcript.Transcript(cols[1])
        for i in range(len(cols)):
            if i in fields:
                if i < 8:
                    setattr(txt, fields[i], cols[i])

        exonStarts = cols[6].rstrip(',').split(',')
        exonEnds = cols[7].rstrip(',').split(',')
        txt.exons = []
        for e in range(len(exonStarts)):
            txt.exons.append([int(exonStarts[e])+1, int(exonEnds[e])])

        #calculate transcript length for coverage
        for exon in txt.exons:
            txt.length += int(exon[1]) - int(exon[0]) + 1

        return txt

    return None

def index(infile, output):
    """Generates an index file"""
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split("\t")
        start = int(int(cols[3]) / 1000)
        end = int(int(cols[4]) / 1000)
        target = cols[1]
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

if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()
    if options.index:
        index(args[0], options.index)
