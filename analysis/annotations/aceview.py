"""
This module provides methods for parsing and indexing a UCSC Acembly file.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser
from analysis import transcript
from utilities import tools

fields = {1:"name", 2:"chrom", 3:"strand", 4:"txStart", 5:"txEnd", 6:"cdsStart", 7:"cdsEnd", 8:"exonCount", 9:"exonStarts", 10:"exonEnds", 11:"alias"}

def parse(infile):
    """Parses a UCSC Acembly file"""
    txts = []
    ff = open(infile, 'r')  
    for line in ff.readlines():
        txt = parse_line(line)
        if txt is not None:
            txts.append(txt)
    ff.close()
    
    return txts

def parse_line(line):
    """Parses an individual record of a UCSC Acembly file"""
    cols = line.rstrip("\n").split("\t")
    if cols[0]:
        txt = transcript.Transcript(cols[1])
        for i in range(len(cols)):
            if i in fields:
                if i < 9 or i == 11:
                    setattr(txt, fields[i], cols[i])

        exonStarts = cols[9].rstrip(',').split(',')
        exonEnds = cols[10].rstrip(',').split(',')
        txt.exons = []
        for e in range(len(exonStarts)):
            txt.exons.append([int(exonStarts[e]) + 1, int(exonEnds[e])])

        # calculates transcript length for coverage
        for exon in txt.exons:
            txt.length += int(exon[1]) - int(exon[0]) + 1
                
        return txt

    return None

def index(infile, output, genome):
    """Generates an index file"""
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split("\t")
        start = int(int(cols[4]) / 1000)
        end = int(int(cols[5]) / 1000)
        target = cols[2]
        target = tools.proper_chrom(target, genome)
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
        
def txt2gene(infile):
    """Generates an index file"""
    genes = {}
    count = 0
    for line in open(infile, 'r'):
        count += 1
        cols = line.rstrip('\n').split('\t')
        if not genes.has_key(cols[11]):
            genes[cols[11]] = []
        genes[cols[11]].append([cols[1], count])
        
    return genes

if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")
    parser.add_option("-g", "--genome", dest="genome", help="genome")
    
    (options, args) = parser.parse_args()

    if options.index:
        index(args[0], options.index, options.genome)
