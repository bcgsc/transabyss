"""
This module provides methods for parsing and indexing a UCSC repeatmasker file.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser

fields = {5:"chrom", 6:"start", 7:"end", 10:"family"}

def parse(file):
    """Parses a UCSC repeatmasker file"""
    repeats = []
    for line in open(file, 'r'):
        repeat = parse_line(line)
        if repeat:
            repeats.append(repeat)
    
    return repeats

def parse_line(line):
    """Parses individual lines of a UCSC repeatmasker file"""
    cols = line.rstrip("\n").split(" ")

    if len(cols) >= 5:
        target = cols[0]
        if not 'chr' in target:
            target = 'chr' + target
        repeat = {'target':target, 'start':int(cols[1]), 'end':int(cols[2]), 'type':cols[4]}
        return repeat

    return None

def index(infile, output):
    """Generates an index file"""
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split(" ")
        start = int(int(cols[1]) / 1000)
        end = int(int(cols[2]) / 1000)
        target = cols[0]
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


if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()
    if options.index:
        index(args[0], options.index)
