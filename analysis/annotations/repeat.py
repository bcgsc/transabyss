"""
This module provides methods for indexing and overlapping UCSC repeats files:
simple_repeats, all_rmsk, and segdups.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser
from utilities.overlap_coord import OverlapCoord
from utilities.intspan import overlap, subsume

def parse(file):
    """Wrapper function for parsing repeats file"""
    repeats = []
    for line in open(file, 'r'):
        repeat = parse_line(line, simple_repeat=simple_repeat)
        if repeat:
            repeats.append(repeat)
    
    return repeats

def parse_line(line):
    """Parses individual line of UCSC simple repeat file"""
    cols = line.rstrip("\n").split()
    target = cols[0]
    repeat = None
    if not 'chr' in target:
        target = 'chr' + target        
    repeat = {'target':target, 'start':int(cols[1]), 'end':int(cols[2]), 'type':cols[3]}
        
    return repeat

def prepare_overlap(genome, annodir, repeat_types=None):
    """Extracts index info into dictionary"""
       
    genome_dir = os.path.join(annodir, genome)
    
    repeats = {'simple_repeats': genome + "_simple_repeats.coords",
               'rmsk': genome + "_all_rmsk.coords",
               'segdup': genome + "_segdups.coords"
               }
    
    overlaps = {}
    if os.path.isdir(genome_dir):
        for repeat_type, filename in repeats.iteritems():
            if repeat_types and not repeat_type in repeat_types:
                continue
            
            repeat_file = genome_dir + '/' + filename
            index_file = repeat_file + '.idx'
            if os.path.exists(repeat_file):
                overlaps[repeat_type] = OverlapCoord(repeat_file, index_file)
                overlaps[repeat_type].extract_index()
            else:
                print 'cannot find file:' + repeat_file
        
        return overlaps
    
    return None

def find_overlaps(test, repeat_overlaps):
    """Overlaps given coordinates with repeats to identify subsuming(simple repeats, segdups)
    or overlaps(rmsk)
    """
    overlaps = {}
    for repeat_type, repeat_overlap in repeat_overlaps.iteritems():
        overlaps[repeat_type] = {}
        repeats = repeat_overlap.overlap(test['chrom'], test['start'], test['end'], parse_line=parse_line)
    
        if repeat_type == 'simple_repeats' or repeat_type == 'segdup':
            for repeat in repeats:
                if subsume([test['start'], test['end']], [repeat['start'], repeat['end']]):
                    overlaps[repeat_type][repeat['type']] = True
        elif repeat_type == 'rmsk':
            for repeat in repeats:
                if overlap([test['start'], test['end']], [repeat['start'], repeat['end']]):
                    overlaps[repeat_type][repeat['type']] = True
                    
    return overlaps
    
def index(infile, output):
    """Indexes repeat file"""
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split()

        start = int(int(cols[1])/1000)
        end = int(int(cols[2])/1000)
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
