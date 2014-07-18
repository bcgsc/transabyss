"""
This module uses the line cache module to extract annotation features
by coordinates

Author: Readman Chiu rchiu@bcgsc.ca
"""
import linecache
import re

class OverlapCoord:
    """Provides methods for overlapping event coordinates with annotations that have been indexed 
    by genomic coordinates and line numbers
    """
    def __init__(self, source_file, index_file):
        self.source_file = source_file
        self.index_file = index_file
            
    def extract_index(self):
        """Extracts index and store in dictionary"""
        self.index = {}
        for line in open(self.index_file, 'r'):
            coord, line_nums = line.rstrip().split()
            self.index[coord] = line_nums
    
    def overlap(self, chrom, start, end, parse_line=None):
        """Extract all annotation objects that overlap given coordinate
        Uses line numbers in indices to extract corresponding line
        by linecache and then parse the line using the given
        parse_line function to convert it to object
        If parse_line is not provided, will return data lines
        """
        start = int(start)
        end = int(end)
        start_index = ':'.join((chrom, str(int(start/1000))))
        end_index = ':'.join((chrom, str(int(end/1000))))
                 
        line_nums = {}
        if self.index.has_key(start_index):
            for line in self.index[start_index].split(','):
                line_nums[line] = True
        if end_index != start_index:
            int(int(start)/1000)
            for coord in range(int(start/1000)+1,int(end/1000)+1):
                idx = ':'.join((chrom, str(coord)))
                if self.index.has_key(idx):
                    for line in self.index[idx].split(','):
                        line_nums[line] = True

        data = []
        for line_num in sorted(line_nums.keys()):
            line = linecache.getline(self.source_file, int(line_num))
            if not re.match('\w', line):
                continue
            if parse_line:
                d = parse_line(line)
                if d:
                    data.append(d)
            else:
                data.append(line)
            
        return data
    
    def finish(self):
        """Clears linecache"""
        linecache.clearcache()
        