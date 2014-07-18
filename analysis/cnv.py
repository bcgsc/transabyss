"""
This module extracts CNV data for given coordinate

Author: Readman Chiu rchiu@bcgsc.ca
"""
from utilities import intspan, tools

class CNV:
    """Extracts CNV info and links to fusion events through genomic coordinates"""
    def __init__(self, cnv_file):
        self.cnv_file = cnv_file
    
    def extract(self):
        """Extracts segment info from cnv file"""
        self.segments = {}
        
        for line in open(self.cnv_file, 'r'):
            (chrom, start, end, cn) = line.rstrip('\n').split()
            
            if not self.segments.has_key(chrom):
                self.segments[chrom] = [[int(start), int(end), cn]]
            else:
                self.segments[chrom].append([int(start), int(end), cn])
    
    def overlap(self, chrom, start, end, buffer=20000):
        """Reports copy numbers in overlapping segments with given coordinate"""
        copy_numbers = []
        span = [max(1, int(start)-buffer), int(end)+buffer]
        
        chr_num = tools.get_chr_number(chrom)
        if self.segments.has_key(chr_num):
            for segment in self.segments[chr_num]:
                if intspan.overlap([segment[0], segment[1]], span):
                    copy_numbers.append(segment[2])
                
        if copy_numbers:
            return ','.join(copy_numbers)
        else:
            return 'na'