"""
This module provides methods for determining whether SNV bubbles are mapped to the same contig.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
from utilities.align_parsers import psl

class BubbleMapper:
    """Provides methods for relating SNV bubble to canonical contigs"""
    def __init__(self, align_file):
        self.align_file = align_file
        self.bubble2contig = {}
        
    @classmethod
    def split_name(cls, name):
        """Splits the name of a bubble into its digit and alphabet components"""
        bubble_name = re.compile(r'(?P<num>\S+)(?P<branch>[A-Z])')
        m = bubble_name.search(name)
        if m:
            return m.group('num'), m.group('branch')
        else:
            return None, None
        
    def map_to_contigs(self):
        """Maps bubbles to contigs given a Blat alignment file"""
        mapping = {}
        
        aligns = psl.parse(self.align_file, noline=True)
        for align in aligns:
            num, branch = BubbleMapper.split_name(align.query)
            if not mapping.has_key(num):
                mapping[num] = {}
            if not mapping[num].has_key(branch):
                mapping[num][branch] = {}
            mapping[num][branch][align.target] = [align.tstart, align.match]
        
        # checks if branches map to the same contig
        for num in sorted(mapping.keys()):        
            branches = sorted(mapping[num].keys())
            for i in range(len(branches)):
                targets = mapping[num][branches[i]]
                for target in targets:
                    matches = {mapping[num][branches[i]][target][1]:True}
                    agreed = 1
                    for j in range(len(branches)):
                        if j != i:
                            if mapping[num][branches[j]].has_key(target):
                                agreed += 1
                                matches[mapping[num][branches[j]][target][1]] = True
                                continue
                    
                    if agreed == len(branches) and len(matches.keys()) > 1:
                        for branch in branches:
                            self.bubble2contig[num + branch] = target
                        break
                                        
    def is_bubble_mapped_to_contig(self, name):
        """Returns contig where bubble is mapped to"""
        if self.bubble2contig.has_key(name):
            return self.bubble2contig[name]
        else:
            return False
        