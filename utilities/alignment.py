"""
This module provides generalized methods for dealing with alignment

Author: Readman Chiu rchiu@bcgsc.ca
"""
import sys
import re
from assembly import Contig
from utilities import tools

class Alignment:
    """Stores alignment info, provides methods for extracting alignment info,
    and adds alignment info that may not be directly accessible from aligner's output (e.g.
    splice sites)
    """
    def __init__(self, method):
        self.method = method
        self.blocks = self.query_blocks = self.splice_sites = []
        self.mismatch = None
        
        self.line = None
        self.method = None
        self.query = target = None
        self.query_len = target_len = None
        self.target_strand = None
        self.query_strand = None
        self.qstart = qend = tstart = tend = None
        self.blocks = None
        self.splice_sites = None
        self.query_blocks = None
        self.match_len = matches = mismatches = None
        self.identity = None
        self.score = None
        self.rank = None
        self.model = None
        self.pairwise = None
        self.psl_str = None
        self.orient = None
        self.contig = None
        self.num_dels = self.bases_dels = self.num_ins = self.bases_ins = None
        self.target_size = None
        
    def get_splice_sites(self, refseq):
        """Extracts splice sites of every gap in alignment"""
        self.splice_sites = []
        for i in range(len(self.blocks)-1):
            if self.blocks[i+1][0] - self.blocks[i][1] < 5:
                ss = "NA"
            else:
                intron = refseq.GetSequence(self.target, self.blocks[i][1] + 1, self.blocks[i+1][0]-1)
                ss = intron[:2].lower() + intron[-2:].lower()
                
            self.splice_sites.append(ss)
        
    def set_orient(self, splice_motifs):
        """Sets orientation(+/-) of contig based on splice motifs"""
        if not splice_motifs:
            return None
    
        counts = {"forward":0, "backward":0, "unknown":0}

        motifs = {'forward':[], 'backward':[]}
        for motif in splice_motifs.keys():
            motifs['forward'].append(motif.lower())
            motifs['backward'].append(tools.reverse_complement(motif).lower())

        for ss in self.splice_sites:
            orient = "unknown"

            if ss in motifs['forward']:
                orient = "forward"
            elif ss in motifs['backward']:
                orient = "backward"

            counts[orient] += 1

        orient = None
        if counts['forward'] > 0 and counts['backward'] == 0 and counts['forward'] > counts['unknown']:
            orient = '+'
        elif counts['backward'] > 0 and counts['forward'] == 0 and counts['backward'] > counts['unknown']:
            orient = '-'

        self.orient = orient
        
    def qpos_to_tpos(self, qpos):
        """Maps query postion to target position"""
        
        # find corresponding block
        block = None
        for i in range(len(self.query_blocks)):
            if (qpos >= self.query_blocks[i][0] and qpos <= self.query_blocks[i][1]) or\
               (qpos <= self.query_blocks[i][0] and qpos >= self.query_blocks[i][1]):
                block = i
                break
        
        tpos = None
        if block is not None:
            if self.query_strand == '+':
                tpos = self.blocks[block][0] + qpos - self.query_blocks[block][0]
            else:
                tpos = self.blocks[block][1] - (qpos - self.query_blocks[block][1])
                
        return tpos
        
    def gff(self, type):
        """Outputs alignment in gff format"""
        if not self.blocks:
            print "no blocks", self.query, self.target
            return ""
        
        #chromosome
        if 'chr' in self.target:
            chrom = self.target
        else:
            chrom = 'chr' + self.target

        #contig
        contig = self.query.split(" ")[0]
        if self.contig and self.contig.num:
            contig = str(self.contig.num)

        #strand
        if self.target_strand:
            strand = "+"
        else:
            strand = "-"

        #gene orientation
        orient = "?"
        if self.orient:
            if self.orient == 'forward':
                orient = "+"
            elif self.orient == 'backward':
                orient = "-"

        group = ":".join([self.target, str(self.tstart), "%s%s" % (strand, orient), contig])   
        
        #do this after recording strand info in group
        if orient != "?":
            strand = orient

        out = ""     
        for i in range(len(self.blocks)):
            qstart, qend = self.query_blocks[i]
            tstart, tend = self.blocks[i]
            
            splice_site = ""
            if len(self.blocks) > 1 and i < len(self.blocks)-1 and len(self.blocks):
                splice_site = self.splice_sites[i]

            #hide query block and splice sites in "program" field
            annot = ":".join([str(qstart), str(qend), splice_site])

            out += "\t".join([chrom, annot, type, tstart, tend, ".", strand, ".", group]) + "\n"

        return out
            
    def psl(self, refseq=None, genome=None):
        """Outputs alignment in UCSC psl format"""
        if self.psl_str:
            psl_str = self.psl_str

            cols = self.psl_str.split("\t")
            # make sure target has chr in it for UCSC
            if not re.match('^(chr|scaffold)', cols[13], re.IGNORECASE):
                cols[13] = 'chr' + cols[13]
            
            # gene orientation
            if self.orient:
                if self.orient == 'forward':
                    orient = "+"
                elif self.orient == 'backward':
                    orient = "-"
                cols[8] = orient
                
            psl_str = "\t".join(cols)
            
        else:
            # first few fields are not accurate
            fields = []                        
            fields.append(str(self.matches))
            fields.append(str(self.mismatches))
            
            # repMatches = not available, use 0
            fields.append('0')
            
            ncount = 0
            if refseq:
                tseq = refseq.GetSequence(self.target, int(self.tstart), int(self.tend))
                ncount = len(re.findall('n', tseq, re.IGNORECASE))
            fields.append(str(ncount))
            
            fields.append(str(self.num_ins))
            fields.append(str(self.bases_ins))
            fields.append(str(self.num_dels))
            fields.append(str(self.bases_dels))
                            
            target = self.target
            if genome is not None:
                target = tools.proper_chrom(target, genome)
                
            target_len = self.target_size
            fields.extend([self.query_strand, self.query, str(self.query_len), str(int(self.qstart)-1), str(self.qend), target, str(target_len), str(int(self.tstart)-1), str(self.tend)])
            
            fields.append(str(len(self.blocks)))
            
            block_sizes = [str(t[1] - t[0] + 1) for t in self.blocks]
            fields.append(','.join(block_sizes) + ',')
            
            if self.query_strand == '+':
                qstarts = [str(q[0] - 1) for q in self.query_blocks]
            else:
                qstarts = [str(self.query_len - q[0]) for q in self.query_blocks]
            fields.append(','.join(qstarts) + ',')
                
            tstarts = [str(t[0] - 1) for t in self.blocks]
            fields.append(','.join(tstarts) + ',')
            
            psl_str = "\t".join(fields)
                        
        return psl_str + "\n"
                           
    def correct_blocks(self, splice_motifs, refseq, query_seq): 
        """Correct blocks and gaps of Blat alignment"""
        if not self.splice_sites or not splice_motifs:
            return False
                
        #fix neighboring gaps
        self.correct_neighbor_gaps(splice_motifs, refseq)
                        
        #fix single gaps
        self.correct_single_gaps(splice_motifs, refseq)
                        
    def correct_single_gaps(self, splice_motifs, refseq):
        """
        Post-process blocks after fix_single_gaps()
        """
        gaps = {}
        
        for i in range(len(self.blocks)-1):
            ss = self.splice_sites[i]
            if ss and not splice_motifs.has_key(ss) and not splice_motifs.has_key(tools.reverse_complement(ss).lower()):                    
                gaps[i] = []
                if self.query_blocks[i][1] < self.query_blocks[i+1][0]:
                    for j in range(self.query_blocks[i][1]+1, self.query_blocks[i+1][0]):
                        gaps[i].append(j)
                else:
                    for j in range(self.query_blocks[i][1]-1, self.query_blocks[i+1][0], -1):
                        gaps[i].append(j)
                
        if gaps:
            target_blocks = self.blocks[:]
            query_blocks = self.query_blocks[:]
            splice_sites = self.splice_sites[:]
            
            gap_indices = gaps.keys()
            gap_indices.sort(lambda x,y: x-y)
            
            changed_blocks = {}
            
            for i in gap_indices:
                tblock1 = target_blocks[i][:]
                tblock2 = target_blocks[i+1][:]
                qblock1 = query_blocks[i][:]
                qblock2 = query_blocks[i+1][:]
                
                splice_site, new_block = self.fix_single_gap(tblock1, tblock2, qblock1, qblock2, splice_motifs, refseq, self.query_strand, gaps[i])

                if splice_site:          
                    sys.stderr.write("Type2a %s changed blocks %s %s to %s %s\n" % (self.query, self.target, target_blocks[i], tblock1, splice_site))
                    sys.stderr.write("Type2a %s changed blocks %s %s to %s %s\n" % (self.query, self.target, target_blocks[i+1], tblock2, splice_site))
                    target_blocks[i] = tblock1
                    target_blocks[i+1] = tblock2
                    query_blocks[i] = qblock1
                    query_blocks[i+1] = qblock2
                    splice_sites[i] = splice_site
                    
                    if new_block:
                        changed_blocks[i] = new_block
                        
            if changed_blocks:
                changed_blocks_indices = changed_blocks.keys()
                changed_blocks_indices.sort(lambda x,y: y-x)
                
                for i in changed_blocks_indices:
                    if changed_blocks[i] != -1:
                        sys.stderr.write("Type2b %s add block %s %s %s %s\n" % (self.query, self.target, changed_blocks[i]['target'], changed_blocks[i]['query'], 'NA'))
                        if abs(target_blocks[i][1] - changed_blocks[i]['target'][0]) < abs(target_blocks[i+1][0] - changed_blocks[i]['target'][1]):
                            splice_sites.insert(i, 'NA')
                        else:
                            splice_sites.insert(i+1, 'NA')
                        
                        query_blocks.insert(i+1, changed_blocks[i]['query'])
                        target_blocks.insert(i+1, changed_blocks[i]['target'])
                    else:
                        del query_blocks[i]
                
    
            if target_blocks != self.blocks and self.check_corrections(query_blocks, target_blocks, self.query_strand, '+', self.query):
                self.blocks = target_blocks[:]
                self.query_blocks = query_blocks[:]
                self.splice_sites = splice_sites
                if not self.mismatch or int(self.mismatch) == 0:
                    self.mismatch = 1
                  
                    
    def correct_neighbor_gaps(self, splice_motifs, refseq):
        """Post-process blocks after fix_neighbor_gaps()"""
        gaps = {}
        
        for i in range(len(self.blocks)-1):
            ss = self.splice_sites[i]

            if ss and not splice_motifs.has_key(ss) and not splice_motifs.has_key(tools.reverse_complement(ss).lower()):
                if abs(self.query_blocks[i+1][0] - self.query_blocks[i][1]) == 1:
                    gaps[i] = 0

        if gaps:
            target_blocks = self.blocks[:]
            query_blocks = self.query_blocks[:]
            splice_sites = self.splice_sites[:]
            
            gap_indices = gaps.keys()
            gap_indices.sort(lambda x,y: x-y)
            
            # fix by moving exon and then shuffle
            replaced = {}
            replaced_ordered = []
            for i in range(len(gap_indices)):
                i1 = gap_indices[i]
                i2 = i1 + 1
                i0 = i1 - 1
                
                if gaps[i1] == 1:
                    continue
                
                if i2 + 1 < len(target_blocks):
                    tblock1 = target_blocks[i1][:]
                    tblock2 = target_blocks[i2][:]
                    tblock3 = target_blocks[i2+1][:]
                    qblock1 = query_blocks[i1][:]
                    qblock2 = query_blocks[i2][:]
                    qblock3 = query_blocks[i2+1][:]
                    splice_site, new_block = self.fix_neighbor_gaps(tblock1, tblock2, tblock3, qblock1, qblock2, qblock3, splice_motifs, refseq, self.query_strand)
                    if splice_site:
                        idx = ' '.join((str(i1), str(i2), str(i2+1)))
                        if new_block:
                            replaced[idx] = tblock1, tblock3, qblock1, qblock3, splice_site, new_block
                        else:
                            replaced[idx] = tblock1, tblock3, qblock1, qblock3, splice_site
                            
                        gaps[i1] = 1
                        if gaps.has_key(i2):
                            gaps[i2] = 1
                        replaced_ordered.append(idx)
                        
                # if not fixed, try backward
                if gaps[i1] == 0 and i0 >= 0 and i2 <= len(target_blocks) and (not gaps.has_key(i0) or gaps[i0] == 0):
                    tblock1 = target_blocks[i0][:]
                    tblock2 = target_blocks[i1][:]
                    tblock3 = target_blocks[i2][:]
                    qblock1 = query_blocks[i0][:]
                    qblock2 = query_blocks[i1][:]
                    qblock3 = query_blocks[i2][:]                    
                    splice_site, new_block = self.fix_neighbor_gaps(tblock1, tblock2, tblock3, qblock1, qblock2, qblock3, splice_motifs, refseq, self.query_strand)
                    if splice_site:
                        idx = ' '.join((str(i0), str(i1), str(i2)))
                        if new_block:
                            replaced[idx] = tblock1, tblock3, qblock1, qblock3, splice_site, new_block
                        else:
                            replaced[idx] = tblock1, tblock3, qblock1, qblock3, splice_site
                            
                        gaps[i1] = 1
                        if gaps.has_key(i0):
                            gaps[i0] = 1
                        replaced_ordered.append(idx)

            # make sure delete from back to front
            replaced_ordered.reverse()
            for indices in replaced_ordered:
                new_blocks = replaced[indices]

                ok = True
                for index in indices.split(' '):
                    if gaps.has_key(int(index)) and gaps[int(index)] > 1:
                        ok = False
                        break

                if ok:
                    idx = [int(i) for i in indices.split(' ')]
                    sys.stderr.write("Type3 %s changed blocks %s %s to %s %s\n" % (self.query, self.target, self.blocks[idx[0]], new_blocks[0], new_blocks[4]))
                    sys.stderr.write("Type3 %s changed blocks %s %s to %s %s\n" % (self.query, self.target, self.blocks[idx[2]], new_blocks[1], new_blocks[4]))
                    if len(new_blocks) == 5:
                        sys.stderr.write("Type3 %s removed block %s %s\n" % (self.query, self.target, self.blocks[idx[1]]))
                    else:
                        sys.stderr.write("Type3 %s changed blocks %s %s to %s\n" % (self.query, self.target, self.blocks[idx[1]], new_blocks[-1]['target']))
                    target_blocks[idx[0]] = new_blocks[0]
                    target_blocks[idx[2]] = new_blocks[1]
                    query_blocks[idx[0]] = new_blocks[2]
                    query_blocks[idx[2]] = new_blocks[3]
                    
                    if len(new_blocks) == 5:
                        del target_blocks[idx[1]]
                        del query_blocks[idx[1]]
                        del splice_sites[idx[1]]
                        splice_sites[idx[0]] = new_blocks[4]
                    else:
                        target_blocks[idx[1]] = new_blocks[-1]['target']
                        query_blocks[idx[1]] = new_blocks[-1]['query']
                        
                        if target_blocks[idx[1]][0] - target_blocks[idx[0]][1] < target_blocks[idx[2]][0] - target_blocks[idx[1]][0]:
                            splice_sites[idx[0]] = 'NA'
                            splice_sites[idx[1]] = new_blocks[4]
                        else:
                            splice_sites[idx[0]] = new_blocks[4]
                            splice_sites[idx[1]] = 'NA'

            if target_blocks != self.blocks and self.check_corrections(query_blocks, target_blocks, self.query_strand, '+', self.query):
                self.blocks = target_blocks[:]
                self.query_blocks = query_blocks[:]
                self.splice_sites = splice_sites

                if not self.mismatch or int(self.mismatch) == 0:
                    self.mismatch = 1
                    
    def fix_single_gap(self, tblock1, tblock2, qblock1, qblock2, splice_motifs, refseq, query_strand, extra_query=[]):
        """Shuffles sequence from end to end to see if canonical splice sites can be achieved"""
        max_shuffle_size = 10
        
        tgap = [tblock1[1]+1, tblock2[0]-1]
        tsize = tgap[1] - tgap[0] + 1

        min_size = 10
        max_size = 100000
        if tsize < min_size or tsize > max_size:
            return False, None
        
        fixed = False
        new_block = {'query':[], 'target':[]}
        
        possible_shuffles = []
        for left_shuffle in range(-1 * max_shuffle_size, max_shuffle_size + 1):                        
            for right_shuffle in range(-1 * max_shuffle_size, max_shuffle_size + 1):
                #skip shuffle in opposite directions, and cases where one side is 0 but there is no extra sequence to move
                if left_shuffle * right_shuffle < 0 or (left_shuffle * right_shuffle == 0 and len(extra_query) == 0):
                    continue
                                
                coord = tgap[0] + left_shuffle, tgap[1] + right_shuffle
                gap_seq = refseq.GetSequence(self.target, coord[0], coord[1])
                ss = gap_seq[:2] + gap_seq[-2:]
                
                if splice_motifs.has_key(ss.lower()) or splice_motifs.has_key(tools.reverse_complement(ss).lower()):
                    possible_shuffles.append({'motif':ss.lower(), 'left':left_shuffle, 'right':right_shuffle, 'shuffle_size':abs(left_shuffle) + abs(right_shuffle)})
        
        splice_site = None
        left_shuffle = right_shuffle = None
        
        if possible_shuffles:
            possible_shuffles.sort(self.compare_shuffles)        
            splice_site = possible_shuffles[0]['motif']
            left_shuffle = possible_shuffles[0]['left']
            right_shuffle = possible_shuffles[0]['right']
            
            successful_shuffles = []
            for i in range(len(possible_shuffles)):
                success, new_block = self.shuffle(possible_shuffles[i]['left'], possible_shuffles[i]['right'], extra_query[:], qblock1, qblock2, tblock1, tblock2, query_strand)
                if success:
                    splice_site = possible_shuffles[i]['motif']
                    print 'success', splice_site, new_block, qblock1, qblock2, tblock1, tblock2
                    break
                
            if splice_site:
                return splice_site, new_block
            else:
                return False, None
            
        else:
            return False, None
        
    def shuffle(self, left_shuffle, right_shuffle, extra_query, qblock1, qblock2, tblock1, tblock2, query_strand):
        """Shuffle gap coordinates to see if canonical splice-sites can be rescued"""
        fixed = False
        new_block = {'query':[], 'target':[]}
        
        if left_shuffle != None and right_shuffle != None:
            fixed = True
            
            # move bases from right to left
            if left_shuffle >= 0 and right_shuffle >= 0:
                if query_strand == '+':
                    for i in range(0, right_shuffle):
                        extra_query.append(qblock2[0] + i)
                    new_block['query'].append(max(extra_query[0], extra_query[-1] - left_shuffle + 1))
                    new_block['query'].append(extra_query[-1])
                    qblock2[0] += right_shuffle

                else:
                    for i in range(0, right_shuffle):
                        extra_query.append(qblock2[0] - i)
                    new_block['query'].append(min(extra_query[0], extra_query[-1] + left_shuffle - 1))
                    new_block['query'].append(extra_query[-1])
                    qblock2[0] -= right_shuffle
                                                
                tblock2[0] += right_shuffle
                new_block['target'].append(tblock1[1] + max(1, left_shuffle - max(right_shuffle, len(extra_query)) + 1))
                new_block['target'].append(tblock1[1] + max(1, left_shuffle))
                        
            # move bases from left to right
            elif left_shuffle <= 0 and right_shuffle <= 0:
                if query_strand == '+':
                    for i in range(abs(left_shuffle)):
                        extra_query.insert(0, qblock1[1] - i)
                    new_block['query'].append(extra_query[0])
                    new_block['query'].append(min(extra_query[-1], extra_query[0] + abs(right_shuffle) - 1))
                    qblock1[1] += left_shuffle
                    
                else:                        
                    for i in range(abs(left_shuffle)):
                        extra_query.insert(0, qblock1[1] + i)
                                                             
                    new_block['query'].append(extra_query[0])
                    new_block['query'].append(max(extra_query[-1], extra_query[0] - abs(right_shuffle) + 1))
                    qblock1[1] -= left_shuffle
                                                
                tblock1[1] += left_shuffle
                new_block['target'].append(tblock2[0] - max(1, abs(right_shuffle)))
                new_block['target'].append(tblock2[0] - max(1, abs(right_shuffle) - max(abs(left_shuffle), len(extra_query)) + 1))
            else:
                fixed = False                    
                
        #check if query coordinate make sense
        if fixed:
            # make sure nothing is less than zero
            for block in (qblock1, qblock2, new_block['query']):
                if block:
                    if block[0] < 0 or block[1] < 0:
                        fixed = False
                        break
            
            # make sure everything is ascending or descending
            if new_block['query']:
                if query_strand == '+' and not (qblock1[0] <= qblock1[1] and qblock1[1] <= new_block['query'][0] and 
                        new_block['query'][0] <= new_block['query'][1] and new_block['query'][1] <= qblock2[0] and 
                        qblock2[0] <=qblock2[1]):
                    fixed = False
                    
                if query_strand == '-' and not (qblock1[0] >= qblock1[1] and qblock1[1] >= new_block['query'][0] and 
                        new_block['query'][0] >= new_block['query'][1] and new_block['query'][1] >= qblock2[0] and 
                        qblock2[0] >= qblock2[1]):
                    fixed = False
                    
        if fixed:
            return True, new_block
        else:
            return False, None
        
    def compare_shuffles(self, s1, s2):
        """Ranks shuffles based on shuffling size first, and then splice-motif"""
        if s1['shuffle_size'] < s2['shuffle_size']:
            return -1
        elif s1['shuffle_size'] > s2['shuffle_size']:
            return 1
        else:
            if (s1['motif'].lower() == 'gtag' or s1['motif'].lower() == 'ctac') and s2['motif'].lower() != 'gtag' and s2['motif'].lower() != 'ctac':
                return -1
            elif (s2['motif'].lower() == 'gtag' or s2['motif'].lower() == 'ctac') and s1['motif'].lower() != 'gtag' and s1['motif'].lower() != 'ctac':
                return 1
            else:
                return 0
                    
    def fix_neighbor_gaps(self, tblock1, tblock2, tblock3, qblock1, qblock2, qblock3, splice_motifs, refseq, query_strand):
        """Fix neigboring gaps if necessary to see if canonical splice sites can be re-established.
        In case of consecutive gaps that don't have canonical splice sites,
        move middle block to either end to see if canonical splice can be achieved
        (without or without shuffling after movement)
        """        
        middle_block_size = tblock2[1] - tblock2[0] + 1
        max_shuffle_size = 10
        if middle_block_size > max_shuffle_size:
            return False, None
        
        tgap = [tblock1[1]+1, tblock3[0]-1]
        new_block = {'query':[], 'target':[]}
        
        possible_shuffles = []
        for shuffle in (-1 * middle_block_size, middle_block_size):
            if shuffle > 0:
                left_shuffle, right_shuffle = shuffle, 0
            else:
                left_shuffle, right_shuffle = 0, shuffle
            coord = tgap[0] + left_shuffle, tgap[1] + right_shuffle
            gap_seq = refseq.GetSequence(self.target, coord[0], coord[1])
            ss = gap_seq[:2] + gap_seq[-2:]
            
            if splice_motifs.has_key(ss.lower()) or splice_motifs.has_key(tools.reverse_complement(ss).lower()):
                possible_shuffles.append({'motif':ss.lower(), 'left':left_shuffle, 'right':right_shuffle, 'shuffle_size':abs(left_shuffle) + abs(right_shuffle)})
                    
        splice_site = None
        left_shuffle = right_shuffle = None
                
        if possible_shuffles:
            possible_shuffles.sort(self.compare_shuffles)        
            splice_site = possible_shuffles[0]['motif']
            left_shuffle = possible_shuffles[0]['left']
            right_shuffle = possible_shuffles[0]['right']        
                
        if splice_site and left_shuffle != None and right_shuffle != None:
            #move bases from right to left
            if right_shuffle == 0:
                tblock1[1] += left_shuffle
                if query_strand == '+':
                    qblock1[1] += left_shuffle
                else:
                    qblock1[1] -= left_shuffle
            else:
                tblock3[0] += right_shuffle
                if query_strand == '+':
                    qblock3[0] += right_shuffle
                else:
                    qblock3[0] -= right_shuffle
        else:
            splice_site, new_block = self.fix_single_gap(tblock1, tblock3, qblock1, qblock3, splice_motifs, refseq, query_strand, extra_query=range(qblock2[0], qblock2[1]+1))
            
        if new_block != None and new_block['query']:
            return splice_site, new_block
        else:
            return splice_site, None
            
    def check_corrections(self, query_blocks, target_blocks, query_strand, target_strand, query):
        """Checks if after corrections block coordinates make sense"""
        passed = True
        
        strands = (query_strand, target_strand)
        blocks_list = (query_blocks, target_blocks)
        for b in range(len(blocks_list)):
            blocks = blocks_list[b]
            if strands[b] == '+':
                for i in range(len(blocks)-1):
                    if not (blocks[i][0] <= blocks[i][1]) or not (blocks[i][1] <= blocks[i+1][0]) or not (blocks[i+1][0] <= blocks[i+1][1]):
                        passed = False
                        break
                        
            else:
                for i in range(len(blocks)-1):
                    if not (blocks[i][0] >= blocks[i][1]) or not (blocks[i][1] >= blocks[i+1][0]) or not (blocks[i+1][0] >= blocks[i+1][1]):
                        passed = False
                        break
                    
        if not passed:
            sys.stderr.write("Failed correcting blocks:" + query + '\n')
        else:
            sys.stderr.write("Succeeded correcting blocks:" + query + '\n')
            
        return passed
                    
    def merge_blocks(self, query=True):
        """Merge target blocks that are consecutive
        This may not be necessary but seems sensible to do.
        Do this for model_matcher after alignment correction.
        Would modify alignment's target blocks
        Added 'query_too' with default as True as it should
        merge query blocks as well - maintain same number of blocks
        """                
        i = 0
        merged = {}
        removed = {}
        blocks = self.blocks[:]
        query_blocks = self.query_blocks[:]
        
        # identify consecutive target blocks
        while i < len(blocks):
            last_merged = None
            for j in range(i, len(blocks)-1):
                if int(blocks[j+1][0]) - int(blocks[j][1]) == 1:
                    last_merged = j + 1
                else:
                    break
                
            if last_merged != None:
                merged[i] = last_merged
                
                for r in range(i + 1, last_merged + 1):
                    removed[r] = True
                    
                i = last_merged + 1
            else:
                i += 1
                        
        if merged and removed:
            # modify block ending coordinate
            for i in range(len(blocks)):
                if merged.has_key(i):
                    blocks[i][1] = blocks[merged[i]][1]
                    if query:
                        query_blocks[i][1] = query_blocks[merged[i]][1]
            
            # remove blocks
            removed_sorted = removed.keys()
            removed_sorted.sort(lambda x,y: y - x)
        
            for i in removed_sorted:
                del blocks[i]
                if query:
                    del query_blocks[i]
            
            self.blocks = blocks
            if query:
                self.query_blocks = query_blocks
                                                        
