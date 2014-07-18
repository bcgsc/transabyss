"""
This module provides functions for extracting alignments from PSL files

Author: Readman Chiu rchiu@bcgsc.ca
"""
import sys
import os
import re
import math
import subprocess
from utilities import alignment, intspan, tools
from utilities.tools import ucsc_chroms

fields = {1:"match",
          2:"mismatch",
          3:'repmatch',
          4:'ncount',
          5:'qnuminsert',
          6:'qbaseinsert',
          7:'tnuminsert',
          8:'tbaseinsert',
          9:"query_strand",
          10:"query",
          11:"query_len",
          12:"qstart",
          13:"qend",
          14:"target",
          15:"target_len",
          16:"tstart",
          17:"tend",
          18:"block_count",
          19:"block_sizes",
          20:"qstarts",
          21:"tstarts"
          }

def parse(file, filters=None, noblocks=False, noline=False, minimum=False, splice_motif_file=None, refseq=None, genome=None):
    """Parses a PSL file and returns a list of alignment objects"""
    if is_truncated(file):
	sys.stderr.write("Aborted because of truncated file: %s" % (file))
	sys.exit(100)
	
    aligns = []
    prev_query = None
    prev_line = None
    group = []
    cols = []

    splice_motifs = None
    if splice_motif_file:
        splice_motifs = alignment.get_splice_motifs(splice_motif_file)
	    
    for line in open(file, 'r'):
        line = line.rstrip()

        # double lines seen in blat output?
        if line == prev_line:
            continue

        # check if line begins with number - not headers
        if line and ord(line[0]) >= ord('0') and ord(line[0]) <= ord('9'):
            if not filters:
                aligns.append(create_align(line, noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs))
                continue
            
            cols = line.split("\t")

            if filters and filters.has_key('exclude') and filters['exclude'].has_key(cols[9]):
                continue

            #filtering
            if prev_query and not cols[9] == prev_query:
                indices = screen(group, filters)

                #only extract blocks for filtered set
                if indices:
                    for idx in indices:
                        a = create_align(group[idx], noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs)
                        aligns.append(a)
                    
                del group[:]
                
            group.append(line)

            prev_query = cols[9]
            prev_line = line
            del cols[:]

    #last group
    if group and filters:
        indices = screen(group, filters)
        if indices:
            for idx in indices:
                a = create_align(group[idx], noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs)
                aligns.append(a)

    return aligns

def create_align(line, noblocks=False, noline=False, minimum=False, refseq=None, splice_motifs=None):
    """Creates alignment object from PSL line"""
    cols = line.split("\t")

    align = alignment.Alignment(method="blat")
    for field in fields:
        if fields[field] in ('query', 'target', 'query_strand', 'mismatch'):
            setattr(align, fields[field], cols[field-1])

        if not minimum:
            if fields[field] == "tstart" or fields[field] == 'qstart':
                setattr(align, fields[field], int(cols[field-1])+1)
            elif fields[field] in ('qstarts', 'tstarts', 'block_sizes'):
                continue
            else:
                setattr(align, fields[field], cols[field-1])

    if not minimum:
        align.score = calc_score(cols[0], cols[1], cols[4], cols[6])
        align.match_len = int(cols[0]) + int(cols[1]) +  int(cols[2])
        align.identity = calc_identity(int(cols[11])+1, cols[12], int(cols[15])+1, cols[16], cols[4], cols[1], align.match_len)

    if not noline:
        setattr(align, 'psl_str', line)
        
    if not noblocks:
        align.blocks = get_blocks(cols[20], cols[18])
        align.query_blocks = get_blocks(cols[19], cols[18], strand=cols[8], qsize=cols[10])
	
    # splice sites and orientation
    if align.blocks and refseq:
        align.get_splice_sites(refseq)
        if splice_motifs:
            align.set_orient(splice_motifs)
	        
    return align

def get_blocks(starts, block_sizes, strand=None, qsize=None):
    """Extracts block coordinates"""
    blocks = []
    starts = starts.rstrip(',').split(',')
    blk_sizes = block_sizes.rstrip(',').split(',')
    if strand and qsize and strand == '-':
        for i in range(len(starts)):
            end = int(qsize)-int(starts[i])
            start = end - int(blk_sizes[i]) + 1
            blocks.append([end, start])
    else:
        for i in range(len(blk_sizes)):
            block = [int(starts[i])+1, int(starts[i])+int(blk_sizes[i])]
            blocks.append(block)
        
    return blocks

def screen(group, filters):
    """Filters alignments based on filteres"""
    if not filters:
        return range(len(group))
    
    keep = []
    scores = {}
    ranks = {}
    for idx in range(len(group)):
        cols = group[idx].split("\t")
        scores[idx] = calc_score(cols[0], cols[1], cols[4], cols[6])

    indices = scores.keys()
    indices.sort(lambda x,y: scores[y]-scores[x])

    # rank the alignments
    rank = 1
    for ii in range(len(indices)):
        if 'bestn' in filters and rank > int(filters['bestn']):
            break
        idx = indices[ii]      
        if ii == 0:
            ranks[idx] = 1
        elif scores[idx] == scores[indices[ii-1]]:
            ranks[idx] = ranks[indices[ii-1]]
        else:
            rank += 1
            ranks[idx] = rank

    # skips entire group
    # unique means best alignment only has 1 hit
    if 'unique' in filters and len(indices) > 1 and scores[indices[0]] == scores[indices[1]]:
        return keep
    
    count = 1
    for idx in indices:
        cols = group[idx].split("\t")
    
        if filters.has_key('bestn') and ranks[idx] > int(filters['bestn']):
            break

        if filters.has_key('count') and count > int(filters['count']):
            break
                    
        if filters.has_key('qlen') and int(cols[10]) < int(filters['qlen']):
            break

        match_len = int(cols[0]) + int(cols[1]) +  int(cols[2])
        if filters.has_key('match'):
            if float(match_len)*100/float(cols[10]) < float(filters['match']):
                continue

        if filters.has_key('identity'):
            # if 100% identity, just look at number of mismatches
            if float(filters['identity']) == 100.0 and int(cols[1]) > 0:
                continue
            else: 
                identity = calc_identity(int(cols[11])+1, cols[12], int(cols[15])+1, cols[16], cols[4], cols[1], match_len)
                if identity < float(filters['identity']):
                    continue

        if filters.has_key('ungapped') and int(cols[17]) > 1:
            continue

        count += 1
        keep.append(idx)

    return keep

def calc_score(match, mismatch, qnuminsert, tnuminsert):
    """Calcualtes alignment score using UCSC formula"""
    return int(match) - int(mismatch) - int(qnuminsert) - int(tnuminsert)

def calc_identity(qstart, qend, tstart, tend, qnuminsert, mismatch, match_len):
    """Calcualtes percentage of identity using UCSC formula"""
    q_alisize = int(qend) - int(qstart)
    t_alisize = int(tend) - int(tstart)
    alisize = min(q_alisize, t_alisize)

    if alisize <= 0:
        return 0

    sizediff = q_alisize - t_alisize
    if sizediff < 0:
        sizediff = 0

    insert_factor = int(qnuminsert)
    millibad = (1000 * (float(mismatch) + insert_factor + round(3 * math.log(1 + sizediff)))) / float(match_len)

    return round(100.0 - millibad * 0.1, 1)

def is_truncated(infile):
    """Determines if a PSL file is truncated"""
    if os.path.exists(infile):
	p = subprocess.Popen(["tail", "-1", infile], stdout=subprocess.PIPE)
	last_line = p.communicate()[0]

	num_columns = 21
	if num_columns == len(last_line.split()) or re.match('^-+$', last_line):
	    return False

    return True
    