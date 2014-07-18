"""
This module provides functions for extracting alignments from SAM files

Author: Readman Chiu rchiu@bcgsc.ca
"""
from utilities import alignment
from utilities.bam import BAM
from utilities.tools import ucsc_chroms, get_splice_motifs
import re
import sys

fields = {1:"qname",
          2:"flag",
          3:'rname',
          4:'pos',
          5:'mapq',
          6:'cigar',
          7:'rnext',
          8:'pnext',
          9:"tlen",
          10:"seq",
          11:"qual",
          }

cigar_pattern_global = re.compile('\d+[A-Z]')
cigar_pattern_local = re.compile('(\d+)([A-Z])')

def parse(samfile, filters=None, noblocks=False, splice_motif_file=None, minimum=False, noline=None, refseq=None, header=False, original_order=False): 
    """Parses alignments in SAM format, returns alignment objects"""
    splice_motifs = None
    if splice_motif_file:
        splice_motifs = get_splice_motifs(splice_motif_file)
    
    aligns = []
    group = []
    prev_query = None
    
    # parses header for target sizes
    target_sizes = None
    if header:
        target_sizes = parse_header(samfile)

    # body
    f = open(samfile, 'r')
    body = [line.rstrip('\n') for line in f if line.strip() and not line[0] == '@']

    # if original alignment order is to be maintained
    ordered = {}
    if original_order:
        for i in range(len(body)):
            cols = body[i].split('\t')
            key = '-'.join([cols[0], cols[2], cols[3]])
            ordered[key] = i
           
    # in case sorting by query is not done properly
    body.sort()
    
    for line in body: 
        if not is_complete(line):
            sys.exit('incomplete record:%s' % (line.rstrip('\n')))
        
        #no target - skip
        if line.split('\t')[2] == '*':
            continue
        
        #check if line begins with number - not headers
        if not filters:            
            align = create_align(line, noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs, target_sizes=target_sizes)
            if align is not None:
                aligns.append(align)
            continue
            
        cols = line.split("\t")

        if filters and filters.has_key('exclude') and filters['exclude'].has_key(cols[0]):
            continue

        # filtering
        if prev_query and not cols[0] == prev_query:
            indices = screen(group, filters)

            #only extract blocks for filtered set
            if indices:
                for idx in indices:
                    align = create_align(group[idx], noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs, target_sizes=target_sizes)
                    if align is not None:
                        aligns.append(align)
                    
            del group[:]
                
        group.append(line)

        prev_query = cols[0]
        prev_line = line
        del cols[:]

    # last group
    if group and filters:
        indices = screen(group, filters)
        if indices:
            for idx in indices:
                align = create_align(group[idx], noblocks, noline, minimum, refseq=refseq, splice_motifs=splice_motifs, target_sizes=target_sizes)
                if align is not None:
                    aligns.append(align)

    if original_order:
        aligns.sort(lambda x,y: ordered['-'.join([x.query, x.target, x.tstart])] - ordered['-'.join([y.query, y.target, y.tstart])])

    return aligns
            
def parse_header(samfile):
    """Parses header to extract target sizes"""
    f = open(samfile, 'r')
    sizes = {}
    for line in f:
        if '@SQ' in line:
            cols = line.rstrip('\n').split('\t')   
            target = cols[1].split(':')[1]
            target_size = cols[2].split(':')[1]
            sizes[target] = target_size
    f.close()
        
    return sizes
        
def create_align(record, noblocks=False, noline=False, minimum=False, refseq=None, splice_motifs=None, target_sizes=None):
    """Creates alignment object for individual record"""
    cols = record.split("\t")
    
    align = alignment.Alignment(method="sam")
    align.query = cols[0]
    align.target = cols[2]

    if target_sizes is not None and target_sizes.has_key(align.target):
        align.target_size = target_sizes[align.target]
    
    align.tstart = cols[3]
    align.score = cols[4]
    if int(cols[1]) & 16:
        align.query_strand = '-'
    else:
        align.query_strand = '+'
        
    if align.target == '*':
        return None
        
    align.query_len = get_query_len(cols[5])
        
    align.blocks, align.query_blocks = get_blocks_from_cigar(cols[5], int(align.tstart), align.query_strand)
    
    if align.blocks is None or align.query_blocks is None:
        return None
    
    align.qstart = min(align.query_blocks[0][0], align.query_blocks[0][1], align.query_blocks[-1][0], align.query_blocks[-1][1])
    align.qend = max(align.query_blocks[0][0], align.query_blocks[0][1], align.query_blocks[-1][0], align.query_blocks[-1][1])
    align.tend = align.blocks[-1][1]
        
    total = 0
    for i in range(len(align.blocks)):
        tsize = align.blocks[i][1] - align.blocks[i][0] + 1
        if align.query_strand == '+':
            qsize = align.query_blocks[i][1] - align.query_blocks[i][0] + 1
        else:
            qsize = align.query_blocks[i][0] - align.query_blocks[i][1] + 1
            
        total += qsize
    
    # the following is all for UCSC psl track
    align.matches, align.mismatches = parse_MD(record)
    align.num_ins, align.bases_ins, align.num_dels, align.bases_dels, num_skips, bases_skips = get_indels_from_CIGAR(cols[5])
    if align.matches is not None and align.mismatches is not None:
        align.mismatches -= (align.bases_ins + align.bases_dels) 
    # no negative value accepted
    if align.mismatches is not None:
        align.mismatches = max(0, align.mismatches)
    if bases_skips is not None:
        align.bases_dels += bases_skips
    if num_skips is not None:
        align.num_dels += num_skips
    
    align.identity = get_percent_identity(record)
    
    # splice sites and orientation
    if align.blocks and refseq:
        align.get_splice_sites(refseq)
        if splice_motifs:
            align.set_orient(splice_motifs)
        
    return align

def get_num_matches(md_string):
    """Calculates total matches from MD string"""
    total_matches = 0
    matches = re.findall("\d+", md_string)
    for m in matches:
        total_matches += int(m)
            
    return total_matches
    
def get_cigar_ops(cigar_string):
    """Extracts operation and length from CIGAR string"""
    ops_str = cigar_pattern_global.findall(cigar_string)
    ops = []
    for length_op in ops_str:
        m = cigar_pattern_local.match(length_op)
        length, op = m.group(1), m.group(2)         
        ops.append([op, int(length)])
        
    return ops
    
def query_align_bounds(cigar_string):
    """Determines qstart and qend from CIGAR string"""
    ops = get_cigar_ops(cigar_string)
                
    qstart = 1
    i = 0
    while i < len(ops):
        op, length = ops[i]
        if op == 'M':
            break
        #soft or hard clipped
        elif op == 'S' or op == 'H':
            qstart += length
        i += 1
            
    qend = get_query_len(cigar_string)
    i = len(ops) - 1
    while i >= 0:
        op, length = ops[i]
        if op == 'M':
            break
        elif op == 'S' or op == 'H':
            qend -= length
        i -= 1
        
    return qstart, qend

def get_query_len(cigar_string): 
    """Determines query length from CIGAR string"""
    ops = get_cigar_ops(cigar_string)
    
    lens = []
    for op, length in ops:
        if op != 'D' and op != 'N' and op != 'P':
            lens.append(int(length))
    return sum(lens)
    
def get_blocks_from_cigar(cigar_string, pos, strand):  
    """Extracts block coordinates from CIGAR string"""
    ops = get_cigar_ops(cigar_string)
    query_len = get_query_len(cigar_string)
    
    tstart = pos  
    if strand == '+':
        qstart = 1
    else:
        qstart = query_len
        
    tblocks = []
    qblocks = []          
    for i in range(len(ops)):
        op, length = ops[i]
        
        if op == 'D' and i == 0:
            return None, None
        
        if op == 'S' or op == 'H':
            if i == 0:
                if strand == '+':
                    qstart = qstart + length
                else:
                    qstart = qstart - length   
                                        
                continue
        
        tblock = None
        qblock = None
        
        if not tblocks and op != 'M':
            return None, None
        
        #match
        if op == 'M':
            tend = tstart + length - 1
            if strand == '+':
                qend = qstart + length - 1
            else:
                qend = qstart - length + 1

            tblock = [tstart, tend]
            qblock = [qstart, qend]
            
        #intron, skipped reference or deletion in reference
        elif op == 'D' or op == 'N':
            #intron
            if op == 3 and length < 3:
                continue
            
            qend = qend
            tend = tstart + length - 1 
            
        #insertion to reference
        elif op == 'I':
            tend = tend
            if strand == '+':
                qend = qstart + length - 1
            else:
                qend = qstart - length + 1
        
        if tblock:
            tblocks.append(tblock)
        if qblock:
            qblocks.append(qblock)
                    
        tstart = tend + 1
        if strand == '+':
            qstart = qend + 1
        else:
            qstart = qend - 1
                                     
    return tblocks, qblocks
       
def screen(group, filters):
    """Screen alignments in dictionary based on filters given"""
    if not filters:
        return range(len(group))
    
    keep = []
    scores = {}
    ranks = {}
    for idx in range(len(group)):
        cols = group[idx].split('\t')
        scores[idx] = int(cols[4])

    indices = scores.keys()
    indices.sort(lambda x,y: scores[y]-scores[x])

    #rank
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
        cols = group[idx].split('\t')
        query_len = get_query_len(cols[5]) 
        qstart, qend = query_align_bounds(cols[5])
        match_len = qend - qstart + 1
        pc_identity = get_percent_identity(group[idx])
        
        if filters.has_key('bestn') and ranks[idx] > int(filters['bestn']):
            break

        if filters.has_key('count') and count > int(filters['count']):
            break
                    
        if filters.has_key('qlen') and query_len < int(filters['qlen']):
            break

        if filters.has_key('match'):
            if float(match_len) * 100 / float(query_len) < float(filters['match']):
                continue

        if filters.has_key('identity') and pc_identity < float(filters['identity']):
            break
        
        if filters.has_key('ungapped') and not re.match('^\d+M$', cols[5]):
            break

        count += 1
        keep.append(idx)

    return keep
    
def get_percent_identity(record):
    """Determines percent of identity from MD tag"""
    cols = record.split('\t')
    query_len = get_query_len(cols[5]) 
    qstart, qend = query_align_bounds(cols[5])
    match_len = qend - qstart + 1
    
    if len(cols) > 11:
        for i in range(11, len(cols)):
            tag, ttype, value = cols[i].split(':')
            if tag == 'MD':
                num_matches = get_num_matches(value)

                pc_identity = float(num_matches) * 100 / float(match_len)
                return pc_identity
            
    return 0

def get_deletions_from_MD(string):
    """Extracts deletions from MD tag"""
    dels = re.findall('\^\D+', string)
    num_dels = len(dels)
    
    bases_del = 0
    for d in dels:
        bases_del += len(d) - 1
        
    return num_dels, bases_del

def get_indels_from_CIGAR(string):
    """Extracts indels from MD tag"""
    ins = re.findall('(\d+)I', string)
    num_ins = len(ins)
    bases_ins = 0
    for i in ins:
        bases_ins += int(i)
        
    dels = re.findall('(\d+)[D]', string)
    num_dels = len(dels)  
    bases_dels = 0
    for i in dels:
        bases_dels += int(i)
        
    skips = re.findall('(\d+)[N]', string)
    num_skips = len(skips)
    bases_skips = 0
    for i in skips:
        bases_skips += int(i)
        
    return int(num_ins), int(bases_ins), int(num_dels), int(bases_dels), int(num_skips), int(bases_skips)

def parse_MD(record):
    """Extracts number of matches and mismatches from MD tag"""
    cols = record.split('\t')
    num_matches, num_mismatches = None, None
    
    if len(cols) > 11:
        for i in range(11, len(cols)):
            tag, ttype, value = cols[i].split(':')
            if tag == 'MD':
                num_matches = get_num_matches(value)
                #num_dels, bases_del = get_deletions_from_MD(value)
                
            elif tag == 'NM':
                num_mismatches = cols[i].split(':')[2]
       
        if num_matches is not None and num_mismatches is not None:
            return int(num_matches), int(num_mismatches)
    
    return num_matches, num_mismatches

def is_complete(record):
    """Determines if SAM record is complete
    Criteria: whether it has 11 fields
    """
    min_fields = 11
    
    num_fields = 0
    for field in record.split('\t'):
        if field != '':
            num_fields += 1
    
    return num_fields >= min_fields
