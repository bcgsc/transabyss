#!/usr/bin/env python

# written by Ka Ming Nip
# updated on June 19, 2014
# Copyright 2014 Canada's Michael Smith Genome Sciences Centre

import argparse
import gzip
import os
import re
import string
import sys


def is_aln_match(cigar):
    return cigar[-1] == 'M' and cigar[0:-1].isdigit()

def is_aln_match_good(cigar, min_percent_identity=0.95, max_consecutive_edits=2, md=None, no_indels=False):
    if cigar == '*':
        return False
    
    if max_consecutive_edits == 0:
        no_indels = True
    
    m = 0
    i = 0
    d = 0
    s = 0
    
    for item in re.findall('[0-9]+[MIDS]', cigar):
        op = item[-1]
        count = int(item[0:-1])
        if op == 'M':
            m += count
        elif op == 'I':
            if no_indels or count > max_consecutive_edits:
                return False
            i += count
        elif op == 'D':
            if no_indels or count > max_consecutive_edits:
                return False
            d += count
        elif op == 'S':
            if count > max_consecutive_edits:
                return False
            s += count
    
    if md is None:
        seq_matches = m
    else:
        seq_matches = 0
        num_consecutive_zeroes = 0
        for item in re.findall('[0-9]+', md):
            value = int(item)
            
            if value == 0:
                num_consecutive_zeroes += 1
            else:
                num_consecutive_zeroes = 0
                seq_matches += value
            
            if num_consecutive_zeroes > max_consecutive_edits:
                return False
    
    return float(seq_matches)/float(m + i + s) >= min_percent_identity

def is_segment_unmapped(flag):
    # 0x4       segment unmapped
    return int(flag) & 0x4

def is_reverse_complemented(flag):
    # 0x10	SEQ being reverse complemented
    return int(flag) & 0x10

def is_redundant(query, refs):       
    # NOT redundant if query is the longest sequence or the sequence with lexicographically largest id when lengths are equal
    (longest_cid, longest_length) = sorted(refs.items(), key=lambda x: (x[1],x[0]), reverse=True)[0]
    #print str(query) + ' ' + str(longest_cid)
    return query != longest_cid

def extract_nr(sam=None, require_samestrand=False, min_percent_identity=0.95, max_consecutive_edits=2, no_indels=False, report_redundant=False):
    f = gzip.open(sam, 'rb')
    
    lengths = {}
    reject_set = set()
    
    current_query = None
    current_length = None
    
    refs = {}
    
    for line in f:
        items = line.rstrip().split()
        name = items[0]

        if name == '@SQ':
            lengths[items[1][3:]] = int(items[2][3:])
        elif len(items) >= 11 and name not in ['@HD', '@SQ', '@RG', '@PG', '@CO'] and name in lengths:
            # This is an alignment line
            if not current_query:
                # the first alignment in SAM
                current_query = name
                current_length = lengths[current_query]
                refs = {current_query : current_length}
            
            md = None
            
            if len(items) > 11:
                for field in items[11].split():
                    partition = field.split(':')
                    if partition[0] == 'MD':
                        md = partition[2]
                        if len(md) == 0:
                            md = None
            
            flag = items[1]

            if name != current_query:
                # Reached the next set of alignments; must look at the current set of alignments
                if is_redundant(current_query, refs):
                    reject_set.add(current_query)
                
                # update everything
                current_query = name
                current_length = lengths[current_query]
                refs = {current_query : current_length}

            if not is_segment_unmapped(flag) and is_aln_match_good(items[5], min_percent_identity, max_consecutive_edits, md, no_indels):
                if not require_samestrand or not is_reverse_complemented(flag):
                    reference = items[2]
                    if reference != current_query:
                        ref_length = lengths[reference]
                        if ref_length >= current_length:
                            refs[reference] = ref_length
    f.close()
    
    # process the final batch of alignments
    if is_redundant(current_query, refs):
        reject_set.add(current_query)

    if report_redundant:
        return reject_set
    
    return set(lengths.keys()) - reject_set

def extract_unmapped(sam=None, require_samestrand=False, min_percent_identity=0.95, max_consecutive_edits=2, no_indels=False):
    cids = set()
    f = gzip.open(sam, 'rb')
    for line in f:
        items = line.rstrip().split()
        name = items[0]

        if len(items) >= 11 and name not in ['@HD', '@SQ', '@RG', '@PG', '@CO']:
            md = None
            if len(items) > 11:            
                for field in items[11].split():
                    partition = field.split(':')
                    if partition[0] == 'MD':
                        md = partition[2]
                        if len(md) == 0:
                            md = None

            flag = items[1]
            unmapped = is_segment_unmapped(flag)
            good_aln = is_aln_match_good(items[5], min_percent_identity, max_consecutive_edits, md, no_indels)
            rc = is_reverse_complemented(flag)
            
            if unmapped or not good_aln or (require_samestrand and rc):
                cids.add(name)
    
    return cids

def extract_mapped(sam=None, require_samestrand=False, min_percent_identity=0.95, max_consecutive_edits=2, no_indels=False):
    cids = set()
    f = gzip.open(sam, 'rb')
    for line in f:
        items = line.rstrip().split()
        name = items[0]

        if len(items) >= 11 and name not in ['@HD', '@SQ', '@RG', '@PG', '@CO']:
            md = None
            if len(items) > 11:            
                for field in items[11].split():
                    partition = field.split(':')
                    if partition[0] == 'MD':
                        md = partition[2]
                        if len(md) == 0:
                            md = None

            flag = items[1]
            unmapped = is_segment_unmapped(flag)
            good_aln = is_aln_match_good(items[5], min_percent_identity, max_consecutive_edits, md, no_indels)
            rc = is_reverse_complemented(flag)
            
            if not unmapped and good_aln and (not require_samestrand or not rc):
                cids.add(name)
                     
    return cids
    
def __main__():        
    parser = argparse.ArgumentParser(description='Extract sequence names from an input SAM file.')
    parser.add_argument('sam', metavar='SAM', type=str, help='input gzip\'d SAM file of alignments')
    parser.add_argument('mode', metavar='MODE', choices=['nr', 'unmapped', 'mapped'], help='choose one of [nr, unmapped, mapped]')
    parser.add_argument('--SS', dest='samestrand', help='Do not remove redundant sequences on opposite strands.', action='store_true', default=False)
    parser.add_argument('--min-seq-id', dest='min_seq_id', metavar='FLOAT', type=float, help='Min percentage sequence identity required [%(default)s]', default=0.95)
    parser.add_argument('--max-con-edits', dest='max_con_edits', metavar='INT', type=int, help='Max length of edit allowed [%(default)s bp]', default=2)
    parser.add_argument('--no-indels', dest='no_indels', help='Disallow indels. [indels allowed by default]', action='store_true')
    args = parser.parse_args()

    if args.mode == 'nr':
        cids = extract_nr(sam=args.sam, require_samestrand=args.samestrand, min_percent_identity=args.min_seq_id, max_consecutive_edits=args.max_con_edits, no_indels=args.no_indels)
        for cid in cids:
            print cid
    elif args.mode == 'unmapped':
        cids = extract_unmapped(sam=args.sam, require_samestrand=args.samestrand, min_percent_identity=args.min_seq_id, max_consecutive_edits=args.max_con_edits, no_indels=args.no_indels)
        for cid in cids:
            print cid
    elif args.mode == 'mapped':
        cids = extract_mapped(sam=args.sam, require_samestrand=args.samestrand, min_percent_identity=args.min_seq_id, max_consecutive_edits=args.max_con_edits, no_indels=args.no_indels)
        for cid in cids:
            print cid
    else:
        print 'Invalid mode: %s' % args.mode
        sys.exit(1)
    
if __name__ == '__main__':
    __main__()

#EOF
