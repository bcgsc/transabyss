#!/usr/bin/env python

# written by Ka Ming Nip
# Copyright 2014 Canada's Michael Smith Genome Sciences Centre

import argparse
import os
import string
import sys

def eval_blocks(blocksizes, qstarts, tstarts, max_consecutive_edits=1, no_indels=False):
    blocksize_list = blocksizes.split(',')
    qstart_list = qstarts.split(',')
    tstart_list = tstarts.split(',')
    
    #evaluate insertions
    prev_b = None
    prev_q = None
    for b,q in zip(blocksize_list, qstart_list):
        if prev_b is None or prev_q is None:
            prev_b = int(b)
            prev_q = int(q)
        elif b != '' and q != '':
            # number of bases inserted
            n = int(q) - prev_q - prev_b
            if (no_indels and n > 0) or (not no_indels and n > max_consecutive_edits):
                return False
            #endif
            prev_b = int(b)
            prev_q = int(q)
        #endif
    #endfor
    
    #evaluate deletions
    prev_b = None
    prev_t = None
    for b,t in zip(blocksize_list, tstart_list):
        if prev_b is None or prev_t is None:
            prev_b = int(b)
            prev_t = int(t)
        elif b != '' and t != '':
            # number of bases deleted
            n = int(t) - prev_t - prev_b
            if (no_indels and n > 0) or (not no_indels and n > max_consecutive_edits):
                return False
            #endif
            prev_b = int(b)
            prev_t = int(t)
        #endif
    #endfor
    
    return True    
#enddef

def is_redundant(query, refs):
    #print query + ': ' + str(refs)
    # query is NOT redundant if it is the longest sequence
    longest_cid = find_longest(refs)
    return query != longest_cid
#enddef

def find_longest(refs):
    # break ties by chosing the sequence with the lexicographically largest id
    cid, length = sorted(refs.items(), key=lambda x: (x[1],x[0]), reverse=True)[0]
    return cid
#enddef

def extract_cids(psl, samestrand=False, min_percent_identity=0.95, max_consecutive_edits=1, no_indels=False, report_redundant=False):
    queries = set()
    reject_set = set()
    
    prev_query = None
    refs = {}

    with open(psl) as fh:
        for line in fh:
            #match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
            #     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
            
            match, mismatch, repmatch, ns, qgaps, qgapbases, tgaps, tgapbases, strand, qname, qsize, qstart, qend, tname, tsize, tstart, tend, blocks, blocksizes, qstarts, tstarts = line.strip().split('\t')
            
            queries.add(qname)
            
            if qname != tname and \
                float(match)/float(qsize) >= min_percent_identity and \
                int(qgapbases) <= int(qgaps) * max_consecutive_edits and \
                int(tgapbases) <= int(tgaps) * max_consecutive_edits and \
                float(qgapbases) <= float(qsize) * (1.0 - min_percent_identity) and \
                float(tgapbases) <= float(qsize) * (1.0 - min_percent_identity) and \
                (not samestrand or strand == '+') and \
                int(qstart) <= max_consecutive_edits and \
                int(qend) >= int(qsize) - max_consecutive_edits and \
                int(qsize) <= int(tsize) and \
                eval_blocks(blocksizes, qstarts, tstarts, max_consecutive_edits=max_consecutive_edits, no_indels=no_indels):
                
                if prev_query is None:
                    prev_query = qname
                    refs = {qname:int(qsize), tname:int(tsize)}
                elif prev_query == qname:
                    # next 'good' alignment of the same query
                    refs[tname] = int(tsize)
                else:
                    # we have reached the first 'good' alignment of the next query
                    # check whether the previous query is redundant
                    if len(refs) > 0 and is_redundant(prev_query, refs):
                        reject_set.add(prev_query)
                    #endif
                    
                    # reset
                    prev_query = qname
                    refs = {qname:int(qsize), tname:int(tsize)}
                #endif
            #endif
        #endfor
    #endwith
    
    # evaluate the final query
    if len(refs) > 0 and is_redundant(prev_query, refs):
        reject_set.add(prev_query)
    #endif
    
    if report_redundant:
        return reject_set
    #endif
    
    return queries - reject_set
#enddef

def __main__():        
    parser = argparse.ArgumentParser(description='Extract sequence names of redundnat sequences from self-alignments in PSL format.')
    parser.add_argument('psl', metavar='PSL', type=str, help='input headerless PSL file')
    parser.add_argument('--SS', dest='samestrand', help='Do not remove redundant sequences on opposite strands.', action='store_true', default=False)
    parser.add_argument('--min-seq-id', dest='min_seq_id', metavar='FLOAT', type=float, help='Min percentage sequence identity required [%(default)s]', default=0.95)
    parser.add_argument('--max-con-edits', dest='max_con_edits', metavar='INT', type=int, help='Max length of edit allowed [%(default)s bp]', default=1)
    parser.add_argument('--no-indels', dest='no_indels', help='Disallow indels. [indels allowed by default]', action='store_true')
    args = parser.parse_args()

    for cid in extract_redundant(psl=args.psl, require_samestrand=args.samestrand, min_percent_identity=args.min_seq_id, max_consecutive_edits=args.max_con_edits, no_indels=args.no_indels):
        print(cid)
    #endfor
#enddef
    
if __name__ == '__main__':
    __main__()
#endif

#EOF
