#!/usr/bin/env python

# written by Ka Ming Nip
# updated on September 26, 2014
# Copyright 2014 Canada's Michael Smith Genome Sciences Centre

class PairedPartners:
    """Data structure to store the tuples of incoming and outgoing partner id and the number of supporting read pairs."""

    def __init__(self, ins, outs):
        self.ins = ins
        self.outs = outs
    #enddef
#endclass

def parse_dist(dist_file):
    """Parse a dist file and return a dictionary of cid to PairedPathners."""

    cid_partners_dict = {}

    with open(dist_file, 'r') as fh:
        linestripped = line.strip()
        if len(linestripped) > 0:
            cid, partners_str = linestripped.split(' ', 1)
            outs_str, ins_str = partners_str.split(';', 1)
                        
            out_tuples_list = []
            
            for p in outs_str.strip().split():
                pcid, distance, pairs, error = p.split(',', 3)
                out_tuples_list.append((pcid, pairs))
            #endfor
            
            in_tuples_list = []
            
            for p in ins_str.strip().split():
                pcid, distance, pairs, error = p.split(',', 3)
                in_tuples_list.append((pcid, pairs))
            #endfor
            
            cid_partners_dict[cid] = PairedPartners(in_tuples_list, out_tuples_list)
        #endif
    #endwith
    
    return cid_partners_dict
#enddef
