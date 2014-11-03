#!/usr/bin/env python

# written by Ka Ming Nip
# updated on September 29, 2014
# Copyright 2014 Canada's Michael Smith Genome Sciences Centre

import igraph
import operator
import re
from igraph import Graph
from common_utils import log
from operator import itemgetter
from dist_utils import PairedPartners
from dist_utils import parse_dist

KEEP_VERTEX_STATE = 2
VISITED_VERTEX_STATE = 1
DEFAULT_VERTEX_STATE = 0
REMOVE_VERTEX_STATE = -1

GRAPH_ATT_NAME = 'name'

class AdjGraph:    
    def __init__(self, graph, lengths, mkcs, states, distances, cid_indices, strand_specific=False):
        self.graph = graph
        self.lengths = lengths
        self.mkcs = mkcs
        self.states = states
        self.distances = distances
        self.cid_indices = cid_indices
        self.strand_specific = strand_specific
    #enddef
    
    def get_length(self, vid):
        if not self.strand_specific:
            if vid % 2 > 0:
                # odd index
                return self.lengths[(vid-1)/2]
            else:
                # even index
                return self.lengths[vid/2]
            #endif
        #endif
        return self.lengths[vid]
    #enddef

    def get_mkc(self, vid):
        if not self.strand_specific:
            if vid % 2 > 0:
                # odd index
                return self.mkcs[(vid-1)/2]
            else:
                # even index
                return self.mkcs[vid/2]
            #endif
        #endif
        return self.mkcs[vid]
    #enddef
        
    def get_rc_vertex_index(self, vid):
        if not self.strand_specific:
            if vid % 2 > 0:
                # odd index
                return vid - 1
            else:
                # even index
                return vid + 1
            #endif
        #endif
        return None
    #enddef
    
    def get_state(self, vid):
        if not self.strand_specific:
            if vid % 2 > 0:
                # odd index
                return self.states[(vid-1)/2]
            else:
                # even index
                return self.states[vid/2]
            #endif
        #endif
        return self.states[vid]
    #enddef
    
    def set_state(self, vid, state):
        if self.strand_specific:
            self.states[vid] = state
        else:
            if vid % 2 > 0:
                # odd index
                self.states[(vid-1)/2] = state
            else:
                # even index
                self.states[vid/2] = state
            #endif
        #endif
    #enddef
    
    def get_distance(self, eid):
        if not self.strand_specific:
            if eid % 2 > 0:
                # odd index
                return self.distances[(eid-1)/2]
            else:
                # even index
                return self.distances[eid/2]
            #endif
        #endif
        return self.distances[eid]
    #enddef
#endclass

def rc(name):
    cid, sign = get_cid_and_sign(name)
    
    opposite_sign = '-'
    
    if sign == opposite_sign:
        opposite_sign = '+'
    #endif
        
    return str(cid) + opposite_sign
#enddef

def get_cid(name):
    return int(name[:-1])
#enddef

def get_sign(name):
    return name[-1]
#enddef

def get_cid_and_sign(name):
    return get_cid(name), get_sign(name)
#enddef

def parse_adj(adj_file, k, skip_set=set(), strand_specific=False, no_islands=True):
    vertices = []
    
    lengths = []
    mean_kmer_covs = []
    states = []
    cid_indices = {}
        
    edges = []
    ds = []
    
    # default value of 'd' is k-1
    default_d = k-1
    largest_cid = 0
    
    cursor = 0
    
    with open(adj_file, 'r') as fh:
        for line in fh:
            linestripped = line.strip()
            if len(linestripped) > 0:
                info, outs, ins = linestripped.split(';', 2)
                
                out_items = outs.strip().split()
                in_items = ins.strip().split()
                                
                # info of this contig
                cid, l, c = info.strip().split(' ', 2)
                                
                cid_int = int(cid)
                if cid_int > largest_cid:
                    largest_cid = cid_int
                #endif
                
                if no_islands and len(out_items) == 0 and len(in_items) == 0:
                    # we skip this island contig
                    continue
                #endif
                    
                if cid_int in skip_set:
                    # skip this contig
                    continue
                #endif
                
                l = int(l)
                mkc = float(c) / (l-k+1)
                
                # add the "+" orientation of this contig to the graph
                cid_plus = cid + '+'
                vertices.append(cid_plus)
                
                lengths.append(l)
                mean_kmer_covs.append(mkc)
                states.append(DEFAULT_VERTEX_STATE)
                
                # we keep the index of the "+" orientation
                cid_indices[cid_int] = cursor
                
                if not strand_specific:
                    # "-" orientation of this contig
                    cursor += 1
                    cid_minus = cid + '-'
                    vertices.append(cid_minus)
                #endif
                
                cursor += 1
                                
                # successors
                skipped = False
                for name in out_items:                    
                    if len(name) == 0:
                        if skipped:
                            skipped = False
                        #endif
                    else:
                        m = re.match(r"\[d=([-+]?\d+)\]", name)
                        if m is not None:
                            # eg. [d=-12]
                            if skipped:
                                skipped = False
                            else:
                                # correct the 'd' attribute of the previous edge
                                ds[-1] = int(m.group(1))
                            #endif
                        else:
                            scid = get_cid(name)
                            if scid in skip_set or (cid_int != scid and scid in cid_indices):
                                skipped = True
                            else:
                                skipped = False
                                # edge for the plus orientation
                                edges.append((cid_plus, name))
                                ds.append(default_d)
                                
                                if not strand_specific:
                                    # edge for the minus orientation
                                    edges.append((rc(name), cid_minus))
                                #endif
                            #endif
                        #endif
                    #endif
                #endfor

                # predecessors
                skipped = False
                for name in in_items:                    
                    if len(name) == 0:
                        if skipped:
                            skipped = False
                        #endif
                    else:
                        m = re.match(r"\[d=([-+]?\d+)\]", name)
                        if m is not None:
                            # eg. [d=-12]
                            if skipped:
                                skipped = False
                            else:
                                # correct the 'd' attribute of the previous edge
                                ds[-1] = int(m.group(1))
                            #endif
                        else:
                            pcid = get_cid(name)
                            if pcid in skip_set or (cid_int != pcid and pcid in cid_indices):
                                skipped = True
                            else:
                                skipped = False
                                # edge for the plus orientation
                                edges.append((name, cid_plus))
                                ds.append(default_d)
                                
                                if not strand_specific:
                                    # edge for the minus orientation
                                    edges.append((cid_minus, rc(name)))
                                #endif
                            #endif
                        #endif
                    #endif
                #endfor
            #endif
        #endfor
    #endwith
    
    assert len(lengths) == len(states) == len(mean_kmer_covs)
    if strand_specific:
        assert len(edges) == len(ds)
        assert len(vertices) == len(lengths)
    else:
        assert len(edges) == len(ds) * 2
        assert len(vertices) == len(lengths) * 2
    #endif

    g = Graph(directed=True)
    g.add_vertices(vertices)
    g.add_edges(edges)

    return AdjGraph(g, lengths, mean_kmer_covs, states, ds, cid_indices, strand_specific), largest_cid
#enddef

def walk(adj_file, k, path_file, strand_specific=False, cov_gradient=0.05, dist_file=None):
    # parse adj
    adj_graph, largest_cid = parse_adj(adj_file, k, strand_specific=strand_specific)

    graph = adj_graph.graph
    log('ADJ: %d vertices, %d edges' % (graph.vcount(), graph.ecount()))
    
    # Descending order of mean kmer coverage
    vidx_mkc_tuple_list = []
    if strand_specific:
        for vidx, mkc in enumerate(adj_graph.mkcs):
            vidx_mkc_tuple_list.append((vidx, mkc))
        #endfor
    else:
        for vidx, mkc in enumerate(adj_graph.mkcs):
            vidx_mkc_tuple_list.append((vidx*2, mkc))
        #endfor
    #endif    
    vidx_mkc_tuple_list.sort(key=itemgetter(1), reverse=True)
        
    # initialize the path id
    pathid = largest_cid+1
    
    cid_partners_dict = None
    if dist_file:
        cid_partners_dict = parse_dist(dist_file)
    #endif
    
    with open(path_file, 'w') as fh_path:
        for seed_index, mkc in vidx_mkc_tuple_list:            
            # extend a path from a seed
            path = extend_seed(seed_index, mkc, adj_graph, cov_gradient=cov_gradient)
            
            if len(path) > 0 and cid_partners_dict is not None:
                path = extend_path_with_paired_support(path, adj_graph, cid_partners_dict)
            #endif
                                    
            if len(path) > 1:
                path_as_cids = []
                
                for index in path:
                    path_as_cids.append(graph.vs[index][GRAPH_ATT_NAME])
                #endfor
                
                fh_path.write(str(pathid) + '\t' + ' '.join(path_as_cids) + '\n')
                                
                pathid += 1
            #endif
        #endfor
    #endwith
    
    #return the number of paths walked
    return pathid - largest_cid - 1
#enddef    

def unbraid(adj_file, k, path_file, err_cid_file, strand_specific=False, cov_gradient=0.05, length_diff_tolerance=1):
    # parse adj
    adj_graph, largest_cid = parse_adj(adj_file, k, strand_specific=strand_specific)
    
    graph = adj_graph.graph
    log('ADJ: %d vertices, %d edges' % (graph.vcount(), graph.ecount()))
    
    # Descending order of mean kmer coverage
    vidx_mkc_tuple_list = []
    if strand_specific:
        for vidx, mkc in enumerate(adj_graph.mkcs):
            vidx_mkc_tuple_list.append((vidx, mkc))
        #endfor
    else:
        for vidx, mkc in enumerate(adj_graph.mkcs):
            vidx_mkc_tuple_list.append((vidx*2, mkc))
        #endfor
    #endif
    vidx_mkc_tuple_list.sort(key=itemgetter(1), reverse=True)
    
    # initialize the path id
    pathid = largest_cid+1
    
    num_errors = 0
    with open(path_file, 'w') as fh_path:
        for seed_index, mkc in vidx_mkc_tuple_list:             
            # extend a path from a seed
            path = extend_seed(seed_index, mkc, adj_graph, cov_gradient=cov_gradient)
                                    
            if len(path) > 1:
                path_as_cids = []
                
                for index in path:
                    path_as_cids.append(graph.vs[index][GRAPH_ATT_NAME])
                #endfor
                
                fh_path.write(str(pathid) + '\t' + ' '.join(path_as_cids) + '\n')
                
                # walk along the path to remove erroneous branches
                find_erroneous_branches(path, adj_graph, k, length_diff_tolerance=length_diff_tolerance)
                                
                pathid += 1
            #endif
        #endfor
        
        with open(err_cid_file, 'w') as fh_err:
            if strand_specific:
                for idx, val in enumerate(adj_graph.states):
                    if val == REMOVE_VERTEX_STATE:
                        num_errors += 1
                        name = graph.vs[idx][GRAPH_ATT_NAME]
                        cid = get_cid(name)
                        fh_path.write(str(cid) + '\n')
                        fh_err.write(str(cid) + '\n')
                     #endif
                #endfor
            else:
                for idx, val in enumerate(adj_graph.states):
                    if val == REMOVE_VERTEX_STATE:
                        num_errors += 1
                        name = graph.vs[idx*2][GRAPH_ATT_NAME]
                        cid = get_cid(name)
                        fh_path.write(str(cid) + '\n')
                        fh_err.write(str(cid) + '\n')
                     #endif
                #endfor
            #endif
        #endwith
    #endwith
    
    #return the number of paths walked, number vertices marked for removal
    return pathid - largest_cid - 1, num_errors
#enddef    

def extend_upstream(seed_index, adj_graph, mkc=float('NaN'), cov_gradient=0.05):
    path = []
    
    graph = adj_graph.graph
    
    strand_specific = adj_graph.strand_specific
    
    min_mkc = float('NaN')
    max_mkc = float('NaN')
    
    if cov_gradient > 0:
        min_mkc = float(mkc) * cov_gradient
        max_mkc = float(mkc) / cov_gradient
    #endif
    
    last_best_cov = float(mkc)
        
    predecessors = graph.predecessors(seed_index)
    while len(predecessors) > 0:
        # find the predecessor with the largest mean kmer coverage
        best_index = None
        best_cov = None
        
        for index in predecessors:
            cov = adj_graph.get_mkc(index)
            
            if best_cov is None or cov > best_cov:
                best_index = index
                best_cov = cov
            #endif
        #endfor
        
        best_state = adj_graph.get_state(best_index)
        if (
            best_state == KEEP_VERTEX_STATE or
            best_state == REMOVE_VERTEX_STATE or
            (min_mkc != float('NaN') and cov_gradient != 0 and best_cov < min(min_mkc, last_best_cov*cov_gradient)) or
            (max_mkc != float('NaN') and cov_gradient != 0 and best_cov > max(max_mkc, last_best_cov/cov_gradient))
            ):
            # no more extensions
            break
        #endif
        
        last_best_cov = best_cov
                
        # Mark this vertex as a path vertex
        adj_graph.set_state(best_index, KEEP_VERTEX_STATE)
        
        path.append(best_index)
        predecessors = graph.predecessors(best_index)
    #endwhile
        
    # reverse the path because we were walking backward
    path.reverse()
        
    return path
#enddef
    
def extend_downstream(seed_index, adj_graph, mkc=float('NaN'), cov_gradient=0.05):
    path = []
    
    graph = adj_graph.graph
    
    strand_specific = adj_graph.strand_specific
    
    min_mkc = float('NaN')
    max_mkc = float('NaN')
    
    if cov_gradient > 0:
        min_mkc = float(mkc) * cov_gradient
        max_mkc = float(mkc) / cov_gradient
    #endif
    
    last_best_cov = float(mkc)
    
    successors = graph.successors(seed_index)
    while len(successors) > 0:
        # find the successor with the largest mean kmer coverage
        best_index = None
        best_cov = None
        
        for index in successors:
            cov = adj_graph.get_mkc(index)
            
            if best_cov is None or cov > best_cov:
                best_index = index
                best_cov = cov
            #endif
        #endfor
        
        best_state = adj_graph.get_state(best_index)
        if (
            best_state == KEEP_VERTEX_STATE or
            best_state == REMOVE_VERTEX_STATE or
            (min_mkc != float('NaN') and cov_gradient != 0 and best_cov < min(min_mkc, last_best_cov*cov_gradient)) or
            (max_mkc != float('NaN') and cov_gradient != 0 and best_cov > max(max_mkc, last_best_cov/cov_gradient))
            ):
            # no more extensions
            break
        #endif
        
        last_best_cov = best_cov
        
        # Mark this vertex as a path vertex
        adj_graph.set_state(best_index, KEEP_VERTEX_STATE)
        
        path.append(best_index)
        successors = graph.successors(best_index)
    #endwhile
        
    return path
#enddef

def extend_seed(seed_index, mkc, adj_graph, cov_gradient=0.05):
    # the path to be populated
    path = []
    
    if adj_graph.get_state(seed_index) == DEFAULT_VERTEX_STATE:
        # This vertex has not been visited
        
        # Mark this seed vertex as a path vertex
        adj_graph.set_state(seed_index, KEEP_VERTEX_STATE)
        
        # walk upstream from the seed
        path.extend(extend_upstream(seed_index, adj_graph, mkc=mkc, cov_gradient=cov_gradient))
        
        # add the seed to the path
        path.append(seed_index)
        
        # walk downstream from the seed
        path.extend(extend_downstream(seed_index, adj_graph, mkc=mkc, cov_gradient=cov_gradient))
    #endif
        
    return path
#enddef

def find_erroneous_branches(path, adj_graph, k, length_diff_tolerance=1):
    #log('Finding erroneous branches...')
    
    # adjacency graph in igraph's graph data structure
    graph = adj_graph.graph
    
    # longest braid to remove
    max_braid_length = 2*k - 1 + length_diff_tolerance
    
    # longest tip to remove
    min_tip_length = 2*k - 2
        
    # vertex ids
    predecessors = set()
    successors = set()
    
    path_first_member = path[0]
    path_last_member = path[-1]
        
    # get all immediate neighbors of the members of the path,
    # skipping predecessors of the start of the path and
    # successors of the end of the path
    for member in path:
        if member != path_first_member:
            for p in graph.predecessors(member):
                if adj_graph.get_state(p) == DEFAULT_VERTEX_STATE:
                    # This vertex has not been visited
                    predecessors.add(p)
                #endif
            #endfor
        #endif
        
        if member != path_last_member:
            for s in graph.successors(member):
                if adj_graph.get_state(s) == DEFAULT_VERTEX_STATE:
                    # This vertex has not been visited
                    successors.add(s)
                #endif
            #endfor
        #endif
    #endfor
    
    candidates = predecessors | successors
    
    #log('Found %d candidates' % len(candidates))
    
    strand_specific = adj_graph.strand_specific
    
    for c in candidates:
        length = adj_graph.get_length(c)
                
        if length <= max_braid_length:
            # the candidate's successor on the path
            candidate_path_successor = None
            
            # the candidate's predecessor on the path
            candidate_path_predecessor = None
            
            candidate_successors_ok = True
            candidate_predecessors_ok = True
            
            out_degree = 0
            in_degree = 0
            
            for cs in graph.successors(c):
                out_degree += 1
                if cs in path:
                    candidate_path_successor = cs
                elif not cs in successors:
                    candidate_successors_ok = False
                    break
                #endif
            #endfor

            if not candidate_successors_ok:
                continue
            #endif

            for cp in graph.predecessors(c):
                in_degree += 1
                if cp in path:
                    candidate_path_predecessor = cp
                elif not cp in predecessors:
                    candidate_predecessors_ok = False
                    break
                #endif
            #endfor
            
            if not candidate_predecessors_ok:
                continue
            #endif
            
            if candidate_path_successor and candidate_path_predecessor:
                distance_on_path = get_physical_distance(adj_graph, candidate_path_predecessor, candidate_path_successor, path)
                braid_distance = get_distance(adj_graph, candidate_path_predecessor, c) + length + get_distance(adj_graph, c, candidate_path_successor)
                if distance_on_path and distance_on_path >= braid_distance - length_diff_tolerance and distance_on_path <= braid_distance + length_diff_tolerance:
                    adj_graph.set_state(c, REMOVE_VERTEX_STATE)
                #endif
            elif length <= min_tip_length and (out_degree == 0 or in_degree == 0):
                adj_graph.set_state(c, REMOVE_VERTEX_STATE)
            #endif
        #endif
    #endfor    
#enddef

def get_distance(adj_graph, source, sink):
    e = adj_graph.graph.get_eid(source, sink)
    return adj_graph.get_distance(e)
#enddef

def get_physical_distance(adj_graph, start, stop, path):
    # return the physical distance BETWEEN start and stop
    
    distance = 0
    u = None
    begun = False
    graph = adj_graph.graph
    
    for v in path:
        if v == start:
            begun = True
            u = v
        elif begun:
            distance += get_distance(adj_graph, u, v)
            if v == stop:
                return distance
            else:
                distance += adj_graph.get_length(v)
                u = v
            #endif
        #endif
    #endfor
    
    return None
#enddef

def extend_path(members, adj_graph):
    # extend from both ends of this path
    
    extended = []
    
    graph = adj_graph.graph

    # walk upstream
    seed_cid, seed_sign = get_cid_and_sign(members[0])
    seed_index = adj_graph.cid_indices[seed_cid]
    if seed_sign == '-':
        seed_index = adj_graph.get_rc_vertex_index(seed_index)
    #endif
    for vid in extend_upstream(seed_index, adj_graph):
        name = graph.vs[vid][GRAPH_ATT_NAME]
        extended.append(name)
    #endfor
    
    # add original members
    extended.extend(members)
    
    # walk downstream
    seed_cid, seed_sign = get_cid_and_sign(members[-1])
    seed_index = adj_graph.cid_indices[seed_cid]
    if seed_sign == '-':
        seed_index = adj_graph.get_rc_vertex_index(seed_index)
    #endif
    for vid in extend_downstream(seed_index, adj_graph):
        name = graph.vs[vid][GRAPH_ATT_NAME]
        extended.append(name)
    #endfor
    
    return extended
#enddef

def extend_path_with_paired_support(path_members_list, adj_graph, cid_partners_dict):
    # For each direction, extend to a neighbor that has the most read pair support to the path.
    # Do not extend to a neighbor already a member of the given path.
    
    graph = adj_graph.graph
    
    upstream_extension = []
    downstream_extension = []    
    
    path_indexes_set = set(path_members_list)
        
    path_cids_set = set()
    for index in path_members_list:
        path_cids_set.add(graph.vs[index][GRAPH_ATT_NAME])
    #endfor
            
    # Extend in the backward direction.
    backward_index = path_members_list[0]    
    while backward_index is not None:
        
        index_support_dict = {}
                
        predecessors = graph.predecessors(backward_index)
        num_predecessors = len(predecessors)
                
        if num_predecessors > 0:
            for index in predecessors:
                if not index in path_indexes_set and not index in upstream_extension:
                    # Tally the path's read pair support for this predecessor.
                    name = graph.vs[index][GRAPH_ATT_NAME]
                    cid, sign = get_cid_and_sign(name)
                    
                    support = 0
                    
                    if cid in cid_partners_dict:
                        if sign == '+':
                            for partner_name, partner_pairs in cid_partners_dict[cid].outs:                            
                                if partner_name in path_cids_set:
                                    support += partner_pairs
                                #endif
                            #endfor
                        else:
                            for partner_name, partner_pairs in cid_partners_dict[cid].ins:
                                if rc(partner_name) in path_cids_set:
                                    support += partner_pairs
                                #endif
                            #endfor
                        #endif
                    #endif
                    
                    if support > 0:
                        index_support_dict[index] = support
                    #endif
                #endif
            #endfor
        
            if len(index_support_dict) > 0:
                # Find the predecessor with the most read pair support.
                backward_index = sorted(index_support_dict.items(), key=itemgetter(1), reverse=True)[0][0]
            else:
                backward_index = None
            #endif
        #elif num_predecessors == 1:
        #    # Only ONE predecessor
        #    backward_index = predecessors[0]
        else:
            backward_index = None
        #endif
        
        if backward_index is not None:
            if not backward_index in path_indexes_set and \
                not backward_index in upstream_extension:
                # Mark this vertex as visited if not already
                adj_graph.set_state(backward_index, KEEP_VERTEX_STATE)
                
                # Add this vertex to the path
                upstream_extension.append(backward_index)
            else:
                backward_index = None
            #endif
        #endif
    #endwhile
    upstream_extension.reverse()
       
    upstream_extension_set = set(upstream_extension)
    
    # Extend in forward direction
    forward_index = path_members_list[-1]
    while forward_index is not None:
        
        index_support_dict = {}
        
        successors = graph.successors(forward_index)
        num_successors = len(successors)
        
        if num_successors > 0:
            for index in successors:
                # Tally the path's read pair support for this successor.
                name = graph.vs[index][GRAPH_ATT_NAME]
                cid, sign = get_cid_and_sign(name)
                
                support = 0
                
                if cid in cid_partners_dict:
                    if sign == '+':
                        for partner_name, partner_pairs in cid_partners_dict[cid].ins:
                            if partner_name in path_cids_set:
                                support += partner_pairs
                            #endif
                        #endfor
                    else:
                        for partner_name, partner_pairs in cid_partners_dict[cid].outs:
                            if rc(partner_name) in path_cids_set:
                                support += partner_pairs
                            #endif
                        #endfor
                    #endif
                #endif
                
                if support > 0:
                    index_support_dict[index] = support
                #endif
            #endfor
            
            if len(index_support_dict) > 0:
                # Find the successor with the most read pair support.
                forward_index = sorted(index_support_dict.items(), key=itemgetter(1), reverse=True)[0][0]
            else:
                forward_index = None
            #endif
        #elif num_successors == 1:
        #    # Only ONE successor
        #    forward_index = successors[0]
        else:
            forward_index = None
        #endif
        
        if forward_index is not None:
            if not forward_index in path_indexes_set and \
                not forward_index in upstream_extension_set and \
                not forward_index in downstream_extension:
                # Mark this vertex as visited if not already
                adj_graph.set_state(forward_index, KEEP_VERTEX_STATE)
            
                # Add this vertex to the path
                downstream_extension.append(forward_index)
            else:
                forward_index = None
            #endif
        #endif
    #endwhile
    
    # create the new path
    new_path = []
    new_path.extend(upstream_extension)
    new_path.extend(path_members_list)
    new_path.extend(downstream_extension)
    
    return new_path
#enddef

def remove_redundant_paths(rrefs, adj_file, k, braids_file, paths_file, new_remove_file, strand_specific=False):
    # load braid_cids
    braid_cids = set()
    with open(braids_file, 'r') as fh:
        for line in fh:
            linestripped = line.strip()
            if len(linestripped) > 0:
                braid_cids.add(int(linestripped))
            #endif
        #endfor
    #endwith
    
    # load adj but skip the edges and vertices from braids
    adj_graph, largest_cid = parse_adj(adj_file, k, skip_set=braid_cids, strand_specific=strand_specific)
    
    cid_indices = adj_graph.cid_indices
    graph = adj_graph.graph
    
    # load paths
    nr_path_members = {}
        
    rrefs_nonpaths = rrefs.copy()
    
    # mark path contigs    
    with open(paths_file, 'r') as fh_in_path:        
        for line in fh_in_path:
            items = line.strip().split('\t')
            num_items = len(items)
            if num_items == 2:
                pathid = items[0]
                members = items[1].split()
                
                if pathid in rrefs:
                    rrefs_nonpaths.remove(pathid)
                    
                    # redundant path members
                    for name in members:
                        cid = get_cid(name)
                        index = cid_indices[cid]
                        # mark the vertex as 'visited' for now
                        adj_graph.set_state(index, VISITED_VERTEX_STATE)
                    #endfor
                else:
                    # non-redundant path members
                    nr_path_members[pathid] = members
                    for name in members:
                        cid = get_cid(name)
                        index = cid_indices[cid]
                        adj_graph.set_state(index, KEEP_VERTEX_STATE)
                    #endfor
                #endif
            #endif
        #endfor
    #endwith
    
    rrefs_islands = []
    for cid_str in rrefs_nonpaths:
        cid_int = int(cid_str)
        if cid_int in cid_indices:
            # mark the vertex as 'visited' for now
            adj_graph.set_state(cid_indices[cid_int], VISITED_VERTEX_STATE)
        else:
            rrefs_islands.append(cid_int)
        #endif    
    #endfor
    rrefs_nonpaths = None
    
    num_contigs_marked_for_removal = 0
    
    # for each non-redundant path, extend each path until reaching another nr path
    # extensions can change states from 'visited' to 'keep'
    for pathid in sorted(nr_path_members):
        members = nr_path_members[pathid]
        extended = extend_path(members, adj_graph)
        if len(extended) > len(members):
            nr_path_members[pathid] = extended
        #endif
    #endfor
    
    # write the list of braids and members of redundant paths
    with open(new_remove_file, 'w') as fh_out_remove:
        for cid in braid_cids:
            fh_out_remove.write(str(cid) + '\n')
        #endfor
        
        for cid in rrefs_islands:
            fh_out_remove.write(str(cid) + '\n')
        #endfor
        
        rr_cids = set()
        if strand_specific:
            for idx, val in enumerate(adj_graph.states):
                # the 'visited' vertices are from redundant paths
                # they would have been marked as 'keep' instead of 'visited' if they were part of extensions
                if val == VISITED_VERTEX_STATE:
                    name = graph.vs[idx][GRAPH_ATT_NAME]
                    cid = get_cid(name)
                    rr_cids.add(cid)
                #endif
            #endfor
        else:
            for idx, val in enumerate(adj_graph.states):
                # the 'visited' vertices are from redundant paths
                # they would have been marked as 'keep' instead of 'visited' if they were part of extensions
                if val == VISITED_VERTEX_STATE:
                    name = graph.vs[idx*2][GRAPH_ATT_NAME]
                    cid = get_cid(name)
                    rr_cids.add(cid)
                #endif
            #endfor
        #endif
        
        for cid in rr_cids:
            fh_out_remove.write(str(cid) + '\n')
        #endfor
        
        num_contigs_marked_for_removal += len(rr_cids)
        #endfor
    #endwith
    
    return num_contigs_marked_for_removal
#enddef

def has_edges(adj_file):
    with open(adj_file, 'r') as fh:
        for line in fh:
            linestripped = line.strip()
            if len(linestripped) > 0:
                info, outs, ins = linestripped.split(';', 2)
                                
                if len(outs.strip().split()) > 0:
                    return True
                #endif
                
                if len(ins.strip().split()) > 0:
                    return True
                #endif
            #endif
        #endfor
    #endwith
    
    return False
#enddef
    
#EOF
