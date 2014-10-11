#!/usr/bin/env python

# written by Ka Ming Nip
# updated on July 16, 2014

import glob
import itertools
import math
import operator
import os
import shutil
from common_utils import run_shell_cmd
from common_utils import run_multi_shell_cmds
from utilities import psl_cid_extractor
from utilities import sam_cid_extractor


def concat_fastas(path_prefix_map, cat):
    """Concatenate multiple fasta files together, giving the contigs of each set the specified prefix.
    """
    
    with open(cat, 'w') as fout:
        for prefix, path in path_prefix_map.items():
            line_start = '>' + prefix
            with open(path, 'r') as fin:
                for line in fin:
                    if line[0] == '>':
                        #header                        
                        fout.write(line.replace('>', line_start, 1))
                    else:
                        fout.write(line)
                    #endif
                #endfor
            #endwith
        #endfor
    #endwith
#enddef

def blat_self_align(fasta, outputpsl, percent_id=0.95, max_consecutive_edits=1, min_seq_len=32, threads=1):
    """Align all fasta sequences to each other with BLAT.
    """

    files = bin_by_base_and_length(fasta, 20, fasta)
    
    index = 0
    psls = []
    cmds = []
    
    # bin-to-bin alignments (both inter-bin and intra-bin)
    for a,b in itertools.combinations_with_replacement(files.keys(), 2):
        lower = min(a, b)
        upper = max(a, b)
        #print '%d => %d' % (lower, upper)
        psl = outputpsl + '.' + str(index)
        
        cmd_params = ['blat', '-noHead', '-t=dna', '-q=dna', '-out=psl', '-tileSize=18']
        minscore = int(math.ceil(percent_id * lower))
        if minscore <= 18 * 4:
            cmd_params.append('-oneOff=1')
        #endif
        cmd_params.extend(['-maxGap=%d' % max_consecutive_edits, '-maxIntron=%d' % max_consecutive_edits, '-minScore=%d' % minscore, files[upper], files[lower], psl])
        #run_shell_cmd(' '.join(cmd_params))
        cmds.append(' '.join(cmd_params))
        psls.append(psl)
        index += 1
    #endfor
    
    run_multi_shell_cmds(cmds, max_parallel=threads)
    
    # concatenate PSLs
    num_psls = len(psls)
    assert num_psls > 0
    if num_psls == 1:
        shutil.move(psls[0], outputpsl)
    else:
        with open(outputpsl, 'w') as fout:
            for p in psls:
                with open(p, 'r') as fin:
                    for line in fin:
                        if len(line) > 0:
                            fout.write(line)
                        #endif
                    #endfor
                #endwith
            #endfor
        #endwith
    #endif
    
    # remove the tmp fasta files
    for f in files.values():
        if os.path.isfile(f):
            os.remove(f)
        #endif
    #endfor
    
    # remove the tmp psl files
    for f in psls:
        if os.path.isfile(f):
            os.remove(f)
        #endif
    #endfor
#enddef

def bin_by_base_and_length(fasta, bins, out_prefix):
    
    lengths = []
    seqs = []
    
    currentseq = None
    with open(fasta, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                if currentseq is not None:
                    assert currentseq.header is not None and currentseq.seq is not None
                    # The current fasta sequence is now complete
                    lengths.append(len(currentseq.seq))
                #endif
                currentseq = FastaSeq(header=line.rstrip())
                seqs.append(currentseq)
            else:
                assert currentseq is not None
                
                if currentseq.seq is not None:
                    # concatenate the sequence lines
                    currentseq.seq += line.rstrip()
                else:
                    currentseq.seq = line.rstrip()
                #endif
            #endif
        #endfor     

        # The final fasta seq
        if currentseq is not None:
            assert currentseq.header is not None and currentseq.seq is not None
            # The current fasta sequence is now complete
            lengths.append(len(currentseq.seq))
        #endif
    #endwith
    
    lengths.sort()
    expected_bin_size = sum(lengths)/bins
    
    intervals = []    
    curr_bin_size = 0
    lower = lengths[0]
    upper = lower
    for l in lengths:
        if curr_bin_size < expected_bin_size or l == upper:
            curr_bin_size += l
            upper = l
        else:
            intervals.append((lower, upper+1))
            
            curr_bin_size = l
            lower = l
            upper = lower
        #endif
    #endfor
    intervals.append((lower, upper+1))
    
    files = {}
    fhs = []
    for i in range(len(intervals)):
        outpath = out_prefix + '.' + str(i)
        lower = intervals[i][0]
        
        #print str(intervals[i]) + ' : ' + outpath
        files[lower] = outpath
        fhs.append(open(outpath, 'w'))
    #endfor
        
    # [lower, upper)
    for s in seqs:
        length = len(s.seq)
        for index, (lower, upper) in enumerate(intervals):
            if length >= lower and length < upper:
                fhs[index].write('%s\n%s\n' % (s.header, s.seq))
                break
            #endif
        #endfor
    #endfor
    
    for fh in fhs:
        fh.close()
    #endfor
    
    return files
#enddef

def bin_by_length(fasta, bins, out_prefix):
    
    lengths = []
    seqs = []
    
    currentseq = None
    with open(fasta, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                if currentseq is not None:
                    assert currentseq.header is not None and currentseq.seq is not None
                    # The current fasta sequence is now complete
                    lengths.append(len(currentseq.seq))
                #endif
                currentseq = FastaSeq(header=line.rstrip())
                seqs.append(currentseq)
            else:
                assert currentseq is not None
                
                if currentseq.seq is not None:
                    # concatenate the sequence lines
                    currentseq.seq += line.rstrip()
                else:
                    currentseq.seq = line.rstrip()
                #endif
            #endif
        #endfor     

        # The final fasta seq
        if currentseq is not None:
            assert currentseq.header is not None and currentseq.seq is not None
            # The current fasta sequence is now complete
            lengths.append(len(currentseq.seq))
        #endif
    #endwith
    
    lengths.sort()
    count = len(lengths)    
    expected_bin_size = count/bins
    
    intervals = []
    lower = lengths[0]
    curr_bin_size = 0
    upper = lower
    for l in lengths:
        if curr_bin_size < expected_bin_size or l == upper:
            curr_bin_size += 1
            upper = l
        else:
            intervals.append((lower, upper+1))
            
            lower = l
            upper = lower
            curr_bin_size = 1
        #endif
    #endfor
    intervals.append((lower, upper+1))
    
    files = {}
    fhs = []
    for i in range(len(intervals)):
        outpath = out_prefix + '.' + str(i)
        lower = intervals[i][0]
        
        print str(intervals[i]) + ' : ' + outpath
        files[lower] = outpath
        fhs.append(open(outpath, 'w'))
    #endfor
        
    # [lower, upper)
    for s in seqs:
        length = len(s.seq)
        for index, (lower, upper) in enumerate(intervals):
            if length >= lower and length < upper:
                fhs[index].write('%s\n%s\n' % (s.header, s.seq))
                break
            #endif
        #endfor
    #endfor
    
    for fh in fhs:
        fh.close()
    #endfor
    
    return files
#enddef

def blat_merge_fastas(path_prefix_map, merged_fa, concat_fa=None, concat_fa_selfalign_psl=None, percent_identity=0.95, strand_specific=False, cleanup=False, minoverlap=0, threads=1, indel_size_tolerance=1, min_seq_len=32):
    """Merge fasta files into a single fasta file by removing redundant sequences. Redundancy is determined by BLAT alignments.
    """

    if concat_fa is None:
        concat_fa = merged_fa + '.tmp.concat.fa'
    #endif
    
    if concat_fa_selfalign_psl is None:
        concat_fa_selfalign_psl = merged_fa + '.tmp.concat.psl'
    #endif

    # Concatenate the fastas together and give the contigs of each set a prefix
    concat_fastas(path_prefix_map, concat_fa)
    
    # Self-align concatenated fasta with Bowtie2
    blat_self_align(concat_fa, concat_fa_selfalign_psl, percent_id=percent_identity, max_consecutive_edits=indel_size_tolerance, min_seq_len=min_seq_len, threads=threads)

    # Identify NON-redundant contigs
    nrrefs = psl_cid_extractor.extract_cids(psl=concat_fa_selfalign_psl, samestrand=strand_specific, min_percent_identity=percent_identity, max_consecutive_edits=indel_size_tolerance, report_redundant=False)
    
    tmpfiles = []
    nr_fa_long = merged_fa + '.tmp.long.fa'
    nr_fa_short = None
    if minoverlap > 0:
        nr_fa_short = merged_fa + '.tmp.short.fa'
        tmpfiles.append(nr_fa_short)
    #endif
    
    # Gather the non-contained sequences and split into 2 partitions:
    # 1. shorter than (min overlap + 1)
    # 2. longer than or equal to (min overlap + 1)
    filter_fasta(concat_fa, nr_fa_long, min_length=minoverlap+1, keep_set=nrrefs, fasta_out_st=nr_fa_short)
        
    # overlap-layout the long sequences
    if minoverlap > 0:
        # generate the sequence overlap graph
        overlap_dot = merged_fa + '.tmp.long.dot'
        overlap_cmd_params = ['abyss-overlap', '--threads=%d' % threads, '--min=%d' % minoverlap]
        if strand_specific:
            overlap_cmd_params.append('--SS')
        #endif
        overlap_cmd_params.append(nr_fa_long)
        overlap_cmd_params.append('>' + overlap_dot)
        run_shell_cmd(' '.join(overlap_cmd_params))
        
        # layout contigs using the overlap graph
        layout_path = merged_fa + '.tmp.long.path'
        layout_cmd_params = ['abyss-layout', '--kmer=%d' % (minoverlap+1), '--out=%s' % layout_path]
        if strand_specific:
            layout_cmd_params.append('--SS')
        #endif
        layout_cmd_params.append(overlap_dot)
        run_shell_cmd(' '.join(layout_cmd_params))
        
        # generate fasta for O-L
        overlap_fa = merged_fa + '.incomplete'
        mergecontigs_cmd_params = ['MergeContigs', '--kmer=%d' % (minoverlap+1), '--out=%s' % overlap_fa, nr_fa_long, overlap_dot, layout_path]
        run_shell_cmd(' '.join(mergecontigs_cmd_params))
        
        # append the short sequences to the same fasta
        with open(overlap_fa, 'a') as fout:
            with open(nr_fa_short, 'r') as fin:
                for line in fin:
                    fout.write(line)
                #endfor
            #endwith
        #endwith
        
        shutil.move(overlap_fa, merged_fa)
        tmpfiles.extend([concat_fa, concat_fa_selfalign_psl, nr_fa_long, nr_fa_short, overlap_dot, layout_path])
    else:
        shutil.move(nr_fa_long, merged_fa)
        tmpfiles.extend([concat_fa, concat_fa_selfalign_psl])
    #endif
    
    if cleanup and tmpfiles is not None:
        for t in tmpfiles:
            if t is not None and os.path.isfile(t):
                os.remove(t)
            #endif
        #endfor
    #endif
#enddef

def bowtie2_self_align(fasta, outputsam, threads=1, strand_specific=False, path_strip_sam_seq_qual=None, preset='--sensitive', k=2):
    """Align all fasta sequences to each other with Bowtie2.
    """

    # Build index files for the concatenated fasta
    bt2_index_cmd_params = ['bowtie2-build --quiet', fasta, fasta]
    run_shell_cmd(' '.join(bt2_index_cmd_params))
        
    # Self-align concatenated fasta with Bowtie2
    bt2_align_cmd_params = ['set -euo pipefail && bowtie2']
    
    if strand_specific:
        bt2_align_cmd_params.append('--norc')
    #endif
        
    bt2_align_cmd_params.extend([preset, '-k %d' % k, '--omit-sec-seq --end-to-end -f', '-p %d' % threads, fasta, fasta])
    
    if path_strip_sam_seq_qual:
        bt2_align_cmd_params.append('|' + path_strip_sam_seq_qual)
    #endif
    
    bt2_align_cmd_params.append('|gzip -c >' + outputsam)
    
    run_shell_cmd(' '.join(bt2_align_cmd_params))
#enddef

def bowtie2_merge_fastas(path_prefix_map, merged_fa, concat_fa=None, concat_fa_selfalign_sam=None, bt2_threads=1, percent_identity=0.95, strand_specific=False, path_strip_sam_seq_qual=None, cleanup=False, bowtie2_preset='--sensitive', bowtie2_k=2, threads=1, indel_size_tolerance=1):
    """Merge fasta files into a single fasta file by removing redundant sequences. Redundancy is determined by Bowtie2 alignments.
    """

    tmpfiles = []

    if concat_fa is None:
        concat_fa = merged_fa + '.tmp.concat.fa'
    #endif
        
    if concat_fa_selfalign_sam is None:
        concat_fa_selfalign_sam = concat_fa + '.selfalign.sam.gz'
    #endif

    # Concatenate the fastas together and give the contigs of each set a prefix
    concat_fastas(path_prefix_map, concat_fa)
    
    # Self-align concatenated fasta with Bowtie2
    bowtie2_self_align(concat_fa, concat_fa_selfalign_sam, threads=bt2_threads, strand_specific=strand_specific, path_strip_sam_seq_qual=path_strip_sam_seq_qual, preset=bowtie2_preset)

    # Identify NON-redundant contigs
    nrrefs = cid_extractor.extract_nr(sam=concat_fa_selfalign_sam, require_samestrand=strand_specific, min_percent_identity=percent_identity, max_consecutive_edits=indel_size_tolerance, report_redundant=False)
    
    nr_fa = merged_fa + '.tmp.nr.fa'
        
    # Gather the non-contained sequences
    filter_fasta(concat_fa, nr_fa, keep_set=nrrefs)
    
    shutil.move(nr_fa, merged_fa)
    tmpfiles.extend([concat_fa, concat_fa_selfalign_sam])
    
    if cleanup and tmpfiles is not None:        
        for t in tmpfiles:
            if t is not None and os.path.isfile(t):
                os.remove(t)
            #endif
        #endfor
    #endif
#enddef

def abyssmap_rmdups(in_fa, out_fa, strand_specific=False, cleanup=False, threads=1):
    ids_file = in_fa + '.dup_ids' 

    # run abyssmap
    cmd_params = ['abyss-map', '--dup']
    
    if strand_specific:
        cmd_params.append('--SS')
    #endif
    
    if threads > 1:
        cmd_params.append('--threads=%d' % threads)
    #endif
    
    cmd_params.extend([in_fa, in_fa])
    cmd_params.append('> %s' % ids_file)
    
    run_shell_cmd(' '.join(cmd_params))
    
    cids_set = set()
    with open(ids_file, 'r') as fh:
        for line in fh:
            line_stripped = line.strip()
            if len(line_stripped) > 0:
                cids_set.add(line_stripped)
            #endif
        #endfor
    #endwith
        
    filter_fasta(in_fa, out_fa, remove_set=cids_set)
    
    if cleanup:
        os.remove(ids_file)
    #endif
#enddef

def append_with_prefix(core_fa, additional_fa, prefix=''):
    """Append sequences from additional_fa to core_fa and give each a prefix
    """
    
    with open(core_fa, 'a') as fout:
        with open(additional_fa, 'r') as fin:
            line_start = '>' + prefix
            for line in fin:
                if line[0] == '>':
                    #header                        
                    fout.write(line.replace('>', line_start, 1))
                else:
                    fout.write(line)
                #endif
            #endfor
        #endwith
    #endwith
#enddef

def abyssmap_merge_fastas(path_prefix_map, merged_fa, concat_fa=None, strand_specific=False, cleanup=False, threads=1, iterative=False):
    """Merge fasta files into a single fasta file by removing redundant sequences. Only completely contained sequences with exact match are removed.
    """
    
    tmpfiles = []
    
    if concat_fa is None:
        concat_fa = merged_fa + '.tmp.concat.fa'
    #endif
    
    tmpfiles.append(concat_fa)
    
    if iterative:
        tmp_merged_fa = merged_fa + '.tmp.fa'
    
        # get the file size for each fasta        
        prefix_mem_tuples = []
        for prefix, path in path_prefix_map.iteritems():
            prefix_mem_tuples.append( (prefix, os.path.getsize(path)) )
        #endfor
        
        # Sort by file size in ascending order
        prefix_mem_tuples.sort(key=operator.itemgetter(1))
        
        # Add sequences from one additional fasta file to the merge pool in each iteration
        iteration = 0
        for prefix, filesize in prefix_mem_tuples:
            if os.path.isfile(concat_fa):
                os.remove(concat_fa)
            #endif
            
            if iteration > 0:
                shutil.move(tmp_merged_fa, concat_fa)
            #endif
                
            append_with_prefix(concat_fa, path_prefix_map[prefix], prefix=prefix)
            
            abyssmap_rmdups(concat_fa, tmp_merged_fa, strand_specific=strand_specific, cleanup=cleanup, threads=threads)
            
            iteration += 1
        #endfor
        
        shutil.move(tmp_merged_fa, merged_fa)
        
        tmpfiles.append(tmp_merged_fa)
    else:
        # Concatenate all fastas together and give the contigs of each set a prefix
        concat_fastas(path_prefix_map, concat_fa)
        
        # remove duplicates
        abyssmap_rmdups(concat_fa, merged_fa, strand_specific=strand_specific, cleanup=cleanup, threads=threads)
    #endif
    
    if cleanup and tmpfiles is not None:        
        for t in tmpfiles:
            if t is not None and os.path.isfile(t):
                os.remove(t)
            #endif
        #endfor
    #endif
#enddef

class FastaSeq:
    """Data structure for one FASTA sequence.
    """
    def __init__(self, header=None, seq=None):
        self.header = header
        self.seq = seq
    #enddef
#endclass

def filter_fasta(fasta_in, fasta_out, min_length=0, keep_set=None, fasta_out_st=None, remove_set=None):
    """Filter a FASTA file by selecting contigs by a length cutoff and/or a set of ids of contigs to keep/remove.
    """
    
    currentseq = None
    
    fin = open(fasta_in, 'r')
    fout = open(fasta_out, 'w')
    
    fout_st = None
    if fasta_out_st is not None:
        # FASTA file for short sequences
        fout_st = open(fasta_out_st, 'w')
    #endif
    
    count = 0
    count_st = 0
    
    for line in fin:
        if line[0] == '>':
            if currentseq is not None:
                assert currentseq.header is not None and currentseq.seq is not None
                # The current fasta sequence is now complete
                
                if (keep_set is None or currentseq.header.split()[0][1:] in keep_set) and (remove_set is None or not currentseq.header.split()[0][1:] in remove_set):
                    current_seq_length = len(currentseq.seq)
                    
                    if min_length is None or current_seq_length >= min_length:
                        fout.write('%s\n%s\n' % (currentseq.header, currentseq.seq))
                        count += 1
                    elif fout_st is not None:
                        fout_st.write('%s\n%s\n' % (currentseq.header, currentseq.seq))
                        count_st += 1
                    #endif
                #endif
            #endif
            currentseq = FastaSeq(header=line.rstrip())
        else:
            assert currentseq is not None
            
            if currentseq.seq is not None:
                # concatenate the sequence lines
                currentseq.seq += line.rstrip()
            else:
                currentseq.seq = line.rstrip()
            #endif
        #endif
    #endfor     

    # The final fasta seq
    if currentseq is not None:
        assert currentseq.header is not None and currentseq.seq is not None
        # The current fasta sequence is now complete
        
        if (keep_set is None or currentseq.header.split()[0][1:] in keep_set) and (remove_set is None or not currentseq.header.split()[0][1:] in remove_set):
            current_seq_length = len(currentseq.seq)
            
            if min_length is None or current_seq_length >= min_length:
                fout.write('%s\n%s\n' % (currentseq.header, currentseq.seq))
                count += 1
            elif fout_st is not None:
                fout_st.write('%s\n%s\n' % (currentseq.header, currentseq.seq))
                count_st += 1
            #endif
        #endif
    #endif
                
    fin.close()
    fout.close()
    
    if fout_st is not None:
        fout_st.close()
    #endif
    
    return count, count_st
#enddef

def partition_fasta(fasta_file, outdir, format='seq.%s.fa', max_bases=None, max_contigs=None):
    """Partition a FASTA into smaller chunks based on number of bases and/or number of contigs.
    """
    
    path_format = os.path.join(outdir, format)
    
    for fa in glob.glob(path_format % '*'):
        #This could be a re-run; remove the existing output file.
        os.remove(fa)
    #endfor

    total_num_contigs = 0
    total_num_bases = 0

    currentseq = None
    output_num = 1
    num_contigs = 0
    num_bases = 0

    fin = open(fasta_file, 'r')    
    fout = open(path_format % str(output_num), 'w')
    
    for line in fin:
        if line[0] == '>':
            if currentseq is not None:
                assert currentseq.header is not None and currentseq.seq is not None
                # The current fasta sequence is now complete
                
                # write this fasta seq to file                            
                current_seq_length = len(currentseq.seq)
                
                if (max_contigs is not None and num_contigs + 1 > max_contigs) or (max_bases is not None and num_bases + current_seq_length > max_bases):
                    # adding the current sequence to this batch would exceed the thresholds
                    
                    # this batch is done
                    total_num_contigs += num_contigs
                    total_num_bases += num_bases
                    fout.close()
                    
                    # start the next batch
                    output_num += 1
                    fout = open(path_format % str(output_num), 'w')
                    num_contigs = 0
                    num_bases = 0
                #endif
                
                fout.write('%s\n%s\n' % (currentseq.header, currentseq.seq))
                currentseq = FastaSeq(header=line.rstrip())
                
                num_contigs += 1
                num_bases += current_seq_length
            else:
                currentseq = FastaSeq(header=line.rstrip())
            #endif
        else:
            assert currentseq is not None
            
            if currentseq.seq is not None:
                # concatenate the sequence lines
                currentseq.seq += line.rstrip()
            else:
                currentseq.seq = line.rstrip()
            #endif
        #endif
    #endfor
    fin.close()
    
    # Write the final fasta seq to file
    current_seq_length = len(currentseq.seq)
    
    if (max_contigs is not None and num_contigs + 1 > max_contigs) or (max_bases is not None and num_bases + current_seq_length > max_bases):
        # adding the last sequence to this batch would exceed the thresholds
        
        # this batch is done
        total_num_contigs += num_contigs
        total_num_bases += num_bases
        fout.close()
        
        # start the next batch
        output_num += 1
        fout = open(path_format % str(output_num), 'w')
        num_contigs = 0
        num_bases = 0
    #endif
    
    fout.write('%s\n%s\n' % (currentseq.header, currentseq.seq))
    num_contigs += 1
    num_bases += current_seq_length
    
    total_num_contigs += num_contigs
    total_num_bases += num_bases

    fout.close()
    
    return total_num_contigs, total_num_bases, output_num
#enddef

#EOF
