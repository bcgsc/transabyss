"""
This module provides methods for collecting and reporting read support for individual splicing events.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import sys
import os
from utilities.bam import BAM
from event import Event
from transcript import Transcript
from utilities.tools import set_attrs

PACKAGE_DIR = "/".join(os.path.abspath(sys.argv[0]).split("/")[:-2])

def set_read_support(event, genome=False, min_spanning_reads=2, max_coverage_diff=100):
    """Sets read support for event"""
    read_support = True
    
    if event.event_type in ('novel_exon', 'novel_utr', 'AS53'):
        if not event.coverage:
            read_support = False
        else:
            if type(event.coverage).__name__ == 'list':
                contig_coverage = min(event.coverage)
            elif event.coverage == 'na':
                contig_coverage = event.coverage
            else:
                contig_coverage = int(event.coverage)

            if contig_coverage < min_spanning_reads:
                read_support = False
            
    elif event.event_type == 'AS5' or event.event_type == 'AS3':
        if event.spanning_reads < min_spanning_reads:
            read_support = False

    elif event.event_type == 'retained_intron':
        if not event.coverage:
            read_support = False
        else:
            if type(event.coverage).__name__ == 'list':
                contig_coverage = min(event.coverage)
            else:
                contig_coverage = int(event.coverage)

            if contig_coverage < min_spanning_reads:
                read_support = False
                      
    elif event.event_type == 'skipped_exon' or event.event_type == 'read-through':
        # don't check genome coverage because there may be other alleles that don't have skipped exon
        if int(event.spanning_reads) < min_spanning_reads:
            read_support = False

    elif event.event_type == 'novel_intron':
        if int(event.spanning_reads) < min_spanning_reads:
            read_support = False
        
    elif event.event_type == 'novel_transcript':
        valid_links = [l for l in event.coverage if int(l) > min_spanning_reads]
        if len(valid_links) < len(event.coverage):
            read_support = False

    event.read_support = read_support

def find_support_contig(event, contig_bam):
    """Wrapper function for finding read support from reads-to-contig BAM"""
    align = event.align
    target = align.target
    blocks = event.align_blocks
    exons = event.exons
    event.coverage = event.spanning_reads = None
    
    if event.event_type in ('skipped_exon', 'novel_intron'):
        pos = [align.query_blocks[blocks[0]-1][1], align.query_blocks[blocks[1]-1][0]]
        pos.sort(key = int)
        event.contig_coord = '-'.join([str(pos[0]), str(pos[1])])
        spanning_reads = read_support(contig_bam, event.contig, pos, junction=True, seq=True)
        event.spanning_reads = len(spanning_reads)
            
    elif event.event_type == 'read-through':
        if 1 in blocks:
            pos = [align.query_blocks[blocks[-1]-1][1], align.query_blocks[blocks[-1]][0]]
        else:
            pos = [align.query_blocks[blocks[-1]-1][0], align.query_blocks[blocks[-1]-2][1]]
        event.contig_coord = '-'.join([str(pos[0]), str(pos[1])])   
        spanning_reads = read_support(contig_bam, event.contig, pos, junction=True, seq=True)
        event.spanning_reads = len(spanning_reads)
        
    elif event.event_type in ('AS53', 'novel_exon', 'novel_utr'):
        event.spanning_reads = []
        event.coverage = []
        all_spanning_reads = []
                  
        block = blocks[0]
        num_spanning_reads = [0,0]
        contig_coords = ['na', 'na']
        unique_reads = {}
        if block - 2 >= 0:
            pos = [align.query_blocks[block-2][1], align.query_blocks[block-1][0]]
            spanning_reads = read_support(contig_bam, event.contig, pos, junction=True, seq=True)
            num_spanning_reads[0] = len(spanning_reads)
            for read in spanning_reads:
                unique_reads[read.qname + read.seq] = read
            contig_coords[0] = str(pos[0]) + '-' + str(pos[1])
        
        if block <= len(align.query_blocks)-1:
            pos = [align.query_blocks[block-1][1], align.query_blocks[block][0]]
            spanning_reads = read_support(contig_bam, event.contig, pos, junction=True, seq=True)
            num_spanning_reads[1] = len(spanning_reads)
            for read in spanning_reads:
                unique_reads[read.qname + read.seq] = read
            contig_coords[1] = str(pos[0]) + '-' + str(pos[1])
                
        event.spanning_reads = str(num_spanning_reads[0]) + ',' + str(num_spanning_reads[1])
        event.coverage = read_support(contig_bam, event.contig, align.query_blocks[block-1], block=True, seq=False)
        event.contig_coord = ','.join(contig_coords)

    elif event.event_type == 'retained_intron':
        block = blocks[0]
        offset1 = event.txt.exons[exons[0]-1][1] - align.blocks[block-1][0]
        offset2 = align.blocks[block-1][1] - event.txt.exons[exons[1]-1][0]
        if align.query_blocks[block-1][0] < align.query_blocks[block-1][1]:
            pos = [align.query_blocks[block-1][0]+offset1+1, align.query_blocks[block-1][1]-offset2-1]
        else:
            pos = [align.query_blocks[block-1][0]-offset1-1, align.query_blocks[block-1][1]+offset2+1]
        pos.sort(key = int)
        event.contig_coord = '-'.join([str(pos[0]), str(pos[1])])
        spanning_reads = read_support(contig_bam, event.contig, pos, block=True, seq=True)
        event.coverage = len(spanning_reads)

    elif event.event_type == 'AS5' or event.event_type == 'AS3':
        #left
        if (event.event_type == 'AS5' and event.txt.strand == '+') or (event.event_type == 'AS3' and event.txt.strand == '-'):
            if blocks[0] > 1:
                pos = [align.query_blocks[blocks[0]-2][1], align.query_blocks[blocks[0]-1][0]]
            else:
                pos = [align.query_blocks[blocks[0]-1][0], align.query_blocks[blocks[0]-1][0]]
        #right
        else:
            if blocks[0] < len(align.query_blocks):
                pos = [align.query_blocks[blocks[0]-1][1], align.query_blocks[blocks[0]][0]]
            else:
                pos = [align.query_blocks[blocks[0]-1][1], align.query_blocks[blocks[0]-1][1]]

        pos.sort(key = int)
        event.contig_coord = '-'.join([str(pos[0]), str(pos[1])])
        spanning_reads = read_support(contig_bam, event.contig, pos, junction=True, seq=True)
        event.spanning_reads = len(spanning_reads)
        
    elif event.event_type == "novel_transcript":
        spanning_reads = []
        coverage = []
        contig_coords = []
        for i in range(len(align.query_blocks)):
            j = i + 1
            reads = read_support(contig_bam, event.contig, align.query_blocks[i], block=True, seq=True)
            coverage.append(str(len(reads)))        
            contig_coords.append(str(align.query_blocks[i][0]) + '-' + str(align.query_blocks[i][1]))
            if j < len(align.query_blocks):
                pos = [align.query_blocks[i][1], align.query_blocks[j][0]]
                reads = read_support(contig_bam, event.contig, pos, junction=True, seq=True)
                spanning_reads.append(str(len(reads)))
                
        event.spanning_reads = ','.join(spanning_reads)
        event.coverage = ','.join(coverage)
        event.contig_coord = ','.join(contig_coords)

def screen(event):
    """Sets filtering result"""
    if event.read_support:
        event.filter_result = 'passed'
    else:
        event.filter_result = 'failed'
                    
    if event.filter_result == 'passed' and event.event_type in ('AS5', 'AS3', 'novel_exon', 'novel_utr', 'novel_intron', 'AS53', 'novel_transcript') and event.splice and '?' in event.splice:
        event.filter_result = 'failed'

    if event.filter_result == 'passed' and event.event_type == 'retained_intron' and not 'True' in event.multi_3:
        event.filter_result = 'failed'

    elif event.event_type == 'novel_transcript':
        if event.splice and not '?' in event.splice:
            event.filter_result = 'passed'
        else:
            event.filter_result = 'failed'
        
def read_support(contig_bam, contig, pos=None, junction=False, block=None, min_span=4, seq=False):
    """Finding read support from reads-to-contig BAM"""
    reads = []
    if pos is None:
        return reads
    
    pos.sort(lambda x,y: x-y)
    if junction:
        pos_span = [pos[0]-min_span+1, pos[1]+min_span-1]
        reads = contig_bam.confirm(contig, pos_span, feature='breakpoint', seq=seq)
    elif block:
        reads = contig_bam.coverage_simple(contig, pos, total=True, seq=seq)
        
    return reads    
