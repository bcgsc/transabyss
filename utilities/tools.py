"""
This module provides general functions used in the analysis modules

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
import sys
import re
import grp
from string import maketrans
import intspan

def set_attrs(obj, attrs):
    """Given a dictionary of attribute values, set object's attributes"""
    if attrs:
        for attr,value in attrs.iteritems():
            setattr(obj, attr, value)
            
def reverse_complement(seq):
    """Reverse complements sequence string"""
    complement = maketrans("agtcAGTC", "tcagTCAG")
    return seq[::-1].translate(complement)

def is_homopolymer(seq):
    """Checks if sequence is a homopolymer run"""
    result = True

    if len(seq) < 2:
        result = False
    else:
        for i in range(1, len(seq)):
            if seq[0].upper() != seq[i].upper():
                result = False
                break

    return result

def get_refseq_from_2bit(annodir, genome=None, two_bit_file=None):
    """Extracts sequence from 2-bit file"""
    from two_bit import TwoBitFileCls
    if genome:
        genome_dir = os.path.join(annodir, genome)
        two_bit = genome_dir + '/' + genome + '.2bit'
    elif two_bit_file:
        two_bit = two_bit_file
        
    refseq = None
    if os.path.exists(two_bit):
        refseq = TwoBitFileCls(genome_dir + '/' + genome + '.2bit')
    
    return refseq

def proper_chrom(chrom, genome=None, chrom_proper=None):
    """Returns proper chromosome name
    UCSC-format if available
    """
    
    if chrom_proper:
        if chrom_proper.has_key(chrom):
            chrom = chrom_proper[chrom]
        elif chrom[:3].lower() == 'chr' and chrom_proper.has_key(chrom[3:]):
            chrom = chrom_proper[chrom[3:]]
    
    if not re.match('^(chr|scaffold)', chrom, re.IGNORECASE):
        chrom = 'chr' + chrom
            
    return chrom

def get_splice_motifs(infile):
    """Extracts splice motifs from splice motif file
    returns a dictionary where key=sequence and value=name
    """
    motifs = {}
    if os.path.isfile(infile):
        for line in open(infile, 'r'):
            cols = line.rstrip("\n").split(' ')
            if len(cols) == 2:
                motifs[cols[0]] = cols[1]
            elif len(cols) == 1:
                motifs[cols[0]] = True

    return motifs

def get_chr_number(chrom):
    """Extracts chromosome name from 'chrxxx' string"""
    chr_pat = re.compile(r'chr(\S+)')
    
    m = chr_pat.search(chrom)
                
    if m and m.group:
        return m.group(1)
    else:
        return chrom
        
def ucsc_chroms(genome, annodir):
    """Extracts conversion of UCSC chromosome names
    eg. hg19"""
    conversion_file = os.path.join(annodir, genome, "ucsc_chr.txt")

    conversions = {}
    if os.path.exists(conversion_file):
        for line in open(conversion_file, 'r'):
            chr_from, chr_to = line.rstrip('\n').split()
            conversions[chr_from] = chr_to
            
    return conversions

def output_log(version, argv, params_list, outdir):
    """Outputs log file with run time, host, command and parameters"""
    import datetime, socket
    logfile = outdir + '/LOG'
    log = open(logfile, 'a')
    log.write('[%s]\n' % (str(datetime.datetime.now())))
    log.write('command=%s %s %s\n' % (sys.executable, os.path.abspath(argv[0]), ' '.join(argv[1:])))
    log.write('host=%s\n' % (socket.gethostname()))
    log.write('version=%s\n' % (version))
    #log.write('%s %s %s\n' % (sys.executable, os.path.abspath(argv[0]), ' '.join(argv[1:])))
    for params in params_list:
        for param, value in params.iteritems():
            log.write('%s=%s\n' % (param, value))
    log.write('\n')
    log.close()
    
def set_permissions(directory):
    """Recursively set permission to allow group read and write access"""
    gid = grp.getgrnam('trans-abyss').gr_gid
    for root, dirs, files in os.walk(directory):
        os.chown(root, -1, gid)
        os.chmod(root, 02775)
        
def seq_overlap(seq1, seq2):
    """Overlaps end of seq1 to start of seq2
    Returns number of bases matched"""
    x = min(len(seq1), len(seq2))
    while x > 0:
        if seq1[-x:].upper() == seq2[:x].upper():
            break
        x -= 1
        
    return x
        
def compare_chr(chr1, chr2):
    """For sorting chromosome names ignoring 'chr'"""
    if chr1[:3].lower() == 'chr':
        chr1 = chr1[3:]
    if chr2[:3].lower() == 'chr':
        chr2 = chr2[3:]
    
    if re.match('^\d+$', chr1) and not re.match('^\d+$', chr2):
        return -1
    
    elif not re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
        return 1
    
    else:
        if re.match('^\d+$', chr1) and re.match('^\d+$', chr2):
            chr1 = int(chr1)
            chr2 = int(chr2)
            
        if chr1 < chr2:
            return -1
        elif chr1 > chr2:
            return 1
        else:
            return 0
        
def extract_cytobands(annodir, cytoband_file=None, genome=None):
    """Extracts cytoband information from file
    Returns dictionary where key=chr value=(start, end, band, pos/neg)
    """
    if cytoband_file is None and genome is not None:
        cytoband_file = os.path.join(annodir, genome, "cytoBand.txt")
    data = {}
    if os.path.exists(cytoband_file):
        for line in open(cytoband_file, 'r'):
            cols = line.rstrip('\n').split()
        
            if not data.has_key(cols[0]):
                data[cols[0]] = []         
            data[cols[0]].append([int(cols[1]), int(cols[2]), cols[3], cols[4]])

    return data
            
def overlap_cytobands(cytobands, target1, coord1, target2, coord2):
    """Returns cytobands overlapping given 2 coordinates"""
    band1 = None
    band2 = None
    if cytobands.has_key(target1):
        for band in cytobands[target1]:
            if intspan.overlap([coord1, coord1], band[:2]):
                band1 = band[2]
                break
                
    if cytobands.has_key(target2):
        for band in cytobands[target2]:
            if intspan.overlap([coord2, coord2], band[:2]):
                band2 = band[2]
                break
    
    return band1, band2

def formulate_event(event_type, chr1, chr2, band1, band2):
    """Formulates event in terms of cytobands"""
    chr1 = chr1.replace('chr', '')
    chr2 = chr2.replace('chr', '')
    if event_type == 'translocation':
        return "t(%s;%s)(%s;%s)" % (chr1, chr2, band1, band2)
    
    else:
        if band1 == band2:
            return "%s%s%s" % (event_type[:3].lower(), chr1, band1)
        
        else:
            return "%s%s%s-%s" % (event_type[:3].lower(), chr1, band1, band2)
        
def write_header(columns, out):
    """Writes tab-delimited line of given columns to output"""
    out.write('\t'.join(columns) + '\n')
    
def concat_tsv(files, out_file):
    """Concatenates .tsv files"""
    concat_out = open(out_file, 'w')
                    
    #just for check if header is to be written
    count = 0
    for ff in files:
        infile = open(ff, 'r')
        if count == 0:
            concat_out.writelines(infile.readlines())
        else:
            concat_out.writelines(infile.readlines()[1:])
        infile.close()
                        
        count += 1
                        
    concat_out.close()
    
def remove_from_list(items, indices):
    """Removes items from list given indices"""
    indices.sort(key=int)
    indices.reverse()
    for i in indices:
        del items[i]
                
def substr_search_with_consecutive_mismatches(target, query, seed=2, max_mismatches=1):
    """Searches string with consecutive mismatches
    Used for checking open reading frame in gene fusion"""
    # 2 consecutive matches generate a seed
    start = None
    for i in range(len(target)):
        if target[i:i + seed] == query[:seed]:
            start = i
            break
                        
    if start is not None:
        i = start
        j = 0
        finished = False
        mismatches = 0
        last_match = None                            
        while i < len(target) and j < len(query) and not finished:        
            if target[i] != query[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    finished = True                                
            else:
                last_match = i, j
                mismatches = 0
                                                            
            i += 1
            j += 1
                                            
        return last_match

    else:
        return None
    
