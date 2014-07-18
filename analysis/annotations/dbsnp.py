"""
This module provides methods for indexing and overlapping UCSC dbSNP entries.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser
from utilities.overlap_coord import OverlapCoord
from utilities.tools import reverse_complement

fields = ['bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
          'refNCBI', 'refUCSC', 'observed', 'molType', 'class', 'valid', 'avHet',
          'avHetSE', 'func', 'locType', 'weight']

def parse_line(line):
    """Parses individual line of UCSC dbSNP file"""
    cols = line.rstrip('\n').split('\t')
    
    data = {}
    for i in range(len(cols)):
	if i == len(fields):
	    break
	
        data[fields[i]] = cols[i]
        
    if data['class'] == 'single':
        data['type'] = 'snv'
    # sometimes 'class' is 'named' and '(LARGEDELETION)' is reported instead 
    # of actual allele
    elif data['class'] == 'deletion' or 'deletion' in data['observed'].lower():
        data['type'] = 'del'
    elif data['class'] == 'insertion':
        data['type'] = 'ins'
    else:
        data['type'] = 'NA'
        
    if data['type'] == 'ins':
	data['start'] = int(data['chromStart'])
    else:
	data['start'] = int(data['chromStart']) + 1
    data['end'] = int(data['chromEnd'])
    data['allele'] = {}
    data['size'] = 0
    
    for a in data['observed'].split('/'):
        if a != '-':
	    if data['strand'] == '+':
		data['allele'][a.lower()] = True
	    elif data['strand'] == '-':
		data['allele'][reverse_complement(a).lower()] = True
		
            if data['size'] == 0:
                data['size'] = len(a)

    # make sure deletion size is correct, as sometimes '(LARGEDELETION)'
    # will be put as allele
    if data['type'] == 'del':
	data['size'] = int(data['end']) - int(data['start']) + 1

    if data['observed'] == 'lengthTooLong':
	data = {}
	
    return data

def index(infile, output):
    """Indexes dbSNP file"""
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split("\t")
        start = int(int(cols[2])/1000)
        end = int(int(cols[3])/1000)
        target = cols[1]

        if not 'chr' in target:
            target = 'chr' + target

        for n in range(start,end+1):
            index = ':'.join((target,str(n)))
            value = str(line_num)

            if not indices.has_key(index):
                indices[index] = [value]
            else:
                indices[index].append(value)

        line_num += 1

    index_file = open(output, 'w')
    for index in sorted(indices.keys()):
        index_file.write(' '.join((index, ','.join(indices[index]))) + "\n")
        
def prepare_overlap(genome, chrom):
    """Extracts index info into dictionary"""
    package_dir = "/".join(os.path.abspath(__file__).split("/")[:-3])
    genome_dir = package_dir + '/annotations/' + genome
    if not os.path.isdir(genome_dir):
	return None
    
    snp_file = genome_dir + '/snp/' + chrom + '-snp.txt'
    index_file = snp_file + '.idx'
    if not os.path.exists(snp_file) or not os.path.exists(index_file):
	return None
    snp_overlap = OverlapCoord(snp_file, index_file)
    snp_overlap.extract_index()
    
    return snp_overlap
    
def find_concordants(test, snp_overlap, refseq=None, exact=True, buf=20, size_diff=0.01, target=None):
    """Wrapper function for finding concordant dbSNP entries"""
    concordants = []
        
    knowns = snp_overlap.overlap(test['chrom'], test['start'], test['end'], parse_line=parse_line)
    for known in [k for k in knowns if k['type'].lower() == test['type'].lower() and k['size'] == test['size']]:
	if is_concordant(known, test, refseq, target=target):
	    concordants.append(known['name'])
	    
    if refseq and (test['type'] == 'ins' or test['type'] == 'del'):
	for known in [k for k in knowns if k['class'].lower() == 'in-del']:	    
	    if is_part_of_indel(known, test, refseq, target=target):
		concordants.append(known['name'])

    # not exact
    if not exact and not concordants and test['type'] == 'del':
	for known in [k for k in knowns if k['type'].lower() == test['type'].lower() or (k.has_key('observed') and 'deletion' in k['observed'].lower())]:
	    if is_similar(known, test, buf=buf, size_diff=size_diff):
		concordants.append(known['name'])
		
    return concordants
        
def is_concordant(known, test, refseq=None, target=None):
    """Determines if given event is identical as known event"""
    concordant = False
    if test['type'] == known['type']:
        if test['chrom'] == known['chrom']:
            if test['start'] == known['start'] and test['end'] == known['end']:
                if test['type'] == 'del':
                    concordant = True
                elif known['allele'].has_key(test['allele']) and test['allele'] != 'NA':
                    concordant = True
    
    if not concordant and (test['type'] == 'ins' or test['type'] == 'del') and test['size'] == known['size'] and refseq:
        concordant = check_flanking(known, test, refseq, target=target)
            
    return concordant

def is_similar(known, test, buf=0, size_diff=0.0):
    """Determines if given event is the same as known event given an acceptable size differential"""
    similar = False
    
    if test['chrom'] == known['chrom']:
	size_test = test['end'] - test['start'] + 1
	
	# don't do comparison if size_test < 0
	if size_test <= 0:
	    return similar
	
	# can't rely on 'size' attribute in file
	size_known = known['end'] - known['start'] + 1
	size_diff = abs(size_test - size_known)
	size_diff_fraction = float(size_diff) / float(size_test)
	
	if size_diff_fraction < size_diff and abs(test['start'] - known['start']) <= buf and abs(test['end'] - known['end']) <= buf:
	    similar = True
	    
    return similar

def is_part_of_indel(known, test, refseq=None, target=None):
    """Determines if given event is part of an indel event in dbSNP.
    in-del is a non-zero-length ref sequence replaced by a non-zero-length allele
    this is split into a del and an ins from TA
    so this checks whether the TA events is either
    a) a deletion that matches the in-del deletion
    b) an insertion that matches one of the in-del allele
    """
    concordant = False
        
    if test['type'] == 'del':	
	known_start = int(known['chromStart']) + 1
	
	if test['size'] == len(known['refUCSC']):
	    if test['start'] == known_start and test['end'] == known['end']:
		concordant = True
	    else:
		concordant = check_flanking(known, test, refseq, target=target)
    elif test['type'] == 'ins':	
	known_start = int(known['chromStart']) 
	
	if test['start'] == known_start and known['allele'].has_key(test['allele']) and test['allele'] != 'NA':
	    concordant = True
	else:
	    concordant = check_flanking(known, test, refseq, target=target)		
	
    return concordant
                    
def check_flanking(known, test, refseq, flank_size=1000, target=None):
    """Checks flanking sequences of a given event and a dbSNP event to see if they are the same"""
    same = False
        
    # 'target' = chromosome name in reference sequence file
    # which may be different from UCSC annotations (e.g. hg19)
    if target is None:
	target = test['chrom']
    start = min(known['start'], test['start']) - flank_size
    end = max(known['end'], test['end']) + flank_size
    	    
    if test['type'] == 'ins':
	test_seq = refseq.GetSequence(target, start, test['start']) + test['allele'] + refseq.GetSequence(target, test['start'] + 1, end)
	    
	for allele in known['allele'].keys():
	    known_seq = refseq.GetSequence(target, start, known['start']) + allele + refseq.GetSequence(target, known['start'] + 1, end)
	    if known_seq.lower() == test_seq.lower():
		same = True
		break
		    
    elif test['type'] == 'del':
	test_seq = refseq.GetSequence(target, start, test['start'] - 1) + refseq.GetSequence(target, test['end'] + 1, end)
	known_seq = refseq.GetSequence(target, start, known['start'] - 1) + refseq.GetSequence(target, known['end'] + 1, end)
	if known_seq.lower() == test_seq.lower() and known_seq != '':
	    same = True
		    
    return same

if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()
    if options.index:
        index(args[0], options.index)
