"""
This module provides methods for indexing and overlapping UCSC DGV entries.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser
from utilities.overlap_coord import OverlapCoord
from utilities.tools import reverse_complement

fields = ['bin', 'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
          'thickStart', 'thickEnd', 'itemRgb', 'landmark', 'varType', 'reference', 'pubMedId',
          'method', 'sample']

def parse_line(line):
    """Parses individual line of UCSC dbSNP file"""
    cols = line.rstrip('\n').split('\t')
    
    data = {}
    data['chrom'] = cols[1]
    data['start'] = int(cols[2]) + 1
    data['end'] = int(cols[3])
    data['name'] = cols[4]
    data['strand'] = cols[6]
    data['varType'] = cols[11]
    data['method'] = cols[14]

    data['loss'] = False
    data['gain'] = False
    if cols[9] == '13107200':
	data['loss'] = True
    elif cols[9] == '200':
	data['gain'] = True
    elif cols[9] == '9127187':
	data['gain'] = True
	data['loss'] = True
    
    return data

def index(infile, output):
    """Indexes DGV file"""
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
        
def prepare_overlap(genome):
    """Extracts index info into dictionary"""
    package_dir = "/".join(os.path.abspath(__file__).split("/")[:-3])
    genome_dir = package_dir + '/annotations/' + genome
    
    if not os.path.isdir(genome_dir):
	return None
    
    dgv_file = genome_dir + '/dgv.txt'
    index_file = genome_dir + '/dgv.idx'
    if not os.path.exists(dgv_file) or not os.path.exists(index_file):
	return None

    dgv_overlap = OverlapCoord(dgv_file, index_file)
    dgv_overlap.extract_index()
    
    return dgv_overlap
    
def find_concordants(test, dgv_overlap, buf=20, size_diff=0.01):
    """Finds concordant DGV entries"""
    concordants = []
        
    knowns = dgv_overlap.overlap(test['chrom'], test['start'], test['end'], parse_line=parse_line)
    
    if test['type'] == 'del':
	for known in knowns:
	    if known['varType'] in ('InDel', 'CopyNumber') or 'CNV' in known['varType']:
		if known['loss'] and is_similar(known, test, buf=buf, size_diff=size_diff):
		    concordants.append(known['name'])
		        
    elif test['type'] == 'inv':
	for known in knowns:
	    if known['varType'] in ('Inversion', 'InversionBreakpoint') or 'Inversion' in known['varType']:
		if is_similar(known, test, buf=buf, size_diff=size_diff):
		    concordants.append(known['name'])
		
    return concordants


def is_similar(known, test, buf=0, size_diff=0.0): 
    """Determines if given event is the same as known event given an acceptable size differential"""
    similar = False
    if test['chrom'] == known['chrom']:
	size_test = test['end'] - test['start'] + 1
	# don't do comparison if size_test < 0
	if size_test <= 0:
	    return similar
	
	size_known = known['end'] - known['start'] + 1
	size_diff = abs(size_test - size_known)
	size_diff_fraction = float(size_diff) / float(size_test)
	
	if size_diff_fraction < size_diff and abs(test['start'] - known['start']) <= buf and abs(test['end'] - known['end']) <= buf:
	    similar = True

    return similar
    
if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()
    
    if options.index:
        index(args[0], options.index)
