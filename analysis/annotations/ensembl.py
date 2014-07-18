"""
This module provides methods for parsing and indexing a UCSC ensGene file.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
from optparse import OptionParser
from analysis import transcript
from utilities import tools

# for ensGene.txt from UCSC
fields_a = {1:"name", 2:"chrom", 3:"strand", 4:"txStart", 5:"txEnd", 6:"cdsStart", 7:"cdsEnd", 8:"exonCount", 9:"exonStarts", 10:"exonEnds", 12:"alias"}

# for ensGene_ref.txt created in-house
fields_b = {0:"name", 2:"chrom", 3:"strand", 4:"txStart", 5:"txEnd", 6:"cdsStart", 7:"cdsEnd", 8:"exonCount", 9:"exonStarts", 10:"exonEnds", 16:"alias"}

def set_fields(infile=None, line=None):
    """Defines format-related variables related to the above 2 ensGene file formats"""
    sep = name_field = None
    fields = fields_a
    
    # determine which ensGene file it is
    if infile:
        for l in open(infile, 'r'):
            line = l
            break
    if line:
        if line[:3].lower() != 'ens':
            fields = fields_a
            sep = "\t"
            name_field = 1
        else:
            fields = fields_b
            sep = " "
            name_field = 0

    return sep, name_field, fields

def parse(infile):
    """Parses a UCSC ensGene file"""
    txts = []
    ff = open(infile, 'r')  
    for line in ff.readlines():
        txt = parse_line(line)
        if txt is not None:
            txts.append(txt)
    ff.close()
    
    return txts

def parse_line(line):
    """Parses an individual record of a UCSC knownGenes file"""
    sep, name_field, fields = set_fields(line=line)

    cols = line.rstrip("\n").split(sep)
    if sep and len(cols) > 1:
        txt = transcript.Transcript(cols[name_field])
        for i in range(len(cols)):
            if i in fields:
                if i <= 10 or i == 16 or i == 12:
                    setattr(txt, fields[i], cols[i])

        exonStarts = cols[9].rstrip(',').split(',')
        exonEnds = cols[10].rstrip(',').split(',')
        txt.exons = []
        for e in range(len(exonStarts)):
            txt.exons.append([int(exonStarts[e]) + 1, int(exonEnds[e])])

        # calculates transcript length for coverage
        for exon in txt.exons:
            txt.length += int(exon[1]) - int(exon[0]) + 1

        return txt

    return None

def index(infile, output, genome=None):
    """Generates an index file"""
    sep, name_field, fields = set_fields(infile=infile)
    
    indices = {}
    data_file = os.path.abspath(infile)
    line_num = 1
    for line in open(infile, 'r'):
        cols = line.rstrip().split(sep)
        start = int(int(cols[4]) / 1000)
        end = int(int(cols[5]) / 1000)
        target = cols[2]
        if genome is not None:
            target = tools.proper_chrom(target, genome)
        
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

def output(txts, outfile):
    """Outputs in UCSC ensGene format given a list of transcripts"""
    fields = fields_a
    field_idx = {}
    for idx, field in fields.iteritems():
        if field in ('exonStarts', 'exonEnds', 'exonCount'):
            field_idx[field] = idx

    list_size = int(fields.keys()[-1]) + 1
    out = open(outfile, 'w')
    for i in range(len(txts)):
        txt = txts[i]
        data = []
        for idx in range(list_size):
            data.append('NA')
        
        for idx, field in fields.iteritems():
            try:
                value = getattr(txt, field)
            except AttributeError:
                continue
            else:
                data[idx] = str(value)

        data[0] = str(i)
        data[field_idx['exonStarts']] = ','.join([str(int(i[0]) - 1) for i in txt.exons])
        data[field_idx['exonEnds']] = ','.join([str(i[1]) for i in txt.exons])
        data[field_idx['exonCount']] = str(len(txt.exons))
        
        out.write('\t'.join(data) + '\n')
        
    out.close()
    
def txt2gene(infile):
    """Maps a transcript to gene"""
    genes = {}
    count = 0
    for line in open(infile, 'r'):
        count += 1
        if line[:3].lower() != 'ens':
            fields = fields_a
            sep = "\t"
            name_field = 1
        else:
            fields = fields_b
            sep = " "
            name_field = 0

        cols = line.rstrip('\n').split(sep)
        if sep == '\t':
            txt_col, gene_col = 1, 12
        else:
            txt_col, gene_col = 0, 16
           
        if gene_col > len(cols) - 1:
            continue
        
        if not genes.has_key(cols[gene_col]):
            genes[cols[gene_col]] = []
        genes[cols[gene_col]].append([cols[txt_col], count])
        
    return genes
    
if __name__ == '__main__':
    usage = "Usage: %prog annotation-file"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")
    parser.add_option("-g", "--genome", dest="genome", help="genome")

    (options, args) = parser.parse_args()

    if options.index:
        index(args[0], options.index, genome=options.genome)
