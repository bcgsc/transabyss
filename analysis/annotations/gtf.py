"""
This module parses a GTF annotation file and records info into transcript objects.
It also indexes a GTF annotation for line-cache searching.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
from analysis import transcript
from optparse import OptionParser
from ensembl import output

def parse(gtf):
    """Parses GTF file and captures annotation into Transcript objects"""
    txts = {}
    
    for line in open(gtf, 'r'):
        if line[0] == '#':
            continue
        
        cols = line.rstrip('\n').split('\t')
        
        if cols[2] == 'CDS':
            attrs = get_attrs(cols[8])
            cdsStart, cdsEnd = cols[3:5]

            if attrs.has_key('transcript_id'):
                txt_name = attrs['transcript_id']

                if txts.has_key(txt_name):
                    if not txts[txt_name].cdsStart or int(cdsStart) < int(txts[txt_name].cdsStart):
                        txts[txt_name].cdsStart = cdsStart
                    if not txts[txt_name].cdsEnd or int(cdsEnd) > int(txts[txt_name].cdsEnd):
                        txts[txt_name].cdsEnd = cdsEnd

        elif cols[2] == 'exon':
            data = {}
            data['chrom'] = cols[0].lower()
            data['strand'] = cols[6]
            attrs = get_attrs(cols[8])
            if attrs.has_key('gene_name'):
                data['alias'] = attrs['gene_name']
            elif attrs.has_key('gene_id'):
                data['alias'] = attrs['gene_id']
            data['attrs'] = attrs

            start, end = cols[3:5]

            if attrs.has_key('transcript_id'):
                txt_name = attrs['transcript_id']

                if not txts.has_key(txt_name):
                    txt = transcript.Transcript(txt_name)
                    txts[txt_name] = txt
                    txt.exons = []

                    for field, value in data.iteritems():
                        setattr(txt, field, value)

                txts[txt_name].exons.append([int(start), int(end)])

    for txt in txts.values():
        txt.exons.sort(lambda x,y: int(x[0])-int(y[0]))
        txt.txStart = int(txt.exons[0][0])
        txt.txEnd = int(txt.exons[-1][1])

    return txts.values()

def get_attrs(text):
    """Extracts attributes(gene, transcript, etc) from attribute string in GTF file"""
    attrs = {}
    for pair in text.rstrip(';').split(';'):
        name, value = pair.rstrip().lstrip().split(' ')
        value = value.strip('"')
        attrs[name] = value
    return attrs

if __name__ == '__main__':
    usage = "Usage: %prog gtf-file outfile"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")

    (options, args) = parser.parse_args()

    if len(args) == 2:
        output(parse(args[0]), args[1])

