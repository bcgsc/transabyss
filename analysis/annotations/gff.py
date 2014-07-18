"""
This module parses a GFF annotation file and records info into transcript objects.
It also indexes a GFF annotation for line-cache searching.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
import re
from analysis import transcript
from optparse import OptionParser
from ensembl import output

def parse(gff, use_cds=False):
    """Parses GTF file and captures annotation into Transcript objects"""
    txts = {}
    
    for line in open(gff, 'r'):
        if line[0] == '#':
            continue

        cols = line.rstrip('\n').split('\t')
        
        if cols[2] == 'mRNA':
            data = {}
            data['chrom'] = cols[0].lower()
            data['txStart'], data['txEnd'] = cols[3:5]
            data['strand'] = cols[6]
            attrs = get_attrs(cols[8])

            if attrs.has_key('ID'):
                data['name'] = attrs['ID']
                txt = transcript.Transcript(data['name'])
                
            if attrs.has_key('Parent'):
                data['alias'] = attrs['Parent']

            txt.exons = []
            for field, value in data.iteritems():
                setattr(txt, field, value)
            txts[data['name']] = txt

        elif cols[2] == 'CDS':
            if use_cds:
                attrs = get_attrs(cols[8])
                start, end = cols[3:5]
                parent = attrs['Parent']

                for pp in parent.split(','):
                    if txts.has_key(pp):
                        txts[pp].exons.append([int(start), int(end)])
            else:
                attrs = get_attrs(cols[8])
                cdsStart, cdsEnd = cols[3:5]
                parent = attrs['Parent']
            
                for pp in parent.split(','):
                    if txts.has_key(pp):
                        if not txts[pp].cdsStart or int(cdsStart) < int(txts[pp].cdsStart):
                            txts[pp].cdsStart = cdsStart
                        if not txts[pp].cdsEnd or int(cdsEnd) > int(txts[pp].cdsEnd):
                            txts[pp].cdsEnd = cdsEnd

        elif not use_cds and cols[2] == 'exon':
            attrs = get_attrs(cols[8])
            start, end = cols[3:5]
            parent = attrs['Parent']

            for pp in parent.split(','):
                if txts.has_key(pp):
                    txts[pp].exons.append([int(start), int(end)])

    for txt in txts.values():
        txt.exons.sort(lambda x,y: int(x[0])-int(y[0]))
        if not txt.cdsStart:
            txt.cdsStart = txt.exons[0][0]
        if not txt.cdsEnd:
            txt.cdsEnd = txt.exons[-1][1]

    return txts.values()

def get_attrs(text):
    """Extracts attributes(gene, transcript, etc) from attribute string in GFF file"""
    attrs = {}
    for pair in text.split(';'):
        if '=' in pair:
            name, value = pair.split('=')
            attrs[name] = value
    return attrs

if __name__ == '__main__':
    usage = "Usage: %prog gff-file outfile"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--index", dest="index", help="index output file")
    parser.add_option("-c", "--use_cds", dest="use_cds", help="use CDS as exon", action="store_true", default=False)

    (options, args) = parser.parse_args()

    if len(args) == 2:
        output(parse(args[0], use_cds=options.use_cds), args[1])
