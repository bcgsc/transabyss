"""
This module extracts alignments from .psl or .sam files,
and outputs into UCSC track with/without filtering

Author: Readman Chiu rchiu@bcgsc.ca
"""
import sys
import os
import re
import glob
from align_parsers import psl, sam
from optparse import OptionParser
import assembly
from track import Track
from utilities import tools

def main(args, options):
    filters = {}
    if len(args) == 2:
        if options.bestn:
            filters['bestn'] = int(options.bestn)
        if options.identity:
            filters['identity'] = float(options.identity)
        if options.match_percent:
            filters['match'] = float(options.match_percent)
        if options.unique:
            filters['unique'] = True
        if options.query_len:
            filters['qlen'] = int(options.query_len)
        if options.target:
            if os.path.exists(options.target):
                filters['target'] = get_targets(options.target)
            else:
                filters['target'] = [options.target]

        # extracts contigs
        contigs = None
        if options.fasta:
            contigs = get_contigs(options.fasta, options.k)
            
	# extracts reference sequence
        refseq = None
	if options.ref:
	    refseq = tools.get_refseq_from_2bit(options.annodir, options.genome)

        # contig sequences are in batches
        inputs = []
        if options.batch:
            inputs = batch(args[1], args[0])
        else:
            inputs.append(args[0])

        if options.batch:
            if options.out and os.path.exists(options.out):
                os.remove(options.out)
                
        for i in range(len(inputs)):
            infile = inputs[i]
            # for deciding if header needs to be output
            first = False
            if i == 0:
                first = True

	    aligns = parse(infile, args[1], filters)
            if options.out:
                output_aligns(options.out, aligns, options.track_name, contigs=contigs, 
		              append=options.batch, header=first, 
		              genome=options.genome, refseq=refseq, 
		              color=options.color, by_fasta=options.by_fasta, annodir=options.annodir)

	# zip up track
	if options.gzip:
	    os.system('gzip --force --best ' + options.out)
	    
    else:
        parser.error("incorrect number of arguments")

def parse(infile, format, filters=None, splice_motif_file=None, noline=None, refseq=None):   
    """Parse psl or sam file and returns alignment objects"""
    aligns = []
    if format == 'psl':
        aligns = psl.parse(infile, filters, splice_motif_file=splice_motif_file, noline=noline, refseq=refseq)
    elif format == 'sam':
        aligns = sam.parse(infile, filters, splice_motif_file=splice_motif_file, noline=noline, refseq=refseq, header=True, original_order=True)            
    return aligns

def output_aligns(outfile, aligns, track_name=None, contigs=None, append=False, 
                  header=True, genome=None, refseq=None, color=None, by_fasta=False, annodir=None):
    """Outputs alignments to file"""
    ext = os.path.splitext(outfile)[1]
    if append:
        out = open(outfile, 'a')
    else:
        out = open(outfile, 'w')
	
    # desc == name
    track_desc = track_name
    if aligns:
        if track_name and genome:
            Track.ucsc_targets(genome, aligns, annodir)
        
        if header and track_name:
	    header_line = "track name=\"%s\" description=\"%s\" visibility=%d itemRgb=\"On\"" % (track_name, track_desc, 3)
	    if color is not None:
		header_line += ' color=\"%s\"' % (color)
            out.write("%s\n" % (header_line))

        # map contig to contig name
        contig_dict = None
        if contigs:
            contig_dict = dict((contig.num, contig) for contig in contigs)
	    
	# sort alignments by contig order in original input fasta file
	if by_fasta:
	    count = 0
	    ordered = {}
	    for contig in contigs:
		ordered[contig.num] = count
		count += 1
	    aligns.sort(lambda x,y: ordered[x.query] - ordered[y.query])
            
        for align in aligns:
            contig = align.query
            if ' ' in contig:
                contig = align.query.split(" ")[0]
            
            if contig_dict and contig_dict.has_key(contig):
                align.contig = contig_dict[contig]
                                
            if ext == ".gff":
                out.write(align.gff("exon"))

            elif ext == ".psl":
                out.write(align.psl(refseq=refseq, genome=genome))
                
    out.close()

def get_contigs(fasta, kmer=None, lib=None):
    """Extracts contigs from assembly given fasta file
    kmer is given only to update contig's kmer attribute,
    """
    sys.stderr.write("extracting contigs...")

    if kmer:
        ass = assembly.Assembly(lib, k=kmer)
    else:
        ass = assembly.Assembly(lib, k=1)
        
    ass.fasta = fasta
    contigs = ass.get_contigs()
    for contig in contigs:
        if ":" in contig.num:
            k, ctg = contig.num.split(':')
            contig.k = k.lstrip('k')
        elif kmer:
            contig.k = kmer

    sys.stderr.write("done\n")
    return contigs

def batch(ext, indir):
    """Extracts file from given directory with given extension"""
    inputs = glob.glob(os.path.join(indir + '/*.'+'[0-9]*.' + ext))
    return sorted(inputs)

if __name__ == '__main__':
    usage = "Usage: %prog file format [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--annodir", dest="annodir", help="Trans-ABySS annotations directory")
    parser.add_option("-n", "--bestn", dest="bestn", help="best hit filter")
    parser.add_option("-i", "--identity", dest="identity", help="identity filter")
    parser.add_option("-m", "--match", dest="match_percent", help="match length percentage")
    parser.add_option("-u", "--unique", dest="unique", help="unique", action="store_true", default=False)
    parser.add_option("-q", "--query_len", dest="query_len", help="query length filter")
    parser.add_option("-t", "--target", dest="target", help="target e.g. chr1:35422000-35425800, chr13, or file containing targets")
    parser.add_option("-o", "--out", dest="out", help="output file, with proper file extension e.g. x.gff, x.psl, x.snv, or x.gff,x.snv")
    parser.add_option("-k", "--track", dest="track_name", help="track name e.g. \"HS1202 k45\"")
    parser.add_option("-f", "--fasta", dest="fasta", help="fasta file, for extracting contig info")
    parser.add_option("-K", "--k", dest="k", help="kmer value, if not present in contig name")
    parser.add_option("-d", "--batch", dest="batch", help="batch mode", action="store_true", default=False)
    parser.add_option("-g", "--genome", dest="genome", help="genome for UCSC chromosome conversion")
    parser.add_option("-r", "--ref", dest="ref", help="use reference sequence for deducing number of Ns in genomic sequence", action="store_true", default=False)
    parser.add_option("-c", "--color", dest="color", help="RGB color for track e.g. 255,0,0")
    parser.add_option("-z", "--zip", dest="gzip", help="gzip track", action="store_true", default=False)
    parser.add_option("-s", "--by_fasta", dest="by_fasta", help="order output by fasta inputs. Default: original alignment order", action="store_true", default=False)
    
    (options, args) = parser.parse_args()
    main(args, options)
    
