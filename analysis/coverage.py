"""
This module provides methods for calculating transcript and gene coverage.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import linecache
import re
from utilities.bam import BAM
from utilities import tools, intspan
from annotations import ensembl, knownGene, refGene, aceview, ensg

class Coverage:
    """Calculates and reports transcript and gene coverages"""
    report_fields = ['feature',
                     'model',
                     'transcript',
                     'gene',
                     'exon',
                     'strand',
                     'coord',
                     'feature_size',
                     'bases_reconstructed',
                     'reconstruction',
                     'num_reads',
                     'bases_reads',
                     'depth',
                     'contigs',
                     'num_contigs',
                     'best_contig',
                     'best_contig_reconstruction',
                     'align_blocks',
                     'exons'
                     ]

    def __init__(self, matches={}, models=[], annot_files={}):
        self.matches_by_model = matches
	self.models = models
	self.annot_files = annot_files
	self.results = []
    
    def group_matches(self, matches):
	"""Group matches by gene and transcript"""
	matches_by_transcript = {}
	matches_by_gene = {}
	gene2transcript = {}
	for match in matches:			    
	    if not matches_by_transcript.has_key(match.txt.name):
		matches_by_transcript[match.txt.name] = []
	    matches_by_transcript[match.txt.name].append(match)
	    
	    if match.txt.alias is not None:
		if not matches_by_gene.has_key(match.txt.alias):
		    matches_by_gene[match.txt.alias] = []
		matches_by_gene[match.txt.alias].append(match)
		
		if not gene2transcript.has_key(match.txt.alias):
		    gene2transcript[match.txt.alias] = {}
		gene2transcript[match.txt.alias][match.txt.name] = True
		
	return matches_by_transcript, matches_by_gene, gene2transcript
    
    def get_all_transcripts(self, gene, model, gene2txt):
	"""Extracts all transcripts of the given gene of the given model"""
	all_txts = []
	if gene2txt.has_key(gene):
	    for txt_id, line_num in gene2txt[gene]:
		line = linecache.getline(self.annot_files[model], int(line_num))
		txt = {
		    'e': ensembl.parse_line,
		    'r': refGene.parse_line,
		    'k': knownGene.parse_line,
		    'a': aceview.parse_line,
		    'x': ensg.parse_line,
		    'n': ensembl.parse_line,
		    't': ensembl.parse_line,
		    'g': ensembl.parse_line,
		    }[model](line)
		all_txts.append(txt)
		
	return all_txts
    
    def process(self, contig_coverage=None):	
	self.results = []
	for model in self.models:
	    if not self.matches_by_model.has_key(model):
		continue
	    
	    gene2txt = {
                'e': ensembl.txt2gene,
                'r': refGene.txt2gene,
                'k': knownGene.txt2gene,
                'a': aceview.txt2gene,
	        't': ensembl.txt2gene,
	        'g': ensembl.txt2gene,
                }[model](self.annot_files[model])
            #print len(gene2txt.keys())
	    	    		
	    results = self.process_single_model(self.matches_by_model[model].values(), gene2txt, contig_coverage=contig_coverage)
	    if results:
		self.results.extend(results)
		
    
    def process_single_model(self, matches, gene2txt, contig_coverage=None):
	model = matches[0].model
	matches_by_transcript, matches_by_gene, gene2transcript = self.group_matches(matches)
			    
	results = []
	
	for gene in gene2transcript.keys():
	    matches_gene = matches_by_gene[gene]	    
	    all_txts = self.get_all_transcripts(gene, model, gene2txt)
	    
	    if not all_txts:
		continue
	    
	    collapsed_exons = self.collapse_txts(all_txts)
	    collapsed_blocks = self.collapse_align_blocks([match.align for match in matches_gene])
	    exons_size = intspan.cardinality_multi(collapsed_exons)
	    bases_reconstructed, fraction_reconstructed = self.calc_reconstruction(collapsed_exons, collapsed_blocks, size=exons_size)
	    	    
	    contigs = ','.join([match.align.query for match in matches_gene])
	    
	    if contig_coverage:
		nreads, nbases = self.calc_coverage([match.align.query for match in matches_gene], contig_coverage)
	    else:
		nreads, nbases = '-', '-'
	    	    
	    results.append({'feature': 'gene',
	                    'model': model,
	                    'gene': gene,
	                    'strand': all_txts[0].strand,
	                    'coord': '%s:%s-%s' % (all_txts[0].chrom, collapsed_exons[0][0], collapsed_exons[-1][1]),
	                    'feature_size': exons_size,
	                    'bases_reconstructed': bases_reconstructed,
	                    'reconstruction': "%.3f" % (fraction_reconstructed),
	                    'contigs': contigs,
	                    'num_contigs': len(matches_gene),
	                    'align_blocks': self.blocks_as_string(collapsed_blocks),
	                    'exons': self.blocks_as_string(collapsed_exons),
	                    'num_reads': str(nreads),
	                    'bases_reads': str(nbases),
	                    #'depth': "%.3f" % (float(nbases) / float(exon_len)),
	                    })
	
	    transcripts = gene2transcript[gene].keys()
	    for transcript in sorted(transcripts):	    
		matches_txt = matches_by_transcript[transcript]
		txt = matches_txt[0].txt
	    		
		transcript_len = txt.length
	    
		#best_match
		best_match = None
		for match in matches_txt:
		    if best_match is None or match.coverage > best_match.coverage:
			best_match = match
			
		#reconstruction
		collapsed_blocks = self.collapse_align_blocks([match.align for match in matches_txt])
		exons_size = intspan.cardinality_multi(txt.exons)
		bases_reconstructed, fraction_reconstructed = self.calc_reconstruction(txt.exons, collapsed_blocks, size=exons_size)
		
		contigs = ','.join([match.align.query for match in matches_txt])
	    
		if contig_coverage:
		    nreads, nbases = self.calc_coverage([match.align.query for match in matches_txt], contig_coverage)
		else:
		    nreads, nbases = '-', '-'
	    
		results.append({'feature': 'transcript',
		                'model': txt.model,
		                'transcript': txt.name,
		                'gene': txt.alias,
		                'strand': txt.strand,
		                'coord': '%s:%s-%s' % (txt.chrom, int(txt.txStart) + 1, txt.txEnd),
		                'feature_size': txt.length,
		                'bases_reconstructed': bases_reconstructed,
		                'reconstruction': "%.3f" % (fraction_reconstructed),
		                'contigs': contigs,
		                'num_contigs': len(matches_txt),
		                'best_contig': best_match.align.query,
		                'best_contig_reconstruction': "%.3f" % (best_match.coverage),
		                'align_blocks': self.blocks_as_string(collapsed_blocks),
		                'exons': self.blocks_as_string(txt.exons),
		                'num_reads': str(nreads),
		                'bases_reads': str(nbases),
		                #'depth': "%.3f" % (float(nbases) / float(txt.length)),
		                })
			
	return results
    
    def collapse_align_blocks(self, aligns):
	"""Collapses given list alignments blocks into single list of blocks"""
	all_blocks = []
	# important to copy blocks to all_blocks instead of using it directly!
	for align in aligns:
	    all_blocks.append(align.blocks[:])	
	return self.collapse_blocks(all_blocks)
    
    def collapse_txts(self, txts):
	"""Collpases given list of transcripts into single list of blocks"""
	all_exons = []
	for txt in txts:
	    all_exons.append(txt.exons[:])	
	return self.collapse_blocks(all_exons)
    
    def collapse_blocks(self, blocks):
	"""Collpases given list of blocks into single list of blocks"""
	return intspan.union(blocks)
    
    def calc_reconstruction(self, exons, blocks, size=None):
	"""Calcuates reconstruction given list of alignment blocks and exons"""
	bases_reconstructed = intspan.intersect(exons, blocks)
	if size is None:
	    size = intspan.cardinality_multi(exons)
	fraction_reconstructed = float(bases_reconstructed)/float(size)
	
	return bases_reconstructed, fraction_reconstructed
    
    def calc_coverage(self, contig_coverage, contigs, add_reads=True):
	"""Calculates coverage int terms of number of reads and bases"""
	nreads, nbases = 0, 0
	for contig in contigs:
	    if contig_coverage.has_key(contig) and contig_coverage[contig] is not None:
		if self.add_contig_reads:
		    nreads += contig_coverage[contig][0]
		    nbases += contig_coverage[contig][1]
				
		else:
		    if contig_coverage[contig][0] > nreads:
			nreads = contig_coverage[contig][0]
			nbases = contig_coverage[contig][1]

	return nreads, nbases
    
    def output(self, outfile, no_blocks=False):
	"""Outputs coverage results"""
	out = open(outfile, 'w')
	tools.write_header(self.report_fields, out)
	for info in self.results:
	    cols = []
	    for field in self.report_fields:
		if no_blocks and field in ('exons', 'align_blocks'):
		    continue
		
		if info.has_key(field):
		    cols.append(str(info[field]))
		else:
		    cols.append('-')
		    
	    out.write('\t'.join(cols) + '\n')
	    
    def parse_results(self, coverage_file):
	"""Parses coverage results"""
	records = []
	for line in open(coverage_file):
	    if re.search('feature', line):
		continue
	    cols = line.rstrip('\n').split('\t')
	    record = {}
	    for i in range(len(self.report_fields)):
		record[self.report_fields[i]] = cols[i]
	    records.append(record)
	    
	return records

    def group_results(self, records):
	"""Groups match results by genes and transcripts"""
	grouped = {}
	for record in records:
	    feature = record['feature']
	    model = record['model']
	    gene = record['gene']
	    transcript = record['transcript']
	    
	    if not grouped.has_key(model):
		grouped[model] = {'genes':{}, 'transcripts':{}}
		
	    if feature == 'gene':
		if not grouped[model]['genes'].has_key(gene):
		    grouped[model]['genes'][gene] = []
		grouped[model]['genes'][gene].append(record) 
	    
	    elif feature == 'transcript':
		if not grouped[model]['transcripts'].has_key(gene):
		    grouped[model]['transcripts'][gene] = {}
		if not grouped[model]['transcripts'][gene].has_key(transcript):
		    grouped[model]['transcripts'][gene][transcript] = []
		grouped[model]['transcripts'][gene][transcript].append(record)
	    
	return grouped
    
    def combine_results(self, records, add_reads=False):
	"""Combines match results to calculate coverage"""	
	pool_headers = ['num_reads', 'depth', 'contigs', 'num_contigs', 'best_contig', 'best_contig_reconstruction']
	
	grouped_records = self.group_results(records)
	pooled = {}
	for model in grouped_records.keys():
	    pooled[model] = []
	    
	    genes = grouped_records[model]['genes'].keys()
	    for gene in genes:	
		records = grouped_records[model]['genes'][gene]
		pooled_gene = {'contigs_list': [],
		               'best_contig': None,
		               'best_contig_reconstruction': None,
		               'num_reads': None,
		               'depth': None,
		               'align_blocks': None,
		               'exons': None
		               }
		for field in self.report_fields:
		    if field not in pool_headers:
			pooled_gene[field] = records[0][field]
		
		pooled_transcripts = []
		exons_gene = self.string_as_blocks(records[0]['exons'])
		blocks_gene = []
		
		transcripts = grouped_records[model]['transcripts'][gene].keys()		
		for transcript in transcripts:
		    records = grouped_records[model]['transcripts'][gene][transcript]
		    
		    pooled_transcript = {'contigs_list': [],
		                         'best_contig': None,
		                         'best_contig_reconstruction': None,
		                         'num_reads': None,
		                         'depth': None,
		                         'align_blocks': None,
		                         'exons': None
		                         }
		    for field in self.report_fields:
			if field not in pool_headers:
			    pooled_transcript[field] = records[0][field]
		    
		    exons_transcript = self.string_as_blocks(records[0]['exons'])
		    blocks_transcript = []
		    for record in records:
			# contigs
			contigs = record['contigs'].split(',')
			pooled_transcript['contigs_list'].extend(contigs)
			pooled_gene['contigs_list'].extend(contigs)
			
			# best contig
			if pooled_transcript['best_contig'] is None or float(record['best_contig_reconstruction']) > float(pooled_transcript['best_contig_reconstruction']):
			    pooled_transcript['best_contig'] = record['best_contig']
			    pooled_transcript['best_contig_reconstruction'] = record['best_contig_reconstruction']
			    
			# num_reads
			if pooled_transcript['num_reads'] is None:
			    pooled_transcript['num_reads'] = record['num_reads']
			    pooled_transcript['depth'] = record['depth']
			elif record['num_reads'] != '-' and record['num_reads'] != 'na':
			    pooled_transcript['num_reads'] = int(pooled_transcript['num_reads'])
			    if add_reads:
				pooled_transcript['num_reads'] += int(record['num_reads'])
			    elif pooled_transcript['num_reads'] < int(record['num_reads']):
				pooled_transcript['num_reads'] = int(record['num_reads'])
				
			# reconstruction
			blocks = self.string_as_blocks(record['align_blocks'])
			blocks_transcript.append(blocks)
			blocks_gene.append(blocks)
			
		    # reconstruction transcript
		    union_blocks = intspan.union(blocks_transcript)	    
		    bases_reconstructed = intspan.intersect(exons_transcript, union_blocks)
		    fraction_reconstructed = float(bases_reconstructed)/float(pooled_transcript['feature_size'])
		    
		    # update transcript
		    pooled_transcript['contigs'] = ','.join(pooled_transcript['contigs_list'])
		    pooled_transcript['num_contigs'] = len(pooled_transcript['contigs_list'])
		    if pooled_transcript['num_reads'] != '-' and pooled_transcript['num_reads'] != 'na':
			pooled_transcript['depth'] = "%.3f" % (float(pooled_transcript['num_reads']) / float(pooled_transcript['feature_size']))
		    pooled_transcripts.append(pooled_transcript)
		    pooled_transcript['bases_reconstructed'] = bases_reconstructed
		    pooled_transcript['reconstruction'] = "%.3f" % (fraction_reconstructed)
		    pooled_transcript['align_blocks'] = self.blocks_as_string(union_blocks)
		    
		# update gene	
		for pt in pooled_transcripts:
		    if pooled_gene['best_contig'] is None or float(pt['best_contig_reconstruction']) > float(pooled_gene['best_contig_reconstruction']):
			pooled_gene['best_contig']  = pt['best_contig']
			pooled_gene['best_contig_reconstruction'] = pt['best_contig_reconstruction']
			
		    if pooled_gene['num_reads'] is None:
			pooled_gene['num_reads'] = pt['num_reads'] 
			pooled_gene['depth'] = pt['depth'] 
		    elif pt['num_reads'] != '-' and pt['num_reads'] != 'na':
			pooled_gene['num_reads'] = int(pt['num_reads'])
			if add_reads:
			    pooled_gene['num_reads'] += int(pt['num_reads'])
			elif pooled_gene['num_reads'] < int(pt['num_reads']):
			    pooled_gene['num_reads'] = int(pt['num_reads'])
			    
		# coverage gene
		union_blocks = intspan.union(blocks_gene)
		bases_reconstructed = intspan.intersect(exons_gene, union_blocks)
		fraction_reconstructed = float(bases_reconstructed)/float(pooled_gene['feature_size'])
			    
		pooled_gene['contigs'] = ','.join(pooled_gene['contigs_list'])
		pooled_gene['num_contigs'] = len(pooled_gene['contigs_list'])
		if pooled_gene['num_reads'] != '-' and pooled_gene['num_reads'] != 'na':
		    pooled_gene['depth'] = "%.3f" % (float(pooled_gene['num_reads']) / float(pooled_gene['feature_size']))
		pooled_gene['bases_reconstructed'] = bases_reconstructed
		pooled_gene['reconstruction'] = "%.3f" % (fraction_reconstructed)
		pooled_gene['align_blocks'] = self.blocks_as_string(union_blocks)
		
		pooled[model].append(pooled_gene)
		for pt in pooled_transcripts:
		    pooled[model].append(pt)
		    
	self.results = []
	for model in self.models:
	    if not pooled.has_key(model):
		continue
	    self.results.extend(pooled[model])
	        		        
    @classmethod
    def blocks_as_string(cls, blocks, sep=';'):
	"""Outputs blocks as 'start,end' separated by 'sep'"""
	blocks_str = []
	for block in blocks:
	    blocks_str.append(str(block[0]) + ',' + str(block[1]))
    
	return sep.join(blocks_str)

    @classmethod
    def string_as_blocks(cls, blocks_str, sep=';'):
	"""Converts string of blocks into list""" 
	blocks = []
	for block in blocks_str.split(sep):
	    start, end = block.split(',')
	    blocks.append([int(start), int(end)])
	                  
	return blocks