"""
This module provides method for annotating events with genes, transcripts, feature
(exon, intron, utr, etc) based on coordinates.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
import ConfigParser
import linecache
import annotations
from annotations import ensembl, knownGene, refGene, aceview, ensg
from utilities.intspan import overlap, subsume, cardinality
from pep_change import translate, pep_change
from utilities import tools

PACKAGE_DIR = "/".join(os.path.abspath(__file__).split("/")[:-2])

class FeatureFinder:
    """Localizes events in genes and transcripts and identifies potential
    effects on coding sequences
    """
    def __init__(self, genome, model):
	self.genome = genome
        self.model = model
        self.index = None
        self.pep = {}
        
        config = ConfigParser.ConfigParser()
        config.read(PACKAGE_DIR + "/configs/model_matcher.cfg")
        self.annot_file = None

        if genome in config.sections():
            for field in config.options(genome):
                if field.lower() == model.lower():
                    self.annot_file = PACKAGE_DIR + "/annotations/" + genome + "/" + config.get(genome, field)

        if self.annot_file and os.path.exists(self.annot_file):
            idx_file = self.annot_file.replace('.txt', '.idx')
            self.extract_index(idx_file)

    def extract_index(self, idx_file):
	"""Extracts annotation index into dictionary"""
        self.index = {}
        for line in open(idx_file, 'r'):
            coord, lines = line.rstrip().split(" ")
            self.index[coord] = lines
	    
    def get_txt_lines(self, coord):
	"""Extracts annotation lines of transcript overlapping given coordinate"""
	if not self.index:
            return None
	
	target, start, end = coord[:]
        tstart = int(start)
        tend = int(end)
        tstart_index = ':'.join((target, str(int(tstart/1000))))
        tend_index = ':'.join((target, str(int(tend/1000))))

        txts = []
        txt_lines = {}
        if self.index.has_key(tstart_index):
            for line in self.index[tstart_index].split(','):
                txt_lines[line] = True
            if tend_index != tstart_index:
            	int(int(tstart)/1000)
            	for coord in range(int(tstart/1000)+1,int(tend/1000)+1):
                    idx = ':'.join((target, str(coord)))
                    if self.index.has_key(idx):
                        for line in self.index[idx].split(','):
                            txt_lines[line] = True
			    
	return txt_lines
	    
    def get_feature_coord(self, coord, feature_type, feature_id, txt_id=None, txt_name=None):
	"""Given coord, feature_type(exon|intron), feature_id(number), transcript_id or transcript_name, 
	finds coordinate(start, end) of given feature
	"""
	if not txt_id and not txt_name:
	    return None
	
	# extracts transcript object
	txt_lines = self.get_txt_lines(coord)
	the_txt = None		    
	for line_num in txt_lines.keys():
            line = linecache.getline(self.annot_file, int(line_num))
            txt = {
                'e': ensembl.parse_line,
                'r': refGene.parse_line,
                'k': knownGene.parse_line,
                'a': aceview.parse_line,
                'x': ensg.parse_line
                }[self.model](line)
	    
	    if txt_id and txt.name == txt_id:
		the_txt = txt
		break
	    elif txt_name and txt.alias == txt_name:
		the_txt = txt
		break
	
	if the_txt:
	    if feature_type == 'intron':
		for i in range(len(the_txt.exons)-1):
		    intron = [int(the_txt.exons[i][1])+1, int(the_txt.exons[i+1][0])-1]
		    if txt.strand == '+':
			intron_id = str(i+1)
		    else:
			intron_id = str(len(txt.exons)-i)
		    if intron_id == feature_id:
			feature_coord = intron    
		    
	    elif feature_type == 'exon':
		for i in range(len(txt.exons)):
		    if txt.strand == '+':
			exon_id = str(i+1)
		    else:
			exon_id = str(len(txt.exons)-i)
		    if exon_id == feature_id:
			feature_coord = txt.exons[i]
		    
	return feature_coord
    
    def get_feature(self, coord, gene_only=False, gene_strand=False, exact=False, refseq=None, variant=None, change=None, txt_obj=False, chrom=None, strand=None, all_overlaps=False):
	"""Finds feature given coordinate.
	Given coordinate, return gene feature
	if gene_only: return a [gene1, gene2]
	if gene_only and gene_strand: return [gene1, strand1, gene2, strand2]
	else return "gene1:txt1:feature1|gene2:txt2:feature2|protein_change"
	where feature = intronX, exonX, utr
	      protein_change only reports when txt1 == txt2
	if coordinate corresponds to single base:
	if gene_only: return a [gene1]
	if gene_only and gene_strand: return [gene1, strand1]
	else return "gene1:txt1:feature1|protein_change"
	"""
	features = []
	overlaps = []
	target, start, end = coord.split()
	# go through both end points of coordinate given
	for base in (start, end):
	    if len(features) == 1 and start == end:
		break
	    	   
	    # the 'best' candidate transcript and feature
	    the_txt, the_feature = None, None
	    	    
	    # go through every overlapping transcript
	    txt_lines = self.get_txt_lines([target, base, base])
	    for line_num in txt_lines.keys():
		line = linecache.getline(self.annot_file, int(line_num))
		txt = {
		    'e': ensembl.parse_line,
		    'r': refGene.parse_line,
		    'k': knownGene.parse_line,
		    'a': aceview.parse_line,
		    'x': ensg.parse_line,
		    'n': ensembl.parse_line,
		    't': ensembl.parse_line,
		    'g': ensembl.parse_line,
		    }[self.model](line)

		if not overlap([txt.txStart, txt.txEnd], [base, base]):
		    continue
		if strand and txt.strand and txt.strand != strand:
		    continue
    
		feature = self.identify_feature(start, end, txt, exact=exact)
		if feature:
		    # keep all transcripts and features if 'all_overlaps' is True
		    if all_overlaps:
			ff = ':'.join((txt.alias, txt.name, feature))
			overlaps.append((ff, txt))
					
		    if features and features[0][0] != None:			    
			if txt.name == features[0][0].name:
			    the_txt = txt
			    the_feature = ':'.join((txt.alias, txt.name, feature))
		    
		    # updates best candidate if
		    # - best transcript not defined
		    # - best feature not defined
		    # - new feature is exonic but best candidate isn't
		    # - current candidate is exonic but current candidate's CDS is longer
		    # - current candidate is not exonic but current candidate's CDS is longer
		    elif not the_txt or \
		         ((not the_feature or not self.is_exon(the_feature)) and self.is_exon(feature)) or\
		         (the_txt and self.is_exon(feature) and int(the_txt.cds_length()) < int(txt.cds_length())) or\
		         (the_txt and int(the_txt.cds_length()) < int(txt.cds_length())):
			the_txt = txt
			the_feature = ':'.join((txt.alias, txt.name, feature))
					
	    # novel feature of known gene
	    if the_txt and not the_feature:
		if the_txt:
		    the_feature = ':'.join((the_txt.alias, the_txt.name, 'novel'))
		else:
		    the_feature = ':'.join(('NA', 'NA'))
		    
	    if the_feature:
		if gene_only:		    
		    if gene_strand:
			if the_txt:
			    strand = the_txt.strand
			else:
			    strand = 'NA'			    
			features.append([the_feature.split(':')[0], strand])
		    else:
			features.append([the_feature.split(':')[0]])					    
		else:
		    if the_txt:
			txt_name = the_txt.name
		    else:
			txt_name = 'NA'		    
		    
		    features.append([the_txt, the_feature])
	    else:
		if gene_only:		    
		    if gene_strand:
			features.append(['NA', 'NA'])
		    else:
			features.append(['NA'])
		else:
		    features.append([None, 'NA:NA:NA'])
		    	
	if all_overlaps:
	    return overlaps
		
	if gene_only:
	    if txt_obj:
		return features[0][0], features[0][1], the_txt
	    else:
		return features[0][0], features[0][1]
	    
	else:
	    pepchange = 'NA'
	    
	    # chromosome name for extracting sequence - maybe be different from target
	    if chrom is not None:
		target = chrom

	    if len(features) == 1 or (features[0][0] != None and features[1][0] != None and features[0][0].name == features[1][0].name):
		if 'exon' in features[0][1] and not gene_only and refseq:
		    pepchange = self.find_pep_change([target, start, end], features[0][0], refseq, variant, change)
	    
	    if len(features) > 1:
		feature = '|'.join([features[0][1], features[1][1], pepchange])
	    else:
		feature = '|'.join([features[0][1], pepchange])
	    
	    if txt_obj:
		return feature, the_txt
	    else:
		return feature

    def is_exon(self, feature):
	"""Detemines if feature is exon"""
	if 'exon' in feature or 'utr' in feature:
	    return True
	else:
	    return False
	
    def cross_exon_bounds(self, pos, direction, txt, exon_num):
	"""Determines if exon boundary is crossed given position and direction.
	"""
	exon = txt.get_exon_bounds(exon_num)
	
	if direction == 'L' and int(pos) >= int(exon[1]):
	    return True
	elif direction == 'R' and int(pos) <= int(exon[0]):
	    return True

	return False
    
    def del_is_intron(self, start, end, txt):
	"""Determines if coordinate is intron"""
	for i in range(len(txt.exons)-1):
	    intron = [int(txt.exons[i][1])+1, int(txt.exons[i+1][0])-1]
	    
	    olap_size = overlap([start, end], intron)
	    	    
	    if float(olap_size) / float(cardinality([start, end])) > 0.9 and\
	       float(olap_size) / float(cardinality(intron)) > 0.9:
		return True
	    
	return False
    
    def localize_feature(self, coord, txt):
        """Identifies feature in given transcript of a given coordinate.
	Features include intron, exon, splice-accepotor, splice-donor, UTR
	"""   
        feature = None	
	# search introns first
	for i in range(len(txt.exons)-1):
	    intron = [int(txt.exons[i][1])+1, int(txt.exons[i+1][0])-1]
	    if subsume([coord, coord], intron):
		if txt.strand == '+':
		    feature = 'intron' + str(i+1)
		else:
		    feature = 'intron' + str(len(txt.exons) - i - 1)
			    
		splice = None
		if subsume([coord, coord], [intron[0], intron[0]+1]):
		    if txt.strand == '+':
			splice = 'splice-donor'
		    else:
			splice = 'splice-acceptor'
		elif subsume([coord, coord], [intron[-1]-1, intron[-1]]):
		    if txt.strand == '+':
			splice = 'splice-acceptor'
		    else:
			splice = 'splice-donor'
			
		if splice:
		    feature += "(%s)" % (splice)  
		break
	    
	# then exons
	if not feature:
	    for i in range(len(txt.exons)):
		if subsume([coord, coord], txt.exons[i]):
		    if txt.strand == '+':
			feature = 'exon' + str(i+1)
		    else:
			feature = 'exon' + str(len(txt.exons)-i)
		    break
		    
	# must be utr if not in intron and not in exons inside cds	    
	if not feature:
	    if (int(coord) < int(txt.cdsStart) and txt.strand == '+') or (int(coord) > int(txt.cdsEnd) and txt.strand == '-'):
		feature = '5utr'
	    else:
		feature = '3utr'
       
        return feature

    def identify_feature(self, start, end, txt, exact):
	"""Finds feature on transcript given start, end coordinates"""
        feature = None

	# identifies feature of start and end coordinate
	# merge them if they are the same
        if not exact:
            start_feature = self.localize_feature(start, txt)
            end_feature = self.localize_feature(end, txt)

            if start_feature == end_feature:
                feature = start_feature
            elif start_feature and end_feature:
                feature = ';'.join((start_feature, end_feature))
            elif start_feature:
                feature = start_feature
            elif end_feature:
                feature = end_feature
            
	# return exon or intron if start or end coordinates land of boundaries
        else:  
	    # exon
            for i in range(len(txt.exons)):
                if int(txt.exons[i][0]) == int(start) and int(txt.exons[i][1]) == int(end):
                    if txt.strand == '+':
                        feature = 'exon' + str(i+1)
                    else:
                        feature = 'exon' + str(len(txt.exons)-i)
                    break           
            if feature:
                return feature

            #intron:
            for i in range(len(txt.exons)-1):
                if int(txt.exons[i][1])+1 == int(start) and int(txt.exons[i+1][0])-1 == int(end):
                    if txt.strand == '+':
                        feature = 'intron' + str(i+1)
                    else:
                        feature = 'intron' + str(len(txt.exons)-i)
                    break     
            return feature

        return feature

    def find_pep_change(self, coord, txt, refseq, variant, change):
	"""Finds peptide change due to variant"""
        cdna_original = self.construct_cdna(coord, txt, refseq)
        cdna_changed = self.construct_cdna(coord, txt, refseq, variant, change)

        if txt.strand == '+':
            pep_original = translate(cdna_original, orient='+', frame=0)
            pep_changed = translate(cdna_changed, orient='+', frame=0)
        else:
            pep_original = translate(cdna_original, orient='-', frame=0)
            pep_changed = translate(cdna_changed, orient='-', frame=0)

        if not pep_changed or not pep_original:
            return 'na'
        
        mutation = pep_change(pep_original, pep_changed)

        return mutation
        
    def construct_cdna(self, coord, txt, refseq, variant=None, change=None):
	"""Constructs transcript sequence given variant.
	Variants: SNV, insertion, duplication, deletion
	"""
        cdna = ""
        remove = False     
        for i in range(len(txt.exons)):
            if not overlap(txt.exons[i], [txt.cdsStart, txt.cdsEnd]):
                continue
                
	    # extracts reference exon sequence
            if subsume([txt.cdsStart, txt.cdsStart], txt.exons[i]):
                start = int(txt.cdsStart) + 1
            else:
                start = txt.exons[i][0]
            if subsume([txt.cdsEnd, txt.cdsEnd], txt.exons[i]):
                end = txt.cdsEnd
            else:
                end = txt.exons[i][1]
	    exon = refseq.GetSequence(coord[0], int(start), int(end))

	    # modifies exon sequence based on variant
            if change:
                if change.lower() == 'snv' and subsume(coord[1:], [int(start), int(end)]):
                    bases_changed = int(coord[1])-int(start), int(coord[2])-int(start)
                    before_change = exon[:bases_changed[0]]
                    after_change = exon[bases_changed[1]+1:]
                    exon = before_change + variant + after_change
                        
                elif change.lower() in ('ins', 'dup', 'ITD', 'PTD') and subsume(coord[1:], [int(start), int(end)]):
                    base_to_insert = int(coord[1])-int(start)+1
                    exon = exon[:base_to_insert] + variant + exon[base_to_insert:]

                elif change.lower() == 'del':
                    if subsume(coord[1:], [int(start), int(end)]):
                        bases_deleted = int(coord[1])-int(start), int(coord[2])-int(start)
                        new_exon = exon[:bases_deleted[0]] + exon[bases_deleted[1]+1:]
                        exon = new_exon
                        
                    elif subsume([coord[1], coord[1]], [int(start), int(end)]):
                        first_base_deleted = int(coord[1])-int(start)
                        exon = exon[:first_base_deleted]
                        remove = True
                        
                    elif subsume([coord[2], coord[2]], [int(start), int(end)]):
                        last_base_deleted = int(coord[2])-int(start)
                        exon = exon[last_base_deleted+1:]
                        remove = False

                    elif i >0 and subsume([coord[1], coord[1]], [int(txt.exons[i-1][1])+1, int(txt.exons[i][0]-1)]):
                        if not remove:
                            remove = True
                            exon = ''
                        else:
                            remove = False

                    elif i >0 and subsume([coord[2], coord[2]], [int(txt.exons[i-1][1])+1, int(txt.exons[i][0]-1)]):
                        if not remove:
                            remove = True
                            exon = ''
                        else:
                            remove = False
                        
                    elif remove:
                        exon = ''
                        
            cdna += exon
                   
        return cdna
    
