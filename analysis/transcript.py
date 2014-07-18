"""
This module stores transcript annotation as object and provides methods
for extracting and manipulating such information

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
import re
import annotations
import ConfigParser
from utilities.overlap_coord import OverlapCoord
from utilities import tools

PACKAGE_DIR = "/".join(os.path.abspath(__file__).split("/")[:-2])

class Transcript:
    """Stores individual transcript info and provides methods relevant to transcripts"""
    def __init__(self, name):
        self.name = name
	self.chrom = None
	self.strand = None
	self.txStart = None
	self.txEnd = None
	self.cdsStart = None
	self.cdsEnd = None
	self.exonCount = None
	self.exons = None
	self.gene = None
	self.alias = None
	self.model = None
	self.weight = None
	self.length = 0

    def full_name(self):
	"""Returns <transcript>(<gene>)"""
        alias = "NA"
        if self.alias:
            alias = self.alias   
        return "%s(%s)" % (self.name, alias)

    def cds_length(self):
	"""Returns CDS length"""
        if self.cdsStart and self.cdsEnd and self.cdsStart.isdigit() and self.cdsEnd.isdigit():
            return int(self.cdsEnd) - int(self.cdsStart)

        return 0
    
    def txt_length(self):
	"""Returns transcript length"""
	if self.txStart and self.txEnd and self.txStart.isdigit() and self.txEnd.isdigit():
            return int(self.txEnd) - int(self.txStart)

        return 0
    
    def coding_type(self):
	"""Returns transcript type: CODING/NONCODING/NA
	CODING when cdsStart != cdsEnd
	"""
        if self.cdsStart == None or self.cdsEnd == None:
            return 'NA'
	elif str(self.cdsStart).isdigit() and str(self.cdsEnd).isdigit() and self.cdsStart != self.cdsEnd:
            return 'CODING'
        else:
            return 'NONCODING'
	
    def get_exon_bounds(self, exon_num):
	"""Returns exon boundaries given exon number"""
	if int(exon_num) <= len(self.exons):
	    if self.strand == '+':
		return self.exons[int(exon_num) - 1]    
	    else:
		return self.exons[::-1][int(exon_num) - 1]
	    
	return []
    
    def get_intron_bounds(self, intron_num):
	"""Returns intron boundaries given intron number"""
	introns = []
	for i in range(1, len(self.exons)):
	    introns.append([int(self.exons[i - 1][1]) + 1, int(self.exons[i][0]) - 1])
	    
	if int(intron_num) <= len(introns):
	    if self.strand == '+':
		return introns[int(intron_num) - 1]
	    
	    else:
		return introns[::-1][int(intron_num) - 1]
	    
	return []
    
    def utr(self, end):
        """Returns UTR coordinate given 5 or 3""" 
        utr3 = utr5 = None
        if self.coding_type() == 'CODING':
            if self.strand == '+':
                utr5 = [int(self.txStart), int(self.cdsStart) - 1]
                utr3 = [int(self.cdsEnd) + 1, int(self.txEnd)]
            else:
                utr3 = [int(self.txStart) + 1, int(self.cdsStart)]
                utr5 = [int(self.cdsEnd) + 1, int(self.txEnd)]
            
            if end == 5:
                return utr5
            elif end == 3:
                return utr3
            else:
                return None
        else:
            return utr3
	
    def get_sequence(self, refseq, fasta=False, chrom=None):
	"""Extracts transcript sequence"""
	sequence = ''
	if chrom is None:
	    chrom = self.chrom
	else:
	    chrom = chrom
		    
	for i in range(len(self.exons)):
	    exon = refseq.GetSequence(chrom, int(self.exons[i][0]), int(self.exons[i][-1]))
	    sequence += exon.upper()
	  
	if self.strand == '-':
	    sequence = tools.reverse_complement(sequence)
	    
	if fasta:
	    sequence = '>%s_%s\n%s' % (self.name, self.alias, sequence)
	
	return sequence
          
    @classmethod
    def find_by_name(cls, model, annot_file, txt_name):
	"""Returns transcript object(s) of given transcript name"""
        txts = []
        if os.path.exists(annot_file):
            for line in open(annot_file, 'r'):
                if txt_name in line:
                    txt = {
                        'e': annotations.ensembl.parse_line,
                        'r': annotations.refGene.parse_line,
                        'k': annotations.knownGene.parse_line,
                        'a': annotations.aceview.parse_line,
                        'x': annotations.ensg.parse_line,
		        't': annotations.ensembl.parse_line,
		        'g': annotations.ensembl.parse_line,
                        }[model](line)
                    if txt.name == txt_name or txt.alias == txt_name:
                        txts.append(txt)

        return txts		
    
    @classmethod
    def prepare_overlap(cls, genome, models):
	"""Prepares transcript overlapping for given genome and models"""
        overlaps = {}
    
        config = ConfigParser.ConfigParser()
        config.read(PACKAGE_DIR + "/configs/model_matcher.cfg")  
        for m in models.split(','):       
            config = ConfigParser.ConfigParser()
            config.read(PACKAGE_DIR + "/configs/model_matcher.cfg")
            annot_file = index = None

            if genome in config.sections():
                for field in config.options(genome):
                    if field.lower() == m.lower():
                        annot_file = PACKAGE_DIR + "/annotations/" + genome + "/" + config.get(genome, field)

                        if annot_file and os.path.exists(annot_file):
                            index_file = annot_file.replace('.txt', '.idx')                
                            txt_overlap = OverlapCoord(annot_file, index_file)
                            txt_overlap.extract_index()                       
			    overlaps[m] = txt_overlap
                        
        return overlaps

    @classmethod
    def find_overlaps(cls, overlaps, models, chrom, start, end, only_coding=False):
	"""Find transcripts overlapping given coordinate"""
        txts = []
        for model in models:
            lines = overlaps[model].overlap(chrom, int(start), int(end))
            for line in lines:
                txt = {
                    'e': annotations.ensembl.parse_line,
                    'r': annotations.refGene.parse_line,
                    'k': annotations.knownGene.parse_line,
                    'a': annotations.aceview.parse_line,
                    'x': annotations.ensg.parse_line,
		    't': annotations.ensembl.parse_line,
		    'g': annotations.ensembl.parse_line,
                    }[model](line)             
                if only_coding and txt.coding_type().upper() != 'CODING':
		    continue      
                txts.append(txt)
                    
        return txts
    
    @classmethod
    def pick_longest(cls, txts, coding=False):
	"""Finds longest transcript of the list given.
	Length is defined as CDS first, transcript length second
	"""
	best = None
	for txt in txts:
	    if coding and not txt.coding_type().upper() == 'CODING':
		continue
	    if best is None or txt.cds_length() > best.cds_length():
		best = txt
	    elif txt.cds_length() == best.cds_length() and txt.txt_length() > best.txt_length():
		best = txt
		
	return best
    
    @classmethod
    def same_family(cls, gene1, gene2):
	"""Checks if 2 genes belong to the same family based on regex"""
	family_patterns = []
	# e.g. KLF1;KLF2
	family_patterns.append(re.compile(r'(.+?)(\d+)$'))
	# e.g. TUBA1B;TUBA1A, DDX19B;DDX19A
	family_patterns.append(re.compile(r'(.+?\d)([A-Z])$'))
	
	same = False
	# e.g. CATSPER2;CATSPER2P1
	if gene1 in gene2 or gene2 in gene1:
	    same = True
	    return same
	
	for family_pattern in family_patterns:
	    family1 = family2 = None
	    m1 = family_pattern.search(gene1)
	    if m1:
		family1 = m1.group(1)
        
	    m2 = family_pattern.search(gene2)
	    if m2:
		family2 = m2.group(1)
            
	    if family1 is not None and family2 is not None and family1 == family2:
		same = True
    
	return same
    
    @classmethod
    def compare_start(cls, txt1, txt2):
	"""For sorting transcripts"""
	chrom1 = txt1.chrom.lstrip('chr')
	chrom2 = txt2.chrom.lstrip('chr')
	
	if chrom1 == chrom2:
	    if chrom1 < chrom2:
		return -1
	    elif chrom1 > chrom2:
		return 1
	    elif int(txt1.txStart) < int(txt2.txStart):
		return -1
	    elif int(txt1.txStart) > int(txt2.txStart):
		return 1
	    elif int(txt1.txEnd) < int(txt2.txEnd):
		return -1
	    elif int(txt1.txEnd) > int(txt2.txEnd):
		return 1
	    else:
		return 0	    
	else:
	    if chrom1.isdigit() and not chrom2.isdigit():
		return -1
	    elif not chrom1.isdigit() and chrom2.isdigit():
		return 1    
	    elif chrom1 < chrom2:
		return -1
	    elif chrom1 > chrom2:
		return 1
	    elif int(txt1.txStart) < int(txt2.txStart):
		return -1
	    elif int(txt1.txStart) > int(txt2.txStart):
		return 1
	    elif int(txt1.txEnd) < int(txt2.txEnd):
		return -1
	    elif int(txt1.txEnd) > int(txt2.txEnd):
		return 1
	    else:
		return 0
		
		