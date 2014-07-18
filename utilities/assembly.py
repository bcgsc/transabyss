"""
This module provides methods for extracting contig sequences from assembly file

Author: Readman Chiu rchiu@bcgsc.ca
"""
import os
import sys
import re
from utilities import cfg

script_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    
def get_contig_num(num_orient):
    number_pattern = re.compile('(\d+)')
    return number_pattern.search(num_orient).group(1)

class Assembly:
    """Stores assembly info and mainly functions to extract contig sequences"""
    def __init__(self, lib, abyss=None, k=None, loc=None, fasta=None):
        self.lib = lib
        self.k = k
        self.loc = loc
        self.fasta = fasta
        self.junctions = []
        self.pet = None
        self.abyss = abyss
        
        if self.loc:
            self.pet = loc + "/" + lib + "-contigs.fa"
            self.overlap = loc + "/" + lib + "-3-overlap.fa"
            self.adj_bubbles = loc + "/" + lib + "-4.adj"
            self.adj_islands = loc + "/" + lib + "-4.adj"
            self.indel = loc + "/" + lib + "-indel.fa"

    def get_contigs(self, ids=[], length=0, len_range=[], sequence=False):
        """Returns contig objects in assembly"""
        contigs = []

        fasta = None
        if self.fasta:
            fasta = self.fasta
        elif self.pet:
            fasta = self.pet

        # convert contig ids from list to dict if list given
        ids_dict = None
        if ids:
            ids_dict = dict((c, True) for c in ids)
            
        if fasta and os.path.exists(fasta):
            seq = ""
            last_contig = None
            for line in open(fasta, 'r'):
                if line[0] == '>':
                    # store last contig
                    if last_contig:
                        if sequence and len(seq) > 0:
                            last_contig.sequence = seq
                            if not last_contig.length:
                                last_contig.length = len(seq)                              
                        contigs.append(last_contig)
                            
                    seq = ""        
                    cols = line.rstrip('\n').lstrip('>').split(' ')

                    good = True
                    if ids and not ids_dict.has_key(cols[0]):
                        good = False
                        
                    if good:
                        if length and int(length) > 0 and int(cols[1]) != int(length):
                            good = False

                        if len_range:
                            #proper range given
                            if int(len_range[1]) > int(len_range[0]):
                                if (int(cols[1]) < int(len_range[0])) or (int(cols[1]) > int(len_range[1])):
                                    good = False
                            #first number less than second number, indicates check minimum
                            elif int(cols[1]) < int(len_range[0]):
                                good = False
                    if good:
                        size, coverage, children = None, None, None
                        if len(cols) >= 3:
                            if cols[1] and not re.search('\D', cols[1]):
                                size = int(cols[1])
                            if cols[2] and not re.search('\D', cols[2]):
                                coverage = int(cols[2])
                        if len(cols) > 3 and cols[3]:
                            children = cols[3]
                        last_contig = Contig(cols[0], length=size, coverage=coverage, children=children)
                    else:
                        last_contig = None
                    
                elif sequence and line[0] != ' ':
                    seq += line.rstrip('\n')

            if last_contig:
                if sequence and len(seq) > 0:
                    last_contig.sequence = seq
                    if not last_contig.length:
                        last_contig.length = len(seq)
                        
                contigs.append(last_contig)

            seq = ""
            
        return contigs

    def contigs2dict(self, contigs):
        """Converts list of contig objects to dictionary"""
        contig_dict = {}
        for contig in contigs:
            contig_dict[contig.num] = {}
            contig_dict[contig.num]['sequence'] = contig.sequence
            contig_dict[contig.num]['coverage'] = contig.coverage
            contig_dict[contig.num]['length'] = contig.length

        return contig_dict

    def output(self, out_file):
        """Outputs contigs in fasta"""
        out = open(out_file, 'w')

        for contig in self.contigs:
            out.write(contig.fasta())

        out.close()

    def get_seqs(self, contigs):
        """Extracts contig sequences given list of contig objects"""
        ids = dict((c.num, c) for c in contigs)

        fasta = None
        if self.pet:
            fasta = self.pet
        elif self.fasta:
            fasta = self.fasta

        if fasta and os.path.exists(fasta):
            seq = ""
            contig = None
            for line in open(fasta, 'r'):
                if line[0] == '>':
                    if contig and len(seq) > 0:
                        contig.sequence = seq
                        if not contig.length:
                            contig.length = len(seq)
                            
                    seq = ""
                    cols = line.rstrip('\n').lstrip('>').split(' ')

                    if ids.has_key(cols[0]):
                        contig = ids[cols[0]]
                    else:
                        contig = None
                        
                elif line[0] != ' ':
                    seq += line.rstrip('\n')
                    
            if contig and len(seq) > 0:
                contig.sequence = seq
                if not contig.length:
                    contig.length = len(seq)
    
    def get_seq(self, contig):
        """Use grep to extract contig sequence"""
        process = subprocess.Popen(["grep", "-1", contig, self.fasta], shell=False, stdout=subprocess.PIPE)

        sequence = process.communicate()[0].rstrip("\n").split("\n")[-1]
        #print contig, sequence
        if len(sequence) > 0:
            return sequence
        else:
            return None
                      
class Contig:
    """Stores and provides methods for accessing individual contig info"""
    def __init__(self, num, set=False, pet=False, length=0, coverage=0, k=None, children=None):
        self.num = num
        self.length = length
        self.coverage = coverage
        self.k = k
        self.children = children
        self.sequence = None

    def is_pet(self):
        """Determine if contig is PET"""
        if self.children:
            return True
        else:
            return False

    def fasta(self):
        """Returns sequence in FASTA format"""
        headers = [str(self.num)]
        if self.length:
            headers.append(str(self.length))
        if self.coverage:
            headers.append(str(self.coverage))
        if self.children:
            headers.append(str(self.children))
            
        line = ">" + ' '.join(headers) + '\n'
        
        if self.sequence:
            seq = self.sequence
        else:
            seq = ""
        line += seq + "\n"

        return line

    def normalized_coverage(self):
        """Returns normalized coverage of contig"""
        if self.coverage and int(self.coverage) >= 0 and self.length and int(self.length) > 1 and self.k:
            return "%.1f" % (float(self.coverage) / (float(self.length) - float(self.k) + 1.0))
        else:
            return None



