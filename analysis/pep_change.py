"""
This module provides methods for translating cDNA sequences and finding effects
of amino acid change
Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def translate(sequence, min_protein_length=1, orient=None, frame=None, full=False, all=False):
    """Translates cdna sequence into protein"""
    seq = Seq(sequence)
    orfs = []
    seq_len = len(seq)

    if orient == '+':
        strand_and_base = [(+1, seq)]
    elif orient == '-':
        strand_and_base = [(-1, seq.reverse_complement())]
    else:
        strand_and_base = [(+1, seq), (-1, seq.reverse_complement())]
                               
    for strand, nuc in strand_and_base:        
        for fm in range(3):
            if frame != None and fm != frame:
                continue
            
            trans = str(nuc[fm:].translate())
            trans_len = len(trans)

            aa_start = 0
            aa_end = 0
            while aa_start < trans_len :
                aa_end = trans.find("*", aa_start)
                if aa_end == -1 :
                    aa_end = trans_len - 1
                if aa_end - aa_start >= min_protein_length :
                    if strand == 1 :
                        start = fm + aa_start * 3
                        end = min(seq_len - 1, fm + aa_end * 3 + 3 - 1)
                    else :
                        end = seq_len - 1 - fm - aa_start * 3
                        start = end + 1 - (aa_end - aa_start) * 3 - 3

                    orfs.append((start, end, strand, trans[aa_start:aa_end]))
                aa_start = aa_end + 1

    if len(orfs) > 0:
        if not all:
            orfs.sort(lambda x,y: len(y[3])-len(x[3]))
            if not full:
                return orfs[0][-1]
            else:
                return orfs[0]
        else:
            return orfs
    else:
        return None

def pep_change(pep_original, pep_changed):
    """Report effect of amino acid change"""
    mutation = None

    if pep_original == pep_changed:
        mutation = 'synon'
        
    elif pep_original in pep_changed:
        m = re.search(pep_original, pep_changed)
        if m.start() == 0:
            new_aa = pep_changed[m.end():]
        else:
            new_aa = pep_changed[:m.start()]
        mutation = "ins%s" % (new_aa)

    else:
        first_aa_diff = None
        for i in range(len(pep_original)):
            if i < len(pep_changed) and pep_original[i] != pep_changed[i]:
                first_aa_diff = i
                break

        aa_original_after = pep_original[first_aa_diff:]
        aa_changed_after = pep_changed[first_aa_diff:]

        # deletion
        if aa_changed_after in aa_original_after:
            try:
                m = re.search(aa_changed_after, aa_original_after)
            except:
                return None
            
            aa_deleted = []
            for i in range(0, m.start()):
                aa = pep_original[first_aa_diff + i] +  str(first_aa_diff + 1 + i)
                aa_deleted.append(aa)

            if len(aa_deleted) > 1:
                mutation = "%s_%sdel" % (aa_deleted[0], aa_deleted[-1])
            elif not aa_deleted:
                aa_deleted.append(pep_original[len(pep_changed)] +  str(len(pep_changed)+1))
                aa_deleted.append(pep_original[-1] +  str(len(pep_original)))
                mutation = "%s_%sdel" % (aa_deleted[0], aa_deleted[-1])
            else:
                mutation = "%sdel" % (aa_deleted[0])

        # insertion
        elif aa_original_after in aa_changed_after:
            try:
                m = re.search(aa_original_after, aa_changed_after)
            except:
                return None

            if not first_aa_diff:
                first_aa_diff = len(pep_original)
                aa_before_insertion = pep_original[first_aa_diff-1] + str(first_aa_diff)

                mutation = "%sins%s" % (aa_before_insertion, pep_changed[len(pep_original):])
            else:
                aa_before_insertion = pep_original[first_aa_diff-1] + str(first_aa_diff)
                aa_after_insertion = pep_original[first_aa_diff] + str(first_aa_diff+1)
            
                aa_inserted = []
                for i in range(0, m.start()):
                    aa = pep_changed[first_aa_diff + i]
                    aa_inserted.append(aa)

                mutation = "%s_%sins%s" % (aa_before_insertion, aa_after_insertion, ''.join(aa_inserted))

        # substitution
        elif aa_original_after[1:] == aa_changed_after[1:]:
            mutation = "%s%s%s" % (aa_original_after[0], first_aa_diff+1, aa_changed_after[0])

        # indel, frameshift
        else:
            pep_original_rev = pep_original[::-1]
            pep_changed_rev = pep_changed[::-1]
        
            for i in range(len(pep_original_rev)):
                if i < len(pep_changed_rev) and pep_original_rev[i] != pep_changed_rev[i]:
                    first_aa_diff_rev = i
                    break

            if first_aa_diff_rev > 0:
                last_aa_same = -1 * first_aa_diff_rev

                aa_deleted = []
                for i in range(first_aa_diff, len(pep_original) + last_aa_same):
                    aa = pep_original[i] + str(i+1)
                    aa_deleted.append(aa)

                aa_inserted = pep_changed[first_aa_diff:last_aa_same]

                if len(aa_deleted) > 1:
                    mutation = "%s_%sdel" % (aa_deleted[0], aa_deleted[-1])
                else:
                    mutation = "%sdel" % (aa_deleted[0])
                mutation += "ins" + aa_inserted

            # frameshift
            else:
                first_aa_changed = "%s%s%s" % (pep_original[first_aa_diff], str(first_aa_diff+1), pep_changed[first_aa_diff])
                mutation = first_aa_changed + 'fs'
                aa_inserted = pep_changed[first_aa_diff:]
                if '*' in aa_inserted:
                    m = re.search('\*', aa_inserted)
                    mutation += 'X' + str(m.start()+1)

    return mutation
