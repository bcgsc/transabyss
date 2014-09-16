"""
This module provides methods for dealing with individual splicing events.

Author: Readman Chiu rchiu@bcgsc.ca
"""
import re
from utilities.intspan import subsume, overlap, subtract, cardinality
from utilities.tools import reverse_complement, set_attrs
from utilities.bam import BAM
from pep_change import translate, pep_change
from transcript import Transcript

class Event:
    """Individual novel splicing events identified by model_matcher.py"""
    splice_motifs = False
    
    headers = [
        'type',
        'contig',
        'transcript',
        'gene',
        'exons',
        'align_blocks',
        'genome_coord',
        'contig_coord',
        'splice',
        'multi_3',
        'size',
        'orf'
    ]

    headers_support = [
        'spanning_reads',
        'coverage',
    ]

    def __init__(self, type):
        self.event_type = type
        self.contig = None
        self.transcript = None
        self.model = None
        self.exons = []
        self.chrom = None
        self.exon_coords = []
        self.align_blocks = []
        self.align_coords = []
        self.novel = True
        self.artefact = False
        self.align = None
        self.txt = None
        self.filter_result = None
        self.orf = None
        self.genome_coord = None
        self.contig_coord = None

    def details(self, refseq=None):
        """Sets relevant attributes of events"""
        if self.event_type == "retained_intron":
            apart = int(self.exon_coords[1][0]) - 1 - int(self.exon_coords[0][1])
            self.multi_3 = self.multiple_3(apart)
        else:
            self.multi_3 = None

        # splice sites
        self.splice = None
        if self.event_type in ("AS5", "AS3"):
            ss = self.get_splice_sites(self.align_blocks[0], left=self.left, right=self.right)
            if ss and len(ss) > 0:
                self.splice = ss[0]
        elif self.event_type in ("novel_exon", "AS53", "novel_transcript", "novel_utr"):
            splice_sites = []
            for block in self.align_blocks:
                ss = ",".join(self.get_splice_sites(block, left=True, right=True))
                if ss:
                    splice_sites.append(ss)
            if splice_sites:
                self.splice = "-".join(splice_sites)        
        elif self.event_type == "novel_intron":
            ss = self.get_splice_sites(self.align_blocks[0], left=False, right=True, whole=True)
            if ss and len(ss) > 0:
                self.splice = ss[0]
            
        # size
        self.size = None
        if self.event_type in ("novel_exon", "AS53", "novel_transcript", "novel_utr", "novel_intron"):
            sizes = []
            for block in self.align_coords:
                sizes.append(str(cardinality(block)))

            self.size = ",".join(sizes)

        # orf
        if refseq:
            self.find_pep_change(refseq)
            
    @classmethod
    def parse(cls, record):
        """Parses tabulated output into object"""
        cols = record.rstrip('\n').split('\t')
        data = {}
        headers = cls.headers[:]
        headers.extend(cls.headers_support)
        for i in range(len(headers)):
            if cols[i+1] != 'na' and (headers[i] == 'spanning_reads' or headers[i] == 'coverage'):
                if ',' in cols[i+1]:
                    data[headers[i]] = [int(n) for n in cols[i+1].split(',')]
                else:
                    data[headers[i]] = int(cols[i+1])
            else:
                data[headers[i]] = cols[i+1]
            
        e = Event(cols[1])
        set_attrs(e, data)
        return e

    def output(self, gff=False, parsed=False):
        """Outputs event in tabulated format"""
        data = {}
        # initialize
        for header in self.headers:
            data[header] = 'na'
        for header in self.headers_support:
            data[header] = 'na'
        
        # direct translation from object into output, no need to process
        if parsed:
            headers = self.headers[:]
            headers.extend(self.headers_support)
            for i in range(len(headers)):
                value = getattr(self, headers[i])
                if type(value) == list:
                    data[headers[i]] = ','.join([str(n) for n in value])
                else:
                    data[headers[i]] = value
                    
        else:
            data['type'] = self.event_type
            data['contig'] = self.contig
            data['align_blocks'] = ",".join([str(b) for b in self.align_blocks])
            data['genome_coord'] = self.coordinate()
            data['contig_coord'] = self.contig_coord
            if self.contig_coord is None:
                data['contig_coord'] = self.get_contig_coord()
            
            if self.event_type != 'novel_transcript':
                if self.event_type == 'read-through' and self.txt2:
                    data['transcript'] = '%s;%s' % (self.txt.name, self.txt2.name)
                else:
                    data['transcript'] = self.txt.name 
                    
                if self.event_type == 'read-through' and self.txt2:
                    genes = []
                    for t in (self.txt, self.txt2):
                        if t.alias:
                            genes.append(t.alias)
                        else:
                            genes.append('na')
                            
                    data['gene'] = ';'.join(genes)
                else:
                    if self.txt.alias:
                        data['gene'] = self.txt.alias
                    else:
                        data['gene'] = 'na'
                
                exons = []                
                for e in self.exons:
                    exons.append(self.convert_gene_exon_number(self.txt, int(e)))
                data['exons'] = ",".join([str(e) for e in exons])

            for item in ('size', 'multi_3', 'splice', 'orf'):
                try:
                    value = getattr(self, item)
                except AttributeError:
                    pass
                else:
                    if value:
                        data[item] = value

            for item in Event.headers_support:
                try:
                    value = getattr(self, item)                    
                except AttributeError:
                    data[item] = 'na'
                else:
                    if value is not None:                            
                        if value != 'na' and type(value).__name__ == 'list':
                            value = [str(v) for v in value]
                            value = ','.join(value)
                                                        
                        data[item] = str(value)
                    else:
                        data[item] = 'na'
                        
        fields = []
        for item in Event.headers:
            fields.append(data[item])
        for item in Event.headers_support:
            fields.append(data[item])
            
        return ('\t'.join(str(f) for f in fields))
        
    def convert_gene_exon_number(self, txt, positive_exon_num):
        """Converts exon number to one taking gene strand into consideration
        'positive_exon_num': assumed 1-based
        """
        exon_num = positive_exon_num
        if txt.strand == '-':
            exon_num = len(txt.exons) - positive_exon_num + 1
            
        return exon_num
    
    def get_splice_sites(self, block, left=False, right=False, whole=False):
        """Reports splice site sequence and motif
        'left'=donor, 'right'=acceptor, 'whole'=donor + acceptor
        """
        ss = []
        if self.align.splice_sites:
            splice_sites = self.align.splice_sites

            orient = None
            if self.txt:
                if self.txt.strand == '-':
                    orient = '-'
                elif self.txt.strand == '+':
                    orient = '+'
            elif self.align.orient:
                orient = self.align.orient
            
            motif = '?'                
            if left and int(block) > 1:
                splice_site = splice_sites[block-2]
                
                if orient == '+' and Event.splice_motifs.has_key(splice_site):
                    motif = Event.splice_motifs[splice_site]
                elif orient == '-' and Event.splice_motifs.has_key(reverse_complement(splice_site)):
                    motif = Event.splice_motifs[reverse_complement(splice_site)]
                    
                if not whole:
                    ss.append('%s%s(%s)' % (splice_site[:2].lower(), splice_site[-2:].upper(), motif))
                else:
                    ss.append('%s%s(%s)' % (splice_site[:2].upper(), splice_site[-2:].upper(), motif))
                    
            if right and int(block) < len(self.align.blocks):
                splice_site = splice_sites[block-1]
                
                if orient == '+' and Event.splice_motifs.has_key(splice_site):
                    motif = Event.splice_motifs[splice_site]
                elif orient == '-' and Event.splice_motifs.has_key(reverse_complement(splice_site)):
                    motif = Event.splice_motifs[reverse_complement(splice_site)]
                  
                if not whole:
                    ss.append('%s%s(%s)' % (splice_site[:2].upper(), splice_site[-2:].lower(), motif))
                else:
                    ss.append('%s%s(%s)' % (splice_site[:2].upper(), splice_site[-2:].upper(), motif))

        else:
            print 'cannot extract splice sites: %s' % (self.align.query)
                        
        return ss
            
    def coordinate(self):
        """Returns UCSC format of event coordinate"""
        coord = ""
        if self.event_type == "skipped_exon" or self.event_type == "retained_intron":
            coord = self.align.target + ":" + str(self.exon_coords[0][0]) + "-" + str(self.exon_coords[-1][1])
        elif self.event_type == "novel_intron":
            coord = self.align.target + ":" + str(int(self.align_coords[0][1])+1) + "-" + str(int(self.align_coords[-1][0])-1)
        elif self.event_type == "novel_exon" or self.event_type == "novel_utr" or self.event_type == "novel_transcript" or self.event_type == "AS53" or self.event_type == 'read-through':
            coord = self.align.target + ":" + str(self.align_coords[0][0]) + "-" + str(self.align_coords[-1][1])
        elif self.event_type == "AS5" or self.event_type == "AS3":
            coord = self.align.target + ":" + str(self.align_coords[0]) + "-" + str(self.align_coords[-1])

        return coord
    
    def get_contig_coord(self):
        """Returns contig coordinate of event"""
        coords = []
        if self.align.query_blocks and self.align_blocks:
            coords.append(int(self.align.query_blocks[self.align_blocks[0] - 1][0]))
            coords.append(int(self.align.query_blocks[self.align_blocks[-1] - 1][1]))  
            coords.sort(key = int)            
            return str(coords[0]) + '-' + str(coords[1])
        else:
            return '-' 

    def multiple_3(self, apart):
        """Determine if size of event is a multiple of 3"""
        if apart % 3 == 0:
            return 'True'
        else:
            return 'False'

    def set_novelty(self, txts, matches=None):
        """Determines if event is novel"""
        novel_events = []
        
        if self.event_type == "novel_exon" or self.event_type == "AS53" or self.event_type == "novel_utr":
            for txt in txts:
                blocks_to_delete = []

                for b in range(len(self.align_blocks)):
                    for e in range(len(txt.exons)):
                        novel = True
                        if subsume(self.align_coords[b], txt.exons[e]):
                            novel = False
                        #novel utr - requires just one edge to align
                        if novel and self.event_type == 'novel_utr':
                            if int(self.align_coords[b][0]) == int(txt.exons[e][0]) or int(self.align_coords[b][1]) == int(txt.exons[e][1]):
                                novel = False
                        if not novel:
                            blocks_to_delete.append(b)
                            break

                if blocks_to_delete:
                    for b in sorted(blocks_to_delete, reverse=True):
                        del self.align_blocks[b]
                        del self.align_coords[b]

            if not self.align_blocks:
                self.novel = False
                
        elif self.event_type == 'read-through':
            if len(self.align_coords) == 1:
                start, end = self.align.blocks[self.align_blocks[0] - 2][1], self.align.blocks[self.align_blocks[0] - 1][0]                
                # see if any single transcript contains the exon junction
                for txt in txts:
                    found_start, found_end = None, None                    
                    for i in range(len(txt.exons) - 1):
                        if int(txt.exons[i][1]) == start and int(txt.exons[i + 1][0]) == end:
                            found_start, found_end = i, i + 1
                            self.novel = False
                            break
                                                              
                    if not self.novel:
                        break
                
        elif self.event_type == "retained_intron":
            multi = False
            if int(self.exons[-1]) - int(self.exons[0]) > 1:
                multi = True
                self.novel = False
                
            for i in range(len(self.exon_coords)-1):
                retained_intron = [int(self.exon_coords[i][1])+1, int(self.exon_coords[i+1][0])-1]
                middle_exons = {}
                for txt in txts:
                    exons_txt = []
                    for j in range(len(txt.exons)):
                        exon = txt.exons[j]
                        #terminal exon, require subsume
                        if j == 0 or j == len(txt.exons)-1:
                            if subsume(retained_intron, exon):
                                exons_txt.append(exon) 
                        #middle exons, require just overlap
                        elif overlap(exon, (retained_intron[0], retained_intron[1])):
                            exons_txt.append(exon)

                    if exons_txt:
                        middle_exons[txt] = exons_txt

                # only time when original event is novel WITHOUT testing is when it's a single ri and it is clear of overlapping exons
                if len(middle_exons.keys()) == 0 and not multi:
                    self.novel = True
                else:
                    self.novel = False

                    # substract overlapping exons
                    if middle_exons.values() and middle_exons.values()[0]:
                        true_retained_intron = subtract(retained_intron, middle_exons.values())
                    else:
                        true_retained_intron = [retained_intron]

                    # if there is still some intron left after subtraction
                    if true_retained_intron:
                        # create new events
                        for ri in true_retained_intron:
                            event = {'contig': self.contig, 
                                     'chrom':self.chrom, 
                                     'align_blocks':self.align_blocks, 
                                     'align_coords':self.align_coords, 
                                     'type':self.event_type, 
                                     'novel':True
                                     }
                            # determine flanking exons, and transcript by frequency of flanking coordinates
                            flanks = {}
                            for txt in txts:
                                for i in range(len(txt.exons)-1):
                                    left = txt.exons[i]
                                    right = txt.exons[i+1]
                                    if left[1]+1 == ri[0] and right[0]-1 == ri[1]:
                                        #print "ri flanks", txt.full_name(), left, right, ri
                                        coord = ",".join((str(left[0]), str(left[1]), str(right[0]), str(right[1])))
                                        if not flanks.has_key(coord):
                                            flanks[coord] = [[txt,i,i+1]]
                                        else:
                                            flanks[coord].append([txt,i,i+1])

                            if len(flanks.keys()) > 0:
                                # use the most commom frequent exons
                                flanks_sorted = flanks.keys()
                                flanks_sorted.sort(lambda x,y: len(flanks[y])-len(flanks[x]))

                                # use the orginally assigned transcript if possible
                                same_txt = False
                                exon_coords = flanks_sorted[0].split(',')
                                exons = []
                                for txt,e1,e2 in flanks[flanks_sorted[0]]:
                                    if txt.full_name() == self.transcript:
                                        same_txt = True
                                        event['transcript'] = txt.full_name()
                                        event['exons'] = [e1+1,e2+1]
                                        event['exon_coords'] = [exon_coords[:2], exon_coords[2:]]
                                        event['txt'] = txt
                                        break

                                if not same_txt:
                                    txt,e1,e2 = flanks[flanks_sorted[0]][0]
                                    event['transcript'] = txt.full_name()
                                    event['exons'] = [e1+1,e2+1]
                                    event['exon_coords'] = [exon_coords[:2], exon_coords[2:]]
                                    event['txt'] = txt

                                novel_events.append(event)
            
        elif self.event_type == "skipped_exon":
            for txt in txts:
                for e in range(len(txt.exons)-1):
                    intron_span = [int(txt.exons[e][1])+1, int(txt.exons[e+1][0])-1]
                    if self.align_coords and self.exon_coords:
                        if subsume([int(self.exon_coords[0][0]), int(self.exon_coords[-1][1])], intron_span):
                            self.novel = False
                            break
                        
                if not self.novel:
                    break

        elif self.event_type == "novel_intron":
            novel_intron_span = [int(self.align_coords[0][1])+1, int(self.align_coords[1][0])-1]
            novel_intron_size = int(self.align_coords[1][0]) - 1 - int(self.align_coords[0][1])
            for txt in txts:
                for e in range(len(txt.exons)-1):
                    intron_span = [int(txt.exons[e][1])+1, int(txt.exons[e+1][0])-1]

                    if novel_intron_span[0] == intron_span[0] and novel_intron_span[1] == intron_span[1]:
                        self.novel = False
                        break

                if not self.novel:
                    break

        elif 'AS' in self.event_type and self.edge == 'left':
            for txt in txts:
                for exon in txt.exons:
                    if int(self.align_coords[0]) == int(exon[0]):
                        self.novel = False
                        break

                if not self.novel:
                    break

        elif 'AS' in self.event_type and self.edge == 'right':
            for txt in txts:
                for exon in txt.exons:
                    if int(self.align_coords[1]) == int(exon[1]):
                        self.novel = False
                        break

                if not self.novel:
                    break

        return novel_events
    
    def last_matched(self):
        """Determines last-matched alignment block and exon
        This method is only used for is_read_through()
        """
        last_matched_block = self.align.blocks[self.last_matched_block - 1]
        if self.exons[0] == 1:
            # the last-matched exon is supposed to be the second last flanking exon
            if len(self.txt.exons) > 1 and overlap(last_matched_block, self.txt.exons[self.exons[0] - 1 + 1]):
                last_matched_exon = self.txt.exons[self.exons[0] - 1 + 1]                
            # but sometimes it will be the last flanking exon if the unmatched block doesn't
            # overlap the second last exon
            else:
                last_matched_exon = self.txt.exons[self.exons[0] - 1]            
        else:
            if len(self.txt.exons) > 1 and overlap(last_matched_block, self.txt.exons[self.exons[0] - 1 - 1]):
                last_matched_exon = self.txt.exons[self.exons[0] - 1 - 1]        
            else:
                last_matched_exon = self.txt.exons[self.exons[0] - 1]
                                
        return last_matched_block, last_matched_exon
    
    def is_read_through(self, txts, mm):
        """Determines if event is read-through"""
        last_matched_block, last_matched_exon = self.last_matched()        
        for txt2 in txts:
            if txt2.strand != self.txt.strand:
                continue
            
            if txt2.model != self.txt.model:
                continue
            
            if txt2.name == self.txt.name or txt2.alias == self.txt.alias:
                continue
                        
            if not overlap([self.align_coords[0][0], self.align_coords[-1][1]], [txt2.txStart, txt2.txEnd]) or\
               overlap([self.txt.txStart, self.txt.txEnd], [txt2.txStart, txt2.txEnd]):
                continue
                            
            if overlap(last_matched_block, [txt2.txStart, txt2.txEnd]):
                continue
                                                        
            result = mm.match_exons(self.contig, txt2.full_name(), self.align_coords, txt2.exons, txt2.chrom, strand=txt2.strand)                        
            if result and len(result.matched_blocks) == len(self.align_blocks):                                
                exon_bounds_matched = True
                for i in range(len(result.matched_blocks)):                                    
                    # only 1 boundary has to be flush if it's terminal block
                    if i == len(self.align_blocks) - 1:
                        if self.txt.txStart < txt2.txStart:
                            if self.align_coords[result.matched_blocks[i] - 1][0] != txt2.exons[result.matched_exons[i] - 1][0]:
                                exon_bounds_matched = False
                                
                        else:
                            if self.align_coords[result.matched_blocks[i] - 1][1] != txt2.exons[result.matched_exons[i] - 1][1]:
                                exon_bounds_matched = False
                        
                    # both boundaries have to be flush if it's not terminal block
                    else:
                        if not(self.align_coords[result.matched_blocks[i] - 1][0] == txt2.exons[result.matched_exons[i] - 1][0] and\
                               self.align_coords[result.matched_blocks[i] - 1][1] == txt2.exons[result.matched_exons[i] - 1][1]):
                            exon_bounds_matched = False
                            
                if not exon_bounds_matched:
                    continue
                
                if self.txt.txStart < txt2.txStart:
                    txt_span = [int(self.txt.txEnd) + 1, int(txt2.txStart) - 1]
                else:
                    txt_span = [int(txt2.txEnd) + 1, int(self.txt.txStart) - 1]
                                
                # make sure there is no transcripts in between the 1st and 2nd transcripts
                has_txt_between = False
                for t in txts:
                    if t.name == self.txt.name or t.name == txt2.name:
                        continue
                                
                    if subsume([t.txStart, t.txEnd], txt_span):
                        has_txt_between = True
                        break
                        
                    if not has_txt_between:                        
                        if self.txt.alias and txt2.alias and type(self.txt.alias) is str and type(txt2.alias) is str:
                            if not Transcript.same_family(self.txt.alias, txt2.alias):
                                self.event_type = 'read-through'
                                self.txt2 = txt2
                                
    def longest_orf(self):
        """Returns longest ORF of all 6-frame translation"""
        orf = translate(self.align.contig.sequence, full=True)
        if orf:
            coverage = "%.2f" % ((float(orf[1]) - float(orf[0]) + 1) / float(len(self.align.contig.sequence)))
            self.orf = "%d-%d;%s;%s..%s;%s" % (orf[0], orf[1], coverage, orf[-1][:3], orf[-1][-3:], orf[2])
        else:
            self.orf = 'na'

        return None

    def find_pep_change(self, refseq):
        """Finds effect on protein sequence given event"""
        if self.event_type == 'novel_transcript':
            return self.longest_orf()
        
        coord = re.split('[:-]', self.coordinate())
        variant = None
        if self.event_type == 'novel_utr':
            # include the block that's 'matching' (partially)
            if min(self.align_blocks) == 1:
                b1 = self.align_blocks[0]
                b2 = self.align_blocks[-1] + 1
            else:
                b1 = self.align_blocks[0] - 1
                b2 = self.align_blocks[-1]
            qcoord1 = self.align.query_blocks[b1-1][0]
            qcoord2 = self.align.query_blocks[b2-1][1]
        else:
            qcoord1 = self.align.query_blocks[self.align_blocks[0]-1][0]
            qcoord2 = self.align.query_blocks[self.align_blocks[-1]-1][1]
        if qcoord1 < qcoord2:
            variant = self.align.contig.sequence[qcoord1-1:qcoord2]
        else:
            variant = self.align.contig.sequence[qcoord2-1:qcoord1]
            variant = reverse_complement(variant)
              
        # constructs cDNA sequence of both reference and sequence with event
        cdna_original = self.construct_cdna(coord, self.txt, refseq)
        cdna_changed = self.construct_cdna(coord, self.txt, refseq, variant=variant, change=self.event_type, exons=self.exons)

        frame = 0
        if self.event_type == 'novel_utr':
            frame = None

        if self.txt.strand == '+':
            pep_original = translate(cdna_original, orient='+', frame=0)
            pep_changed = translate(cdna_changed, orient='+', frame=frame)
        else:
            pep_original = translate(cdna_original, orient='-', frame=0)
            pep_changed = translate(cdna_changed, orient='-', frame=frame)

        if not pep_changed or not pep_original or self.event_type == 'read-through':
            return 'na'
        
        self.orf = pep_change(pep_original, pep_changed)

    def construct_cdna(self, coord, txt, refseq, variant=None, change=None, exons=None):
        """Constructs cDNA seqeunce given event"""
        cdna = ""
        for i in range(len(txt.exons)):
            if txt.coding_type() != 'CODING' or not overlap(txt.exons[i], [txt.cdsStart, txt.cdsEnd]):
                continue
                
            if subsume([txt.cdsStart, txt.cdsStart], txt.exons[i]):
                start = int(txt.cdsStart) + 1
            else:
                start = txt.exons[i][0]

            if subsume([txt.cdsEnd, txt.cdsEnd], txt.exons[i]):
                end = int(txt.cdsEnd)
            else:
                end = int(txt.exons[i][1])

            exon = refseq.GetSequence(coord[0], int(start), int(end))
            if change:
                if change.lower() == 'retained_intron' and i+1 == exons[0]:
                    intron = refseq.GetSequence(coord[0], end+1, txt.exons[i+1][0]-1)
                    exon += intron
                    
                elif change.lower() == 'novel_exon' and i+1 == exons[0]:
                    exon += variant

                elif change.lower() == 'skipped_exon' and i+1 in exons:
                    exon = ''

                elif change.lower() == 'novel_intron' and i+1 == exons[0]:
                    bases_deleted = int(coord[1])-int(start), int(coord[2])-int(start)
                    new_exon = exon[:bases_deleted[0]] + exon[bases_deleted[1]+1:]
                    exon = new_exon

                elif change.lower() in ['as5', 'as3', 'as53'] and i+1 == exons[0]:
                    new_exon = refseq.GetSequence(coord[0], int(coord[1]), int(coord[2]))
                    exon = new_exon

                elif change.lower() == 'novel_utr' and i+1 == exons[0]:
                    exon = variant
            
            cdna += exon

        return cdna
