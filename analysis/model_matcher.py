"""
This module provides captures novel splicing events and determines gene and transcript coverage.

Author: Readman Chiu rchiu@bcgsc.ca
"""
__version__ = '1.5.1'

import os
import sys
import re
import glob
import ConfigParser
import pprint
import linecache
from optparse import OptionParser
import event, event_support
from annotations import ensembl, knownGene, refGene, aceview, ensg
from utilities import track, assembly
from utilities import tools, intspan
from coverage import Coverage
from transcript import Transcript
from utilities.bam import BAM
from utilities.align_parsers import psl, sam

PACKAGE_DIR = "/".join(os.path.abspath(sys.argv[0]).split("/")[:-2])
CONFIG_FILE = PACKAGE_DIR + "/configs/model_matcher.cfg"
pp = pprint.PrettyPrinter(indent=0)

class MatchResult: 
    """Individual contig-to-transcript matching result"""
    mapping_columns = ['contig',
                       'contig_len',
                       'coord',
                       'model',
                       'transcript',
                       'gene',
                       'strand',
                       'coding',
                       'intronic',
                       'num_aln_blocks',
                       'num_exons',
                       'num_matched_blocks',
                       'matched_blocks',
                       'matched_exons',
                       'score',
                       'coverage',
                       'align_blocks'
                       ]

    def __init__(self, align_id, txt_id, num_blocks, num_exons, model=None, txt=None):
        self.align_id = align_id
        self.txt_id = txt_id
        self.num_blocks = num_blocks
        self.num_exons = num_exons
        self.txt = txt
        self.model = model
        self.matched_blocks = []
        self.matched_exons = []
        self.events = []
        self.coverage = 0.0
        self.model_weight = None
        self.first = self.last = self.middle = 0
        self.coord = ""
        self.within_intron = None
        self.overlap_exon = None
        
    def score(self):
        """Calculates score of a match result"""
        score = float(self.middle) + float(self.first) + float(self.last)
        for e in self.events:
            if e.event_type == "AS5" or e.event_type == "AS3":
                score += 0.5
            elif e.event_type == "AS53":
                score += 1.0
        
        return score

    def set_coverage(self, exon_bases, txt_len):
        """Determines transcript coverage given number of exonic bases covered"""
        self.coverage = float(exon_bases)/float(txt_len)

    @classmethod
    def compare_matches(cls, m1, m2):
        """For sorting matches"""
        if m1.score() > m2.score():
            return -1
        elif m1.score() < m2.score():
            return 1
        elif m1.num_blocks == m1.num_exons and m2.num_blocks != m2.num_exons:
            return -1
        elif m1.num_blocks != m1.num_exons and m2.num_blocks == m2.num_exons:
            return 1
        elif int(m1.middle) > int(m2.middle):
            return -1
        elif int(m1.middle) < int(m2.middle):
            return 1
        else:
            if int(m1.first) or int(m1.last) or int(m2.first) or int(m2.last):
                if int(m1.first) + int(m1.last) > int(m2.first) + int(m2.last):
                    return -1
                elif int(m1.first) + int(m1.last) < int(m2.first) + int(m2.last):
                    return 1
                elif m1.events or m2.events:
                    if m1.events and not m2.events:
                        return 1
                    elif not m1.events and m2.events:
                        return -1
                    elif len(m1.events) < len(m2.events):
                        return -1
                    elif len(m2.events) < len(m1.events):
                        return 1
                    else:
                        if m1.model_weight > m2.model_weight:
                            return -1
                        elif m1.model_weight < m2.model_weight:
                            return 1
                        else:
                            if m1.txt.cds_length() > m2.txt.cds_length():
                                return -1
                            elif m1.txt.cds_length() < m2.txt.cds_length():
                                return 1
                            return 0
                else:
                    if m1.model_weight > m2.model_weight:
                        return -1
                    elif m1.model_weight < m2.model_weight:
                        return 1
                    else:
                        if m1.txt.cds_length() > m2.txt.cds_length():
                            return -1
                        elif m1.txt.cds_length() < m2.txt.cds_length():
                            return 1
                        return 0
                    
            else:
                if m1.events and not m2.events:
                    return 1
                elif not m1.events and m2.events:
                    return -1
                elif len(m1.events) < len(m2.events):
                    return -1
                elif len(m2.events) < len(m1.events):
                    return 1
                else:
                    if m1.model_weight > m2.model_weight:
                        return -1
                    elif m1.model_weight < m2.model_weight:
                        return 1
                    else:
                        if m1.txt.cds_length() > m2.txt.cds_length():
                            return -1
                        elif m1.txt.cds_length() < m2.txt.cds_length():
                            return 1
                        return 0 
    
    def report(self, out):
        """Outputs match result"""
        cols = []
        
        report = {}
        report['contig'] = self.align_id
        report['contig_len'] = self.align.query_len
        report['coord'] = self.coord
        report['num_aln_blocks'] = self.num_blocks
        report['num_matched_blocks'] = len(self.matched_blocks)
        
        #alignment blocks
        blocks = []
        for block in self.align.blocks:
            blocks.append("%s,%s" % (block[0], block[1]))
        report['align_blocks'] = ';'.join(blocks)
        
        if self.matched_blocks:
            report['model'] = self.model
            report['transcript'] = self.txt.name
            if self.txt.alias is not None:
                report['gene'] = self.txt.alias
            report['strand'] = self.txt.strand
            report['coding'] = self.txt.coding_type()
            report['num_exons'] = self.num_exons
            report['matched_blocks'] = ','.join([str(i) for i in self.matched_blocks])
            
            if self.txt.strand == '+':
                report['matched_exons'] = ','.join([str(i) for i in self.matched_exons])
            else:
                exons = [len(self.txt.exons) - i + 1 for i in self.matched_exons]
                exons.sort(key=int)
                report['matched_exons'] = ','.join([str(e) for e in exons])
                
            report['score'] = "%.1f" % (self.score())
            report['coverage'] = "%.3f" % (self.coverage)
            
        elif self.within_intron is not None:
            report['model'] = self.model
            report['transcript'] = self.txt.name
            if self.txt.alias is not None:
                report['gene'] = self.txt.alias
            report['strand'] = self.txt.strand
            report['coding'] = self.txt.coding_type()
            report['num_exons'] = self.num_exons

            if self.txt.strand == '+':
                report['intronic'] = self.within_intron
            else:
                report['intronic'] = len(self.txt.exons) - self.within_intron

        cols = []
        for col in MatchResult.mapping_columns:
            if report.has_key(col):
                cols.append(str(report[col]))
            else:
                cols.append('-')
                        
        output = '\t'.join(cols) + '\n'
        if out:
            out.write(output)
        else:
            return output

class ModelMatcher:   
    """Matches contigs to transcripts through alignments to reference genome,
    identifies novel splicing events, and reports on coverage of genes and transcripts
    """
    def __init__(self, aligns=[], genome=None, model_order=[], 
                 assembly=None, outdir=None, refseq=None, coverage_model=None,
                 contigs_bam=None, genome_bam=None, min_spanning_reads=None, max_coverage_diff=None,
                 debug=False, annodir=None, configfile=None 
                 ):
        self.aligns = aligns
        self.genome = genome
        self.model_order = model_order
        self.debug = debug
        self.outdir = outdir
        self.min_spanning_reads = min_spanning_reads
        self.max_coverage_diff = max_coverage_diff
        self.assembly = assembly
        self.annodir = annodir
        self.configfile = configfile

        self.log_file = None
        if self.outdir is not None and os.path.isdir(self.outdir):
            self.log_file = outdir + "/log.txt"
        if self.log_file is not None:
            self.log = open(self.log_file, 'w')
        else:
            self.log = None

        self.chrom_proper = None
        if self.genome:
            self.get_models()
<<<<<<< HEAD
            self.chrom_proper = tools.ucsc_chroms(self.genome, self.annodir)
=======
            self.chrom_proper = tools.ucsc_chroms(self.genome)
>>>>>>> 76eb20f9024b0f5a3c90fc9f3aeb66a9e63c7b96

        # for event support
        self.contigs_bam = self.genome_bam = self.lib = None
        if contigs_bam:
            self.contigs_bam = BAM(contigs_bam, min_mapq=0)
        if genome_bam:
            self.use_genome = True
            self.genome_bam = BAM(genome_bam, min_mapq=1)
        else:
            self.use_genome = False

        self.refseq = refseq
        self.add_contig_reads = False
        
        # coverage
        if coverage_model and coverage_model in self.model_order:
            self.coverage_model = [coverage_model]
        else:
            self.coverage_model = self.model_order 
            
        self.matches = []
        self.events = []

    def get_models(self):
        """Extracts gene model info from config file info"""
        self.annots = {}
        self.indices = {}
        
        self.config = ConfigParser.ConfigParser()
        self.config.read(self.configfile)
        if self.genome in self.config.sections():
            for field in self.config.options(self.genome):
                value = self.config.get(self.genome, field)

                if field == "order" and not self.model_order:
                    self.model_order = value.split(',')
                else:
                    self.annots[field] = value

            if self.annots:
                for m,f in self.annots.iteritems():
                    full_path = os.path.join(self.annodir, self.genome, f)
                    if self.log:
                        self.log.write(m + ":" + full_path + "\n")
                    if os.path.exists(full_path):
                        self.annots[m] = full_path
                        self.indices[m] = full_path.replace(".txt", ".idx")

        self.splice_file = os.path.join(self.annodir, self.genome, 'splice_motifs.txt')
        self.splice_motifs = tools.get_splice_motifs(self.splice_file)

        if self.log:
            self.log.write("model order: %s\n" % ','.join(self.model_order))
            
    def match_models(self):
        """Main wrapper function of identifying events and calculating coverage"""
        sys.stdout.write("extracting transcript indices...")
        models = {}
        for i in range(len(self.model_order)):
            model = self.model_order[i]
            model_file = self.annots[model]
            if not self.annots.has_key(model) or not os.path.exists(model_file):
                if self.log:
                    self.log.write("doesn't have model: %s %s\n" % (model, model_file))
                continue
            models[model] = {}
            models[model]['source'] = self.annots[model]
            models[model]['weight'] = len(self.model_order) - i
            # extracts transcipt index info
            if os.path.exists(self.indices[model]):
                models[model]['index'] = self.extract_index(self.indices[model])
        sys.stdout.write("done\n")
        
        sys.stdout.write("starts matching %d alignments...\n" % len(self.aligns)) 
        # stores transcript objects - so it won't creates new objects for same trancript
        tobjs = {}
        events = []
        count = 1
        for align in self.aligns:
            # corrects alignment blocks if desired
            if self.correct_blocks and self.refseq and align.contig and align.contig.sequence:
                align.correct_blocks(self.splice_motifs, self.refseq, align.contig.sequence)
                align.merge_blocks()
            
            # finds overlapping txts
            overlapping_txts = []
            for model in self.model_order:
                if not tobjs.has_key(model):
                    tobjs[model] = {}
                
                index = models[model]['index']
                source = models[model]['source']
                txts = self.overlap_txts(tools.proper_chrom(align.target, chrom_proper=self.chrom_proper), align.blocks[0][0], align.blocks[-1][1], 
                                         index, source, model)
                for txt in txts:
                    if not txt:
                        continue
                    txt.model = model
                    txt.weight = models[model]['weight']
                    if not tobjs[model].has_key(txt.name):
                        tobjs[model][txt.name] = txt
                    overlapping_txts.append(tobjs[model][txt.name])

            matches = self.compare_model(align, overlapping_txts)            
            self.matches.extend(matches)
                        
            # gets best matches for each model and best match among all models
            best_matches = self.get_model_best_match(matches)
            best_match = None
            for model in self.model_order:
                if not best_matches.has_key(model):
                    continue
                best_match = best_matches[model][align.query]
                break
                     
            #capture events of best match
            if best_match and best_match.events:
                if len(best_match.matched_blocks) > 0:
                    events.extend(best_match.events)
                # novel transcript
                else:
                    events.append(best_match.events[0])
                remove_start = 1
            else:
                remove_start = 0

            # cleans up memory
            for i in range(remove_start, len(matches)):
                del matches[i].events[:]
            del matches[remove_start:]

            count += 1
            if count % 1000 == 0:
                sys.stderr.write("\rfinished %d alignments\n" % count)

        # get sequence for orf
        if self.assembly:
            contigs_with_events = dict((e.align.contig.num, e.align.contig) for e in events)
            self.assembly.get_seqs(contigs_with_events.values())

        # checks novelty and populates events
        event.Event.splice_motifs = self.splice_motifs
        self.events = [e for e in events if e.novel and not e.artefact]
        [e.details(self.refseq) for e in self.events]
        
        # output mapping file
        mapping_outfile = self.outdir + "/mapping.tsv"
        if os.path.exists(mapping_outfile):
            os.remove(mapping_outfile)
        self.output_mapping(mapping_outfile)
                        
        # get read support
        if self.contigs_bam: 
            sys.stderr.write("finding support\n")
            for e in self.events:
                event_support.find_support_contig(e, self.contigs_bam)
                
            if self.add_contig_reads:
                self.add_contig_support(self.events)
                
        # calculates and reports coverage
        coverage = self.calc_coverage(self.matches)
        coverage_outfile = self.outdir + "/coverage.tsv"
        coverage.output(coverage_outfile)

    def add_contig_support(self, all_events):
        """Sums contig read support for all members of same event"""
        grouped_events = self.group_events(all_events)
        
        for event_type in grouped_events.keys():
            # will not handle novel_transcript as number of blocks varies
            if event_type == 'novel_transcript':
                continue
            
            for coord, events in grouped_events[event_type].iteritems():
                if len(events) > 1:
                    if event_type in ('AS53', 'novel_exon', 'novel_utr'):
                        spanning_reads = [0, 0]
                    else:
                        spanning_reads = 0
                    coverage = 0
        
                    for e in events:        
                        if type(e.spanning_reads) == str and e.spanning_reads == 'na':
                            continue
                            
                        if e.spanning_reads != None and e.spanning_reads != 'na':
                            # spanning_reads += int(e.spanning_reads)
                            if event_type in ('AS53', 'novel_exon', 'novel_utr'):
                                spanning_reads[0] += int(e.spanning_reads[0])
                                spanning_reads[1] += int(e.spanning_reads[1])
                            else:
                                spanning_reads += int(e.spanning_reads)
                                
                        if type(e.coverage) == str and e.coverage == 'na':
                            continue
                        
                        if e.coverage != None:
                            coverage += int(e.coverage)                    

                    for e in events:
                        if event_type in ('AS53', 'novel_exon', 'novel_utr'):
                            e.spanning_reads = str(spanning_reads[0]) + ',' + str(spanning_reads[1])
                        else:
                            e.spanning_reads = spanning_reads
                            
                        if event_type in ('AS53', 'novel_exon', 'novel_utr', 'retained_intron'):
                            e.coverage = coverage
                                                
    def extract_index(self, index_file):
        """Extracts index info into dictionary"""
        index = {}
        for line in open(index_file, 'r'):
            coord, lines = line.rstrip().split(" ")
            index[coord] = lines
        return index

    def overlap_txts(self, target, tstart, tend, index, source, model):
        """Finds overlapping transcripts given coordinate, index, model"""
        tstart = int(tstart)
        tend = int(tend)
        tstart_index = ':'.join((target, str(int(tstart/1000))))
        tend_index = ':'.join((target, str(int(tend/1000))))

        txts = []
        txt_lines = {}
        if index.has_key(tstart_index):
            for line in index[tstart_index].split(','):
                txt_lines[line] = True
        if tend_index != tstart_index:
            int(int(tstart)/1000)
            for coord in range(int(tstart/1000)+1,int(tend/1000)+1):
                idx = ':'.join((target, str(coord)))
                if index.has_key(idx):
                    for line in index[idx].split(','):
                        txt_lines[line] = True

        for line_num in txt_lines.keys():
            line = linecache.getline(source, int(line_num))
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
            txts.append(txt)
            
        return txts

    def as_artefacts(self, match):
        """Removes matching AS5,AS3 and small novel exons"""
        novel_as5 = [e for e in match.events if e.novel and e.event_type == "AS5"]
        novel_as3 = [e for e in match.events if e.novel and e.event_type == "AS3"]

        as3_dict = {}
        for e in novel_as3:
            for block in e.align_blocks:
                as3_dict[int(block)] = e 
                
        as5_dict = {}
        for e in novel_as5:
            for block in e.align_blocks:
                as5_dict[int(block)] = e
                
        for e in novel_as5:
            for block in e.align_blocks:
                if as3_dict.has_key(int(block) - 1):
                    as3 = as3_dict[int(block) - 1]
                    e.artefact = True
                    as3.artefact = True

                # depends on transcript orientation, as5 and as3 can mean left or right
                elif as3_dict.has_key(int(block) + 1):
                    as3 = as3_dict[int(block) + 1]                                 
                    e.artefact = True
                    as3.artefact = True
                    
        # removes small novel exons <10 bp and associated with AS
        novel_exons = [e for e in match.events if e.novel and e.event_type == "novel_exon"]
        for e in novel_exons:
            for block in e.align_blocks:
                if intspan.cardinality(e.align.blocks[block-1]) < 10:
                    if block > 1 and block < len(e.align.blocks):
                        if as3_dict.has_key(block-1) and as3_dict[block-1].right:
                            e.artefact = True
                            as3_dict[block-1].artefact = True
                        if as5_dict.has_key(block-1) and as5_dict[block-1].right:
                            e.artefact = True
                            as5_dict[block-1].artefact = True
                        if as3_dict.has_key(block+1) and as3_dict[block+1].left:
                            e.artefact = True
                            as3_dict[block+1].artefact = True
                        if as5_dict.has_key(block+1) and as5_dict[block+1].left:
                            e.artefact = True
                            as5_dict[block+1].artefact = True
                        
                if e.artefact:
                    self.log.write(' '.join(("novel_exon artefact", e.contig, e.transcript, e.event_type, e.coordinate())) + "\n")

    def deletion_artefact(self, match, align):
        """Removes events associated with gaps <= 20bp"""
        # minimum intron size?
        min_intron_size = 20
        for e in match.events:
            if len(e.align_blocks) == 1:
                b = int(e.align_blocks[0])-1
                # check distance between previous block
                if b > 0 and int(align.blocks[b][0]) - int(align.blocks[b-1][1]) <= min_intron_size:
                    e.artefact = True
                    self.log.write(' '.join(("deletion artefact", e.contig, e.transcript, e.event_type, e.coordinate())) + "\n")
                # check distance between next block
                if b < len(align.blocks)-1 and int(align.blocks[b+1][0]) - int(align.blocks[b][1]) <= min_intron_size:
                    e.artefact = True
                    self.log.write(' '.join(("deletion artefact", e.contig, e.transcript, e.event_type, e.coordinate())) + "\n")
                    
            elif len(e.align_blocks) == 2:
                b1 = int(e.align_blocks[0])-1
                b2 = int(e.align_blocks[1])-1
                # check distance between blocks
                if int(align.blocks[b2][0]) - int(align.blocks[b1][1]) <= min_intron_size:
                    e.artefact = True
                    self.log.write(' '.join(("deletion artefact", e.contig, e.transcript, e.event_type, e.coordinate())) + "\n")

    def intron_artefact(self, match):
        """Removes novel_intron event if it's not in an exon"""
        for e in match.events:
            if e.event_type == "novel_intron":
                if re.match('[\D]', str(e.txt.cdsStart)) or re.match('[\D]',  str(e.txt.cdsEnd)):
                    continue

                start, end = e.align_coords[0][0], e.align_coords[-1][1]
                # novel intron
                if len(e.align_coords) > 1:
                    start, end = e.align_coords[0][1], e.align_coords[-1][0]

                if int(start) < int(e.txt.cdsStart) or int(end) > int(e.txt.cdsEnd):
                    e.artefact = True
                    self.log.write(' '.join(("novel_intron artefact", e.contig, e.transcript, e.event_type, e.coordinate())) + "\n")


    def too_many_exons_skipped(self, match):
        """Removes skipped exons if over 10 exons were skipped e.g. antibody gene"""
        max_skipped = 10  
        for e in match.events:
            if e.event_type == "skipped_exon" and len(e.exons) > max_skipped:
                e.artefact = True
                self.log.write(' '.join(("too-many-exons artefact", e.contig, e.transcript, e.event_type, e.coordinate())) + "\n")
                
    def check_highly_rearranged_gene(self, match):
        """Removes events of genes that are known to be highly rearranged"""
        bad_genes = ('abParts')
        if match.txt.alias and match.txt.alias in bad_genes:
            self.log.write('%s is mapped to a member of highly-rearranged genes, \"%s\" (model:%s) - ignore\n' % (match.align_id, match.txt.alias, match.model))
            for e in match.events:
                e.artefact = True
            
    def compare_model(self, align, txts):
        """Compares given alignment against given overlapping transcripts"""
        matches = []

        within_intron = []
        for txt in txts:
            result = self.match_exons(align.query, txt.full_name(), align.blocks, txt.exons, txt.chrom, strand=txt.strand)
            result.align = align
            result.coord = "%s:%s-%s" % (align.target, align.blocks[0][0], align.blocks[-1][1])
            result.model = txt.model
            result.model_weight = txt.weight
            result.txt = txt
                        
            if result.matched_blocks:
                if result.events:
                    for e in result.events:
                        e.txt = txt
                        e.align = align
                
                    if len(result.events) == 1 and result.events[0].event_type == 'novel_utr':
                        result.events[0].is_read_through(txts, self)
                                                                                               
                # get coverage
                exon_bases = intspan.intersect(align.blocks, txt.exons)
                result.set_coverage(exon_bases, txt.length)

            if result.matched_blocks:
                matches.append(result)
                
            if result.within_intron:
                within_intron.append(result)
                        
        if matches:
            if len(matches) > 1:
                matches.sort(MatchResult.compare_matches)

            for m in matches:
                if m.events:
                    novel_events = {}
                    for i in range(len(m.events)):
                        e = m.events[i]                        
                        novel_events_info = e.set_novelty(txts, matches)
                        
                        for ne in novel_events_info:
                            variant = event.Event(ne['type'])
                            variant.align = align
                            tools.set_attrs(variant, ne)

                            if not novel_events.has_key(i):
                                novel_events[i] = []
                            novel_events[i].append(variant)

                    replaced_events = novel_events.keys()
                    replaced_events.sort(lambda x,y: y-x)

                    for i in replaced_events:
                        del m.events[i]
                    for i in replaced_events:
                        m.events.extend(novel_events[i])

                    # test for validity of events
                    self.as_artefacts(m)
                    self.deletion_artefact(m, align)
                    self.intron_artefact(m)
                    self.too_many_exons_skipped(m)
                    self.check_highly_rearranged_gene(m)
                    
        elif within_intron:
            matches.extend(within_intron)

        else:
            result = MatchResult(align.query, None, len(align.blocks), 0, align)
            result.align = align
            result.coord = "%s:%s-%s" % (align.target, align.blocks[0][0], align.blocks[-1][1])
            result.model = 'no_match'
            print 'no match', align.query
            
            if len(align.blocks) >= 3:
                e = event.Event("novel_transcript")
                e.contig = align.query
                e.align = align
                e.chrom = align.target
                e.align_coords = align.blocks
                e.align_blocks = range(1,len(align.blocks)+1)
                result.events.append(e)
                
            matches.append(result)

        return matches
      
    def match_exons(self, align_id, txt_id, blocks, exons, chrom=None, strand=None):
        """Wrapper function for matching alignment blocks with exons"""
        if self.log:
            self.log.write("%s %d:\t%s\n" % (align_id, len(blocks), pp.pformat(blocks).replace("\n", " ")))
            self.log.write("%s %d:\t%s\n" % (txt_id, len(exons), pp.pformat(exons).replace("\n", " ")))

        result = MatchResult(align_id, txt_id, len(blocks), len(exons))
        result.perfect_matches = 0

        # finds start (first overlap of alignment block and exon) first
        b = 0
        e = 0
        b_start = None
        e_start = None
        
        # identifies starting matching block 0 -> when they overlap
        while b < len(blocks):
            block = blocks[b]
            while e < len(exons):
                exon = exons[e]
                
                # if not overlap, on to next exon
                if int(block[0]) != int(exon[0]) and int(block[1]) != int(exon[1]):
                    e += 1
                    continue
                # last block (but not only block) overlap first exon, not real, skip
                elif b == len(blocks) - 1 and e == 0 and len(blocks) > 1 and len(exons) > 1:
                    break
                # last exon overlap first block (but not only block), not real, skip
                elif e == len(exons) - 1 and b == 0 and len(blocks) > 1 and len(exons) > 1:
                    break

                #start identified, done get out of exon loop
                b_start = b
                e_start = e
                break

            # start not identified, next alignment block
            if b_start == None and e_start == None:
                e = 0
                b += 1
            # start identified, done get out of block loop
            else:
                break

        if b_start >= 0 and e_start >= 0:
            b = b_start
            e = e_start
            if self.log:
                self.log.write("start %d %d %d %d\n" % (b, e, len(blocks), len(exons)))
            if b_start > 0:
                if e_start == 0:
                    last_novel_utr_exon = None
                    for bb in range(0, b_start):
                        if int(blocks[bb][1]) < int(exons[0][0]):
                            last_novel_utr_exon = bb 
                        else:
                            if self.log:
                                self.log.write("novel intron %d %d\n" % (bb, bb+1))
                            
                            for i in range(bb + 1, bb + 3, 2):
                                variant = event.Event("novel_intron")
                                tools.set_attrs(variant, {'contig': align_id, 
                                                          'transcript': txt_id, 
                                                          'chrom': chrom,
                                                          'align_blocks': [i, i + 1],
                                                          'align_coords': blocks[i - 1 : i + 1],
                                                          'exons': [1], 
                                                          'exon_coords': [exons[0]]
                                                          }
                                                )
                                result.events.append(variant)
                                
                    if last_novel_utr_exon != None:
                        if self.log:
                            self.log.write("novel 5utr %s %s %s-%s\n" % (align_id, txt_id, blocks[0][0], blocks[last_novel_utr_exon][1]))

                        for i in range(1, last_novel_utr_exon + 2):
                            variant = event.Event("novel_utr")
                            tools.set_attrs(variant, {'contig': align_id, 
                                                      'transcript': txt_id, 
                                                      'chrom':chrom,
                                                      'align_blocks': [i],
                                                      'align_coords': [blocks[i - 1]],
                                                      'exons': [1],
                                                      'exon_coords':[exons[0]],
                                                      'last_matched_block': last_novel_utr_exon + 2,
                                                      }
                                            )
                            if strand == '+':
                                variant.prime = '5'
                            else:
                                variant.prime = '3'
                            result.events.append(variant)
                else:
                    for bb in range(1, b_start):
                        if self.log:
                            self.log.write("novel exon %d %d\n" % (bb,bb))
                        variant = event.Event("novel_exon")
                        tools.set_attrs(variant, {'contig': align_id,
                                                  'transcript': txt_id, 
                                                  'chrom': chrom,
                                                  'align_blocks': [bb + 1], 
                                                  'align_coords': [blocks[bb]],
                                                  'exons': [e_start], 
                                                  'exon_coords': [exons[e_start]]
                                                  }
                                        )
                        result.events.append(variant)

            while b < len(blocks):
                block = blocks[b]                        
                while e < len(exons):
                    exon = exons[e]
                    matched = self.match(blocks, exons, b, e)
                    if self.log:
                        self.log.write("%d %d %s %s %s\n" % (b, e, pp.pformat(block), pp.pformat(exon), matched))
                    if matched and (not matched == 'sync' or b == len(blocks)-1):
                        result.matched_exons.append(e+1)
                        result.matched_blocks.append(b+1)

                        if matched == 'perfect':
                            result.perfect_matches += 1

                        if not matched == 'perfect' and len(blocks) > 1:
                            if matched == 'as5' and not e == 0 and not b == 0 and not (b == 0 and int(block[0]) > int(exon[0])) and not (b > 0 and int(exon[0]) < int(blocks[b-1][-1])) and not (e > 0 and int(block[0]) < int(exons[e-1][-1])):
                                if self.log:
                                    self.log.write("5as %s %s %s-%s\n" % (align_id, txt_id, block[0], block[1]))

                                if strand == '-':
                                    variant = event.Event("AS3")
                                else:
                                    variant = event.Event("AS5")
                                tools.set_attrs(variant, {'left':True, 'right':False})
                                tools.set_attrs(variant, {'contig':align_id, 'transcript':txt_id, 'chrom':chrom, 'edge':'left'})
                                tools.set_attrs(variant, {'align_blocks':[b + 1], 'align_coords':blocks[b]})
                                tools.set_attrs(variant, {'exons':[e + 1], 'exon_coords':exons[e]})
                                result.events.append(variant)

                            if matched == 'as3' and not e == len(exons)-1 and not b == len(blocks)-1:
                                if self.log:
                                    self.log.write("3as %s %s %s-%s\n" % (align_id, txt_id, block[0], block[1]))
                                if strand == '-':
                                    variant = event.Event("AS5")
                                else:
                                    variant = event.Event("AS3")
                                tools.set_attrs(variant, {'left':False, 'right':True})
                                tools.set_attrs(variant, {'contig':align_id, 'transcript':txt_id, 'chrom':chrom, 'edge':'right'})
                                tools.set_attrs(variant, {'align_blocks':[b + 1], 'align_coords':blocks[b]})
                                tools.set_attrs(variant, {'exons':[e + 1], 'exon_coords':exons[e]})
                                result.events.append(variant)

                            if matched == 'as53' and not e == 0 and not b == 0 and not (b == 0 and int(block[0]) > int(exon[0])) and not (b > 0 and int(exon[0]) < int(blocks[b-1][-1])) and not (e > 0 and int(block[0]) < int(exons[e-1][-1])) and not e == len(exons)-1 and not b == len(blocks)-1:
                                if self.log:
                                    self.log.write("as53 %d %d\n" % (b, b+1))
                                variant = event.Event("AS53")
                                tools.set_attrs(variant, {'contig':align_id, 'transcript':txt_id, 'chrom':chrom})
                                tools.set_attrs(variant, {'align_blocks':[b + 1], 'align_coords':[blocks[b]]})
                                tools.set_attrs(variant, {'exons':[e + 1], 'exon_coords':[exons[e]]})
                                result.events.append(variant)

                        # record matches
                        if e == 0 or e == len(exons)-1:
                            if matched == 'perfect':
                                score = 2
                            # 'as53' will have no score
                            elif matched == 'as5' or matched == 'as3':
                                score = 1
                            else:
                                score = 0

                            if e == 0:
                                result.first = score
                            else:
                                result.last = score
                        else:
                            if matched == 'perfect':
                                result.middle += 2
                            # start(5') is different, end(3') is the same
                            elif matched == 'as5' and b == 0:
                                result.middle += 1
                            # start(5') is the same, end(3') is different
                            elif matched == 'as3' and b == len(blocks)-1:
                                result.middle += 1
                        b += 1
                        e += 1
                        break

                    # novel exon
                    elif int(block[1]) < int(exon[0]):
                        # increment block, until finding matching block/exon
                        bb = b + 1
                        matched = False
                        while not matched and bb <= len(blocks)-1:
                            matched = self.match(blocks, exons, bb, e)
                            if not matched:
                                bb += 1

                        # if bb == len(blocks)-1, check if loop stops because of match or last index
                        # if loops stops because of finding matching exon, then novel exons is between b and bb
                        if matched:
                            if self.log:
                                self.log.write("novel exons %d %d\n" % (b, bb))
                            for i in range(b + 1, bb + 1):
                                variant = event.Event("novel_exon")
                                tools.set_attrs(variant, {'contig': align_id, 
                                                          'transcript': txt_id, 
                                                          'chrom': chrom,
                                                          'align_blocks': [i],
                                                          'align_coords': [blocks[i - 1]],
                                                          'exons': [e, e + 1],
                                                          'exon_coords':exons[e - 1 : e + 1],
                                                          }
                                                )
                                result.events.append(variant)
                                
                        # set next comparison to be bb, if it's the last index, then handle it outside??, it will handle alt.splicing
                        b = bb
                        break

                    # skipped exon
                    elif int(block[0]) > int(exon[1]):
                        # increment exon, until finding matching block/exon
                        ee = e + 1
                        matched = False
                        while not matched and ee <= len(exons)-1:
                            matched = self.match(blocks, exons, b, ee)
                            if not matched:
                                ee += 1
                        # requires ee to be matched
                        if matched:                       
                            if self.log:
                                self.log.write("skipped exons %d %d\n" % (e, ee))
                            for i in range(e + 1, ee + 1):
                                variant = event.Event("skipped_exon")
                                tools.set_attrs(variant, {'contig' : align_id,
                                                          'transcript' : txt_id,
                                                          'chrom': chrom,
                                                          'exons': [i],
                                                          'exon_coords': [exons[i - 1]],
                                                          'align_blocks': range(b, b + 2),
                                                          'align_coords': blocks[b - 1 : b + 1],
                                                          }
                                                )
                                result.events.append(variant)
                        e = ee
                        break

                    # retained intron
                    elif int(block[1]) > int(exon[1]):
                        #increment exon, until finding matching block/exon
                        ee = e + 1
                        if ee <= len(exons)-1 and int(exons[ee][0]) < int(block[1]):
                            matched = False
                            while (not matched or matched == 'sync') and ee <= len(exons)-1:
                                matched = self.match(blocks, exons, b, ee)
                                if not matched or matched == 'sync':
                                    ee += 1
                            if matched and not matched == "as3":
                                if self.log:
                                    self.log.write("retained intron %d %d\n" % (e, ee))
                                for i in range(e + 1, ee + 1):
                                    variant = event.Event("retained_intron")
                                    tools.set_attrs(variant, {'contig': align_id, 
                                                              'transcript': txt_id, 
                                                              'chrom':chrom,
                                                              'align_blocks': [b + 1], 
                                                              'align_coords': [blocks[b]],
                                                              'exons': [i, i + 1],
                                                              'exon_coords': exons[i - 1 : i + 1],
                                                              }
                                                    )
                                    result.events.append(variant)
                        else:
                            b = b + 1
                        e = ee
                        break

                    # novel intron/indel
                    elif int(block[1]) < int(exon[1]):
                        # increment block, until finding matching block/exon
                        bb = b + 1
                        if bb <= len(blocks)-1 and int(blocks[bb][0]) < int(exon[1]):
                            matched = False
                            while not matched and bb <= len(blocks)-1:
                                matched = self.match(blocks, exons, bb, e)
                                if not matched:
                                    bb += 1
                            if matched:
                                if self.log:
                                    self.log.write("novel intron %d %d\n" % (b, bb))
                                for i in range(b + 1, bb + 2, 2):
                                    variant = event.Event("novel_intron")
                                    tools.set_attrs(variant, {'contig': align_id, 
                                                              'transcript': txt_id, 
                                                              'chrom': chrom,
                                                              'align_blocks': [i, i + 1],
                                                              'align_coords': blocks[i - 1 : i + 1],
                                                              'exons': [e + 1], 
                                                              'exon_coords':[exons[e]]
                                                              }
                                                    )
                                    result.events.append(variant)
                        else:
                            e = e + 1
                        b = bb
                        break

                    else:
                        if self.log:
                            self.log.write("problem?\n")
                        e += 1
                        break
                    break
                
                # make sure the loop ends
                if e > len(exons)-1:
                    break
                                
            # if last matched/synced is not last block and is last exon
            if b <= len(blocks)-1:
                for bb in range(b, len(blocks)):
                    if int(blocks[bb][0]) > int(exons[e-1][1]):
                        if self.log:
                            self.log.write("novel 3utr %s %s %d %d %s-%s\n" % (align_id, txt_id, bb, e, exons[0][0], exons[-1][-1]))
                        for i in range(bb + 1, len(blocks) + 1):
                            variant = event.Event("novel_utr")
                            tools.set_attrs(variant, {'contig': align_id,
                                                      'transcript': txt_id,
                                                      'chrom': chrom,
                                                      'align_blocks': [i],
                                                      'align_coords': [blocks[i - 1]],
                                                      'exons': [len(exons)],
                                                      'exon_coords':[exons[-1]],
                                                      'last_matched_block': bb,
                                                      }
                                            )
                            if strand == '-':
                                variant.prime = '5'
                            else:
                                variant.prime = '3'
                            result.events.append(variant)
                        break
                    else:
                        if self.log:
                            self.log.write("novel intron %d %d\n" % (bb-1, bb))
                            for i in range(bb, bb + 2, 2):
                                variant = event.Event("novel_intron")
                                tools.set_attrs(variant, {'contig': align_id, 
                                                          'transcript': txt_id, 
                                                          'chrom':chrom,
                                                          'align_blocks': [i, i + 1],
                                                          'align_coords': blocks[i - 1 : i + 1],
                                                          'exons':[e], 
                                                          'exon_coords': [exons[e - 1]]
                                                          }
                                                )
                                result.events.append(variant)
                
        else:
            inside_exon = None
            inside_intron = None
            overlapped_exons = []
            for e in range(len(exons)):
                # entire alignment within a single exon
                if intspan.subsume([blocks[0][0], blocks[-1][1]], exons[e]):
                    inside_exon = e + 1
                    
                    for b in range(len(blocks)):
                        result.matched_blocks.append(b + 1)
                    result.matched_exons.append(inside_exon)
                
                    if self.log:
                        self.log.write("%s within exon %s of %s" % (align_id, inside_exon, txt_id))
                    break
                
                # entire alignment within a single intron
                if e < len(exons) - 1 and intspan.subsume([blocks[0][0], blocks[-1][1]], [int(exons[e][1]) + 1, int(exons[e + 1][0]) - 1]):
                    inside_intron = e + 1
                    result.within_intron = inside_intron
                    
                    if self.log:
                        self.log.write("%s within intron %s of %s" % (align_id, inside_intron, txt_id))
                    
                    break
                                    
            if not inside_exon and not inside_intron:        
                self.log.write("cannot find start\n")

        if self.log:
            self.log.write("----\n")

        if result.events and not result.matched_blocks:
            del result.events[:]
                     
        return result

    def match(self, blocks, exons, b, e, move=None):
        """Matches alignment blocks and exons"""
        match = False

        if int(blocks[b][0]) == int(exons[e][0]) and int(blocks[b][1]) == int(exons[e][1]):
            match = "perfect"

        # single block and embedded within intron
        elif len(blocks) == 1:
            # one block embedded within an exon
            if intspan.subsume(blocks[b], exons[e]):
                match = "perfect"
            elif int(blocks[b][1]) > int(exons[e][1]):
                match = "as3"
            elif not e == 0 and int(blocks[b][0]) < int(exons[e][0]):
                match = "as5"

        elif len(exons) == 1:
            if intspan.subsume(blocks[b], exons[e]):
                match = "sync"

        elif e == len(exons)-1:
            # if not overlap, or cut into next exon return false
            if intspan.overlap(blocks[b], exons[e]):
                #last block, last exon
                if b == len(blocks)-1:
                    if int(blocks[b][0]) == int(exons[e][0]):
                        match = "perfect"
                    else:
                        match = "as5"
                # not last block, novel exon
                else:
                    # 2nd condition check if it has retained intron
                    if int(blocks[b][0]) == int(exons[e][0]) or int(blocks[b][0]) < int(exons[e-1][1]):
                        match = "as3"
                    else:
                        match = "as5"

        elif b == len(blocks)-1:
            # if not overlap, or cut into next exon return false
            # first check overlap
            if intspan.overlap(blocks[b], exons[e]):
                # last block vs middle exon
                if e < len(exons)-1:
                    # not retained intron
                    if int(blocks[b][1]) < int(exons[e+1][0]):
                        if int(blocks[b][0]) == int(exons[e][0]) or int(blocks[b][0]) < int(exons[e-1][1]):
                            match = "as3"
                        elif int(blocks[b][-1]) == int(exons[e][-1]):
                            match = "as5"
                        else:
                            match = "as53"
                    # retained intron
                    else:
                        match = "sync1"
                # last block vs last exon (last exon = only check as5)
                elif not (int(blocks[b][0]) == int(exons[e][0])):
                    match = "as5"
                # perfect = 3' end not match
                else:
                    match = "perfect"
        
        elif intspan.subsume(blocks[b], exons[e]) or intspan.subsume(exons[e], blocks[b]):
            if int(exons[e][1]) < int(blocks[b+1][0]) and int(blocks[b][1]) < int(exons[e+1][0]):
                if int(blocks[b][0]) == int(exons[e][0]):
                    match = "as3"
                elif int(blocks[b][1]) == int(exons[e][1]):
                    match = "as5"
                else:
                    if int(blocks[b+1][0]) > int(exons[e][1]):
                        match = "as53"
                    else:
                        match = "sync"
            else:
                match = "sync"

        # as5 or as3
        elif int(blocks[b][1]) < int(exons[e+1][0]) and int(exons[e][1]) < int(blocks[b+1][0]):
            if int(blocks[b][0]) == int(exons[e][0]):
                match = "as3"
            elif int(blocks[b][1]) == int(exons[e][1]):
                match = "as5"
            elif intspan.overlap(blocks[b], exons[e]):
                match = "as53"
            else:
                match = "sync"
                    
        elif int(blocks[b][1]) < int(exons[e+1][0]) and int(exons[e][1]) < int(blocks[b + 1][0]) and\
             int(exons[e + 1][0]) < int(blocks[b + 1][1]) and int(blocks[b + 1][0]) < int(exons[e + 1][1]):
            if int(blocks[b][0]) == int(exons[e][0]):
                match = "as3"
            elif int(blocks[b][1]) == int(exons[e][1]):
                match = "as5"
            else:
                match = "as53"
                
        # move on
        elif int(blocks[b][0]) > int(exons[e][1]) and int(blocks[b][1]) < int(exons[e+1][0]):
            match = "sync"
            
        else:
            pass

        return match
            
    def get_model_best_match(self, matches):
        """Determines best match within same model"""
        best_matches = {}
        for match in matches:                
            if not best_matches.has_key(match.model):
                best_matches[match.model] = {}
                
            if not best_matches[match.model].has_key(match.align_id):
                best_matches[match.model][match.align_id] = match
            
            elif float(best_matches[match.model][match.align_id].coverage) < float(match.coverage):
                best_matches[match.model][match.align_id] = match
                                
        return best_matches
             
    def calc_coverage(self, matches):
        """Calculates total coverage given list of matches"""
        best_matches = self.get_model_best_match(matches)
        
        annot_files = {}
        for model in self.model_order:
            annot_files[model] = self.annots[model]
            
        coverage = Coverage(best_matches, self.coverage_model, annot_files=annot_files)
        coverage.process()
        
        return coverage
        
    def output_mapping(self, outfile):
        """Outputs mapping results"""
        out = open(outfile, 'w')
        tools.write_header(MatchResult.mapping_columns, out)
        best_matches = self.get_model_best_match(self.matches)
        
        ordered_by_contig = {}
        for model in self.model_order:
            if not best_matches.has_key(model):
                continue
            for contig, match in best_matches[model].iteritems():
                if not ordered_by_contig.has_key(contig):
                    ordered_by_contig[contig] = []
                ordered_by_contig[contig].append(match)
                                
        for contig in sorted(ordered_by_contig.keys()):
            for match in ordered_by_contig[contig]:
                match.report(out=out)
                
        if best_matches.has_key('no_match'):
            for contig, match in best_matches['no_match'].iteritems():
                match.report(out=out)
                                
        out.close()
                                            
    def group_events(self, events):
        """Groups list of events by event type and coordinate into dictionary"""
        grouped_events = {}
        for e in events:
            if e.novel and not e.artefact and e.event_type != 'novel_transcript':
                if e.genome_coord is not None:
                    coord = e.genome_coord
                else:
                    coord = e.coordinate()
                    
                if not grouped_events.has_key(e.event_type):
                    grouped_events[e.event_type] = {coord:[e]}
                else:
                    if not grouped_events[e.event_type].has_key(coord):
                        grouped_events[e.event_type][coord] = [e]
                    else:
                        grouped_events[e.event_type][coord].append(e)

        self.group_novel_events([e for e in events if e.event_type == 'novel_transcript'], grouped_events)

        return grouped_events

    def group_novel_events(self, events, grouped_events):
        """Groups novel transcripts"""
        novel_transcripts = {}
        count = 0
        for e in events:
            if e.genome_coord is not None:
                chrom, start, end = re.split('[:-]', e.genome_coord)
            else:
                chrom, start, end = re.split('[:-]', e.coordinate())
                                         
            if not novel_transcripts.has_key(chrom):
                novel_transcripts[chrom] = {}
                novel_transcripts[chrom][start + "-" + end] = [e]
            else:
                found = False
                for start_end in sorted(novel_transcripts[chrom].keys()):
                    ee = novel_transcripts[chrom][start_end][0]

                    if intspan.same_blocks(e.align_coords, ee.align_coords) or intspan.same_blocks(ee.align_coords, e.align_coords):
                        novel_transcripts[chrom][start_end].append(e)
                        found = True
                        break

                if not found:
                    if not novel_transcripts[chrom].has_key(start + "-" + end):
                        novel_transcripts[chrom][start + "-" + end] = [e]
                    else:
                        novel_transcripts[chrom][start + "-" + end + "." + str(count)] = [e]

            count += 1

        event_type = 'novel_transcript'
 
        if len(novel_transcripts.keys()) > 0:
            grouped_events[event_type] = {}
        
        for chrom in novel_transcripts.keys():
            for start_end in novel_transcripts[chrom].keys():
                coord = chrom + ":" + start_end
                
                if not grouped_events[event_type].has_key(coord):
                    grouped_events[event_type][coord] = []
                
                for e in novel_transcripts[chrom][start_end]:
                    grouped_events[event_type][coord].append(e)

    def output_events(self, grouped_events, output_file, summary_file=None):
        """Outputs events"""
        output = open(output_file, 'w')
        event_counts = {}

        headers = event.Event.headers[:]
        headers.insert(0, 'id')
        headers.extend(event.Event.headers_support)
        output.write('\t'.join(headers) + "\n")

        count1 = 1
        for event_type in sorted(grouped_events):
            event_counts[event_type] = 0
            
            for coord in grouped_events[event_type].keys():
                event_counts[event_type] += 1
                count2 = 1
                for e in grouped_events[event_type][coord]:
                    if len(grouped_events[event_type][coord]) == 1:
                        count = count1
                    else:
                        count = "%d.%d" % (count1, count2)

                    if self.debug:
                        self.log.write("outputting event %s for contig %s\n" % (e.event_type, e.contig))

                    if e.genome_coord:
                        output.write(str(count) + "\t" + e.output(parsed=True) + "\n")
                    else:
                        output.write(str(count) + "\t" + e.output() + "\n")

                    count2 += 1
                count1 += 1

        output.close()

        # summary file showing number of events
        if summary_file is not None:
            summary = open(summary_file, 'w')
            for type in sorted(event_counts.keys()):
                summary.write(type + '\t' + str(event_counts[type]) + '\n')
            summary.close()
        
    def parse_events(self, events_file):
        """Wrapper function for parsing events of single events file"""
        for line in open(events_file, 'r'):
            if re.search('type', line):
                continue
            
            e = event.Event.parse(line)
            self.events.append(e)
            
    def parse_events_dir(self, path, outdir):
        """Wrapper function for parsing a directory of events"""
        output_dirs = sorted(glob.glob(os.path.join(path, '*')))
                
        events_files = []
        mapping_files = []
        coverage_files = []
        num_output_dirs = 0
        missing_dirs = []
        
        ok = True
        
        for job_num in range(1, len(output_dirs)+1):
            cluster_outdir = "%s/%s" % (args[0], job_num)
            if os.path.isdir(cluster_outdir):
                num_output_dirs += 1
                    
                events_file = cluster_outdir + '/events.tsv'                    
                if os.path.exists(events_file):
                    events_files.append(events_file)
                else:
                    missing_dirs.append(str(job_num))
                    sys.stdout.write("%s does not have events file\n" % (cluster_outdir))
                    ok = False
                
                mapping_file = cluster_outdir + '/mapping.tsv'
                if os.path.exists(mapping_file):
                    mapping_files.append(mapping_file)
                elif not str(job_num) in missing_dirs:
                    missing_dirs.append(str(job_num))
                    sys.stdout.write("%s does not have mapping file\n" % (cluster_outdir))
                    ok = False
                    
                coverage_file = cluster_outdir + '/coverage.tsv'
                if os.path.exists(mapping_file):
                    coverage_files.append(coverage_file)
                elif not str(job_num) in missing_dirs:
                    missing_dirs.append(str(job_num))
                    sys.stdout.write("%s does not have coverage file\n" % (cluster_outdir))
                    ok = False
                        
        sys.stdout.write("output dirs:%s events files:%s mapping files:%s coverage_files:%s\n" % 
                         (num_output_dirs, len(events_files), len(mapping_files), len(coverage_files)))
                
        if num_output_dirs == len(events_files):
            concat_outfile = outdir + "/mapping.tsv"
            tools.concat_tsv(mapping_files, concat_outfile)
            
            concat_outfile = outdir + "/events.tsv"
            tools.concat_tsv(events_files, concat_outfile)
            
            for events_file in events_files:
                self.parse_events(events_file)
                
            coverage_outfile = outdir + "/coverage.tsv"
            all_results = []
            coverage = Coverage(models=self.coverage_model)
            for coverage_file in coverage_files:
                results = coverage.parse_results(coverage_file)
                all_results.extend(results)
            
            coverage.combine_results(all_results)
            coverage.output(coverage_outfile, no_blocks=True)
            
        else:
            sys.stdout.write("missing (%s):%s\n" % (len(missing_dirs), ','.join(missing_dirs)))
            ok = False

        return ok
            
def main(args, options):
    outdir = args[1]
    # outputs log file for command run, parameters, etc
    tools.output_log(__version__, sys.argv, [vars(options)], outdir)
            
    if len(args) == 2:
        models = None
        if options.models:
            models = options.models.split(',')
            
        if not options.results:
            if options.log:
                log_file = outdir + "/log.txt"

            refseq = None
            if options.ref:
<<<<<<< HEAD
                refseq = tools.get_refseq_from_2bit(options.annodir, options.genome)
=======
                refseq = tools.get_refseq_from_2bit(options.genome)
>>>>>>> 76eb20f9024b0f5a3c90fc9f3aeb66a9e63c7b96

            # do this before extracting alignments because need to get splice file info
            mm = ModelMatcher(genome=options.genome, 
                              model_order=models,
                              refseq=refseq, 
                              contigs_bam=options.contigs_bam, 
                              genome_bam=options.genome_bam,
                              min_spanning_reads=options.spanning_reads, 
                              max_coverage_diff=options.coverage_diff, 
                              coverage_model=options.coverage_model,
                              outdir=outdir, 
                              debug=options.debug,
                              annodir=options.annodir,
                              configfile=options.config_file 
                              )
            mm.add_contig_reads = options.add_contig_reads
            mm.correct_blocks = options.fix_align
        
            # get alignments
            sys.stderr.write("extracting alignments...")
            ext = os.path.splitext(args[0])[1]
            filters = {'unique':True, 'bestn':1, 'match':90, 'identity':0}
            if ext == ".psl":
                mm.aligns = psl.parse(args[0], filters, splice_motif_file=mm.splice_file, refseq=refseq)
            elif ext == ".sam":
                mm.aligns = sam.parse(args[0], filters, splice_motif_file=mm.splice_file, refseq=refseq)
            sys.stderr.write(str(len(mm.aligns)) + '\n')

            # get contig sequences
            contigs = None
            ass = None
            if options.fasta:
                sys.stderr.write("getting contig sequences...")
                contigs_with_aligns = [a.query for a in mm.aligns]
                mm.assembly = assembly.Assembly(None, fasta=options.fasta)
                contigs = mm.assembly.get_contigs(ids=contigs_with_aligns, sequence=True)
                sys.stderr.write("done\n")

                # links alignment to contig
                if contigs:
                    contig_dict = dict((c.num, c) for c in contigs)
                    for align in mm.aligns:
                        if contig_dict.has_key(align.query):
                            align.contig = contig_dict[align.query]

            # match against transcript models
            if mm.model_order:
                mm.match_models()

        # events file(s) given
        else:
            mm = ModelMatcher(genome=options.genome, 
                              model_order=models, 
                              coverage_model=options.coverage_model,
                              annodir=options.annodir,
                              configfile=options.config_file
                              )
            if os.path.isdir(args[0]):
                ok = mm.parse_events_dir(args[0], outdir)
                if not ok:
                    sys.exit(1)
            else:
                mm.parse_events(args[0])
                
            if options.add_contig_reads:
                mm.add_contig_support(mm.events)
                                
        # outputs events
        events_outfile = outdir + "/events.tsv"
        grouped_events = mm.group_events(mm.events)
        mm.output_events(grouped_events, events_outfile)

        if options.filter:
            for e in mm.events:
                event_support.set_read_support(e, min_spanning_reads=options.spanning_reads, max_coverage_diff=options.coverage_diff)
                event_support.screen(e)
            events_passed = [e for e in mm.events if e.filter_result == 'passed']
            events_outfile = outdir + "/events_filtered.tsv"
            events_summary = outdir + "/events_summary.tsv"
            grouped_events_passed = mm.group_events(events_passed)
            mm.output_events(grouped_events_passed, events_outfile, events_summary)
            
    else:
        parser.error("incorrect number of arguments")
  
if __name__ == '__main__':
    usage = "Usage: %prog alignment_file|results output_dir"

    parser = OptionParser(usage=usage, version="%prog " + __version__)
    parser.add_option("-g", "--genome", dest="genome", help="genome") 
    parser.add_option("-m", "--models", dest="models", help="model order e.g. k,e,r,a, if not specified, model order specified in config file will be used")
    parser.add_option("-f", "--sequence", dest="fasta", help="fasta file")
    parser.add_option("-r", "--ref", dest="ref", help="use reference sequence for deducing splice sites for psl", action="store_true", default=False)
    parser.add_option("-C", "--contigs_bam", dest="contigs_bam", help="reads to contigs bam file")    
    parser.add_option("-B", "--genome_bam", dest="genome_bam", help="reads to genome bam file")
    parser.add_option("-s", "--add_contig_reads", dest="add_contig_reads", action="store_true", default=False)
    parser.add_option("-M", "--coverage_model", dest="coverage_model", help="single model used for reporting coverage e.g. e")
    parser.add_option("-X", "--filter", dest="filter", help="filter", action="store_true", default=False)
    parser.add_option("-N", "--results", dest="results", help="results files given", action="store_true", default=False)
    parser.add_option("-l", "--log", dest="log", action="store_true", help="output log file", default=False)
    parser.add_option("-d", "--debug", dest="debug", action="store_true", default=False)
    parser.add_option("-A", "--fix_align", dest="fix_align", help="fix alignment block coords", action="store_true", default=False)
    parser.add_option("-S", "--spanning_reads", dest="spanning_reads", help="minimum number of spanning reads. Default=2", type='int', default=2)
    parser.add_option("-R", "--coverage_diff", dest="coverage_diff", help="maximum coverage diff(ratio) between event block and neighbors. Default=100", type='int', default=100)
    parser.add_option("--annodir", dest="annodir", help="the Trans-ABySS 'annotation' directory")
    parser.add_option("--mmcfg", dest="config_file", help="the path of `model_matcher.cfg'")
    
    (options, args) = parser.parse_args()
    main(args, options)
