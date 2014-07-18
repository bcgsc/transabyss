Trans-ABySS - de novo assembly of RNAseq data using ABySS.

Ka Ming Nip <kmnip@bcgsc.ca>, Readman Chiu <rchiu@bcgsc.ca>, Anthony Raymond <traymond@bcgsc.ca>
Copyright 2014 Canada's Michael Smith Genome Sciences Centre, BC Cancer Agency

Please use our Google Group <trans-abyss@googlegroups.com> for discussions and
support. Existing topics can be viewed at:
  <https://groups.google.com/d/forum/trans-abyss>
  
You may also create issues on our GitHub repository at:
  <https://github.com/bcgsc/transabyss/issues>

                                  ~ README ~
================================================================================
VERSION 1.5.1 (July 18, 2014)

+ Support strand-specific RNAseq data.
+ Major improvement in assembly quality; both rare and common transcripts can be
  assembled well with "small" k-mer sizes (ie. 25~32).
+ Package has been divided into 3 main applications:
  (1) transabyss          - RNAseq assembler at a single k-mer size
  (2) transabyss-merge    - merge multiple assemblies from (1)
  (3) transabyss-analyze  - analyze an assembly, either from (1) or (2), for
                            structural variants and novel splice variants
- Analysis results are NOT screened against dbSNP, DGV, etc. anymore.
- Genome assembly and analyses pipeline has been retracted.
- Support for SGE qmake has been retracted.

Program requirements for 'transabyss' and 'transabyss-merge':
  * ABySS 1.5.1+          <https://github.com/bcgsc/abyss/releases>
  * Python 2.7.6          <https://www.python.org/download/releases/2.7.6/>
  * python-igraph 0.7.0+  <http://igraph.org/python/#downloads>
  * BLAT                  <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat>

Program requirements for 'transabyss-analyze':
  * Python 2.7.6  <https://www.python.org/download/releases/2.7.6/>
  * BLAT          <http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat>
  * Pysam         <http://code.google.com/p/pysam/>
  * BioPython     <http://biopython.org/wiki/Download>
  * Bowtie2       <http://bowtie-bio.sourceforge.net/bowtie2/index.shtml>
  * GMAP/GSNAP    <http://research-pub.gene.com/gmap/>
  * Samtools      <http://sourceforge.net/projects/samtools/files/samtools/>
  * reference genome and annotations for the organism of interest

Required Python packages (python-igraph, Pysam, BioPython) can be installed
easily with pip, ie.

  pip install igraph
  pip install pysam
  pip install biopython

Other required softwares must be accessible from your PATH environment variable.

To test `transabyss' on our sample dataset:

  bash sample_dataset/assemble.sh
  
To test `transabyss-analyze' on our sample dataset:

  bash sample_dataset/analyze.sh

Please see TUTORIAL.txt for more information on the usage of each application.


================================================================================
                                    ~ EOF ~
