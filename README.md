# Trans-ABySS
## *De novo* assembly of RNAseq data using ABySS

[Ka Ming Nip](mailto:kmnip@bcgsc.ca), [Readman Chiu](mailto:rchiu@bcgsc.ca)

Copyright 2015 Canada's Michael Smith Genome Sciences Centre, BC Cancer Agency

Please use our [Google Group](mailto:trans-abyss@googlegroups.com) for [discussions and
support](https://groups.google.com/d/forum/trans-abyss).
  
You may also create [issues](https://github.com/bcgsc/transabyss/issues) on our GitHub repository.

If you use Trans-ABySS, please cite:

[Robertson, G., et al. 2010. De novo assembly and analysis of RNA-seq data. Nature Methods 7, 909-912(2010)](http://www.nature.com/nmeth/journal/v7/n11/full/nmeth.1517.html)

--------------------------------------------------------------------------------

Program requirements for `transabyss` and `transabyss-merge`:
  * [ABySS 1.5.1+](https://github.com/bcgsc/abyss/releases)
  * [Python 2.7.6+](https://www.python.org/download/releases/2.7.6/)
  * [python-igraph 0.7.0+](http://igraph.org/python/#downloads)
  * [BLAT](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat)

Program requirements for `transabyss-analyze`:
  * [Python 2.7.6+](https://www.python.org/download/releases/2.7.6/)
  * [BLAT](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat)
  * [Pysam](https://github.com/pysam-developers/pysam)
  * [BioPython](http://biopython.org/wiki/Download)
  * [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
  * [GMAP/GSNAP](http://research-pub.gene.com/gmap/)
  * [Samtools](https://github.com/samtools/samtools)
  * reference genome and annotations for the organism of interest

Required Python packages (ie. python-igraph, Pysam, BioPython) can be installed
easily with `pip`, ie.

```
  pip install python-igraph
  pip install pysam
  pip install biopython
```

Other required softwares must be accessible from your `PATH` environment variable.

To test `transabyss` on our sample dataset:

```
  bash sample_dataset/assemble.sh
```  
  
To test `transabyss-analyze` on our sample dataset:

```
  bash sample_dataset/analyze.sh
```

Please see [TUTORIAL.txt] for more information on usage of each application.


--------------------------------------------------------------------------------
EOF
