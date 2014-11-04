# This make file would install all require software packages except ABySS.

SHELL:=/bin/bash -eu -o pipefail

.DELETE_ON_ERROR:

WORKDIR:=$(shell pwd)
MAX_KMER:=96

ABYSS_VERSION:=1.5.2
BOWTIE2_VERSION:=2.2.3
GMAP_VERSION:=2014-01-21
IGRAPH_VERSION:=0.7.1
PYTHON_VERSION:=2.7.8
PYSAM_VERSION:=0.7.8
SAMTOOLS_VERSION:=0.1.19


all: transabyss transabyss-merge transabyss-analyze

transabyss: bin blat.COMPLETE python.COMPLETE python-igraph.COMPLETE
#transabyss: bin blat.COMPLETE python.COMPLETE python-igraph.COMPLETE abyss.COMPLETE

transabyss-merge: transabyss

transabyss-analyze: bin blat.COMPLETE bowtie2.COMPLETE python.COMPLETE samtools.COMPLETE biopython.COMPLETE pysam.COMPLETE gmap.COMPLETE

bin:
	mkdir -p ./bin

samtools.COMPLETE: bin
	wget --quiet --no-check-certificate -O samtools-$(SAMTOOLS_VERSION).zip -N https://github.com/samtools/samtools/archive/$(SAMTOOLS_VERSION).zip && \
	unzip ./samtools-$(SAMTOOLS_VERSION).zip && \
	cd ./samtools-$(SAMTOOLS_VERSION) && \
	make --quiet && \
	cd $(WORKDIR) && \
	ln -sf ../samtools-$(SAMTOOLS_VERSION)/samtools ./bin/samtools && \
	touch $@

blat.COMPLETE: bin
	wget --quiet --no-check-certificate -N -P ./bin http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
	chmod +x ./bin/blat && \
	touch $@

python.COMPLETE: bin
	wget --quiet --no-check-certificate -N https://www.python.org/ftp/python/$(PYTHON_VERSION)/Python-$(PYTHON_VERSION).tgz && \
	tar -zxf Python-$(PYTHON_VERSION).tgz && \
	cd Python-$(PYTHON_VERSION) && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@

pip.COMPLETE: bin python.COMPLETE
	wget --quiet --no-check-certificate -N https://bootstrap.pypa.io/get-pip.py && \
	./bin/python get-pip.py && \
	touch $@

igraph.COMPLETE: bin
	wget --quiet --no-check-certificate -N http://igraph.org/nightly/get/c/igraph-$(IGRAPH_VERSION).tar.gz && \
	tar -zxf igraph-$(IGRAPH_VERSION).tar.gz && \
	cd igraph-$(IGRAPH_VERSION) && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@

python-igraph.COMPLETE: bin python.COMPLETE pip.COMPLETE igraph.COMPLETE
	export PKG_CONFIG_PATH=$(WORKDIR)/lib/pkgconfig && \
	./bin/pip --quiet install python-igraph && \
	touch $@

biopython.COMPLETE: bin python.COMPLETE pip.COMPLETE
	./bin/pip --quiet install biopython && \
	touch $@

pysam.COMPLETE: bin python.COMPLETE pip.COMPLETE
	./bin/pip --quiet install pysam==$(PYSAM_VERSION) && \
	touch $@

bowtie2.COMPLETE: bin
	wget --quiet --no-check-certificate -N -O bowtie2-$(BOWTIE2_VERSION)-linux-x86_64.zip http://sourceforge.net/projects/bowtie-bio/files/bowtie2/$(BOWTIE2_VERSION)/bowtie2-$(BOWTIE2_VERSION)-linux-x86_64.zip/download && \
	unzip ./bowtie2-$(BOWTIE2_VERSION)-linux-x86_64.zip && \
	ln -sf ../bowtie2-$(BOWTIE2_VERSION)/{bowtie2-build,bowtie2} ./bin && \
	touch $@

gmap.COMPLETE: bin
	#wget --quiet --no-check-certificate -N http://research-pub.gene.com/gmap/index.html && \
	#gmap_version=`grep -i -m 1 'latest version' ./index.html |sed -e 's/.*version //g' -e 's/<\/a>//g'` && \
	gmap_version=$(GMAP_VERSION) && \
	wget --quiet --no-check-certificate -N http://research-pub.gene.com/gmap/src/gmap-gsnap-$$gmap_version.tar.gz && \
	tar -zxf ./gmap-gsnap-$$gmap_version.tar.gz && \
	cd ./gmap-$$gmap_version && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@

boost.COMPLETE: bin
	wget --quiet --no-check-certificate -N -O boost_1_55_0.tar.gz http://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download && \
	tar -zxf ./boost_1_55_0.tar.gz && \
	cd ./boost_1_55_0 && \
	./bootstrap.sh --prefix=$(WORKDIR) && \
	./b2 install && \
	cd $(WORKDIR) && \
	touch $@

sparsehash.COMPLETE: bin
	wget --quiet --no-check-certificate -N http://sparsehash.googlecode.com/files/sparsehash-2.0.2.tar.gz && \
	tar -zxf ./sparsehash-2.0.2.tar.gz && \
	cd ./sparsehash-2.0.2 && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@	

openmpi.COMPLETE: bin
	wget --quiet --no-check-certificate -N http://www.open-mpi.org/software/ompi/v1.8/downloads/openmpi-1.8.1.tar.gz && \
	tar -zxf ./openmpi-1.8.1.tar.gz && \
	cd ./openmpi-1.8.1 && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet all install && \
	cd $(WORKDIR) && \
	touch $@	
	
abyss.COMPLETE: bin boost.COMPLETE sparsehash.COMPLETE openmpi.COMPLETE
	wget --quiet --no-check-certificate -N https://github.com/bcgsc/abyss/releases/download/$(ABYSS_VERSION)/abyss-$(ABYSS_VERSION).tar.gz && \
	tar -zxf ./abyss-$(ABYSS_VERSION).tar.gz && \
	cd ./abyss-$(ABYSS_VERSION) && \
	./configure --quiet --prefix=$(WORKDIR) --enable-maxk=$(MAX_KMER) --with-boost=$(WORKDIR)/include --with-mpi=$(WORKDIR)/include CPPFLAGS=-I$(WORKDIR)/include && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@

