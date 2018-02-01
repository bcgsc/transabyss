# This make file would install all require software packages except ABySS.

SHELL:=/bin/bash -eu -o pipefail

.DELETE_ON_ERROR:

WORKDIR:=$(shell pwd)

IGRAPH_VERSION:=0.7.1
PYTHON_VERSION:=2.7.10

all: transabyss transabyss-merge
	touch $@.COMPLETE

clean: clean_pip clean_blat clean_python clean_igraph
	rm -f all.COMPLETE && \
	rm -rf ./bin ./lib ./include ./share

transabyss: bin blat.COMPLETE python.COMPLETE python-igraph.COMPLETE

transabyss-merge: transabyss

bin:
	mkdir -p ./bin

blat.COMPLETE: bin
	wget --quiet --no-check-certificate -N -P ./bin http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
	chmod +x ./bin/blat && \
	touch $@
	
clean_blat:
	rm -f blat.COMPLETE ./bin/blat

python.COMPLETE: bin
	wget --quiet --no-check-certificate -N https://www.python.org/ftp/python/$(PYTHON_VERSION)/Python-$(PYTHON_VERSION).tgz && \
	tar -zxf Python-$(PYTHON_VERSION).tgz && \
	cd Python-$(PYTHON_VERSION) && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@

clean_python:
	rm -f python.COMPLETE python-igraph.COMPLETE Python-$(PYTHON_VERSION).tgz && \
	rm -rf Python-$(PYTHON_VERSION)

pip.COMPLETE: bin python.COMPLETE
	wget --quiet --no-check-certificate -N https://bootstrap.pypa.io/get-pip.py && \
	./bin/python get-pip.py && \
	touch $@

clean_pip:
	rm -f pip.COMPLETE get-pip.py

igraph.COMPLETE: bin
	wget --quiet --no-check-certificate -N http://igraph.org/nightly/get/c/igraph-$(IGRAPH_VERSION).tar.gz && \
	tar -zxf igraph-$(IGRAPH_VERSION).tar.gz && \
	cd igraph-$(IGRAPH_VERSION) && \
	./configure --quiet --prefix=$(WORKDIR) && \
	make --quiet && \
	make --quiet install && \
	cd $(WORKDIR) && \
	touch $@

clean_igraph:
	rm -f igraph.COMPLETE igraph-$(IGRAPH_VERSION).tar.gz && \
	rm -rf igraph-$(IGRAPH_VERSION)

python-igraph.COMPLETE: bin python.COMPLETE pip.COMPLETE igraph.COMPLETE
	export PKG_CONFIG_PATH=$(WORKDIR)/lib/pkgconfig && \
	./bin/pip --quiet install python-igraph && \
	touch $@
