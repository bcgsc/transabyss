#!/bin/bash

set -euo pipefail

cd $(dirname $0)

kmer1=25
kmer2=32
name1=test.k${kmer1}
name2=test.k${kmer2}
reads1=./reads/rnaseq_1.fq.gz
reads2=./reads/rnaseq_2.fq.gz
assemblydir1=./${name1}
assemblydir2=./${name2}
finalassembly1=${assemblydir1}/${name1}-final.fa
finalassembly2=${assemblydir2}/${name2}-final.fa
mergedassembly=./merged.fa

# set up the environment
TRANSABYSS_PATH=$(readlink -f ..)
export PATH=${TRANSABYSS_PATH}:${TRANSABYSS_PATH}/bin/:${PATH}

# assemble the test dataset
transabyss -k ${kmer1} --se ${reads1} ${reads2} --outdir ${assemblydir1} --name ${name1} --threads 2 --island 0 -c 1

transabyss -k ${kmer2} --se ${reads1} ${reads2} --outdir ${assemblydir2} --name ${name2} --threads 2 --island 0 -c 1

transabyss-merge --mink ${kmer1} --maxk ${kmer2} --prefixes k${kmer1}. k${kmer2}. --out ${mergedassembly} ${finalassembly1} ${finalassembly2}

#EOF
