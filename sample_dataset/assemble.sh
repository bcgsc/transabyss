#!/bin/bash

set -euo pipefail

cd $(dirname $0)

name='test'
reads1=./reads/rnaseq_1.fq.gz
reads2=./reads/rnaseq_2.fq.gz
assemblydir=./${name}/assembly/
finalassembly=${assemblydir}/${name}-final.fa

# set up the environment
TRANSABYSS_PATH=$(readlink -f ..)
export PATH=${TRANSABYSS_PATH}:${TRANSABYSS_PATH}/bin/:${PATH}

# assemble the test dataset
transabyss --se ${reads1} ${reads2} --outdir ${assemblydir} --name ${name} --threads 2 --island 0 -c 1

# check completeness of the assembly
if [ ! -e ${assemblydir}/${name}.FINAL.COMPLETE ] || [ $(cat ${finalassembly} | wc -l) -lt 2 ]
then
    echo 'Assembly is incomplete!'
    exit 1
fi

#EOF
