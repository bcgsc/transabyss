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

analysesdir=./${name}/analyses/
antisensefusionsfile=${analysesdir}/fusions/antisense_fusion.tsv
indelsfile=${analysesdir}/indels/events_filtered.tsv
splicefile=${analysesdir}/splicing/events_filtered.tsv

# analyze the assembly for structural variants
transabyss-analyze -a ${finalassembly} -1 ${reads1} -2 ${reads2} -n ${name} --outdir ${analysesdir} --ref hg19

# check completeness of fusions stage
if [ ! -e ${analysesdir}/fusions/${name}.FUSIONS.COMPLETE ] || [ $(cat ${antisensefusionsfile} | wc -l) -lt 2 ]
then
    echo 'Fusions stage is incomplete!'
    exit 1
else
    echo "fusion:    " $(grep 'antisense_fusion' ${antisensefusionsfile} |cut -f 13)
fi

# check completeness of indels stage
if [ ! -e ${analysesdir}/indels/${name}.INDELS.COMPLETE ] || [ $(cat ${indelsfile} | wc -l) -lt 2 ]
then
    echo 'Indels stage is incomplete!'
    exit 1
else
    echo "insertion: " $(grep 'ins' ${indelsfile} |cut -f 3,4,5)
    echo "deletion:  " $(grep 'del' ${indelsfile} |cut -f 3,4,5)
fi

# check completeness of splicing stage
if [ ! -e ${analysesdir}/splicing/${name}.SPLICING.COMPLETE ] || [ $(cat ${splicefile} | wc -l) -lt 2 ]
then
    echo 'Splices stage is incomplete!'
    exit 1
else
    echo "skip-exon: " $(grep 'skipped_exon' ${splicefile} |cut -f 8)
fi



#EOF
