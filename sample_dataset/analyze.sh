#!/bin/bash

set -eo pipefail

cd $(dirname $0)

name='test'
reads1=./reads/rnaseq_1.fq.gz
reads2=./reads/rnaseq_2.fq.gz
assemblydir=./${name}/assembly
finalassembly=${assemblydir}/${name}-final.fa

# set up the environment
TRANSABYSS_PATH=$(readlink -f ..)
export PATH=${TRANSABYSS_PATH}:${TRANSABYSS_PATH}/bin/:${PATH}
if [[ -z "$PYTHONPATH" ]]
then
    # Set the Trans-ABySS code base directory as PYTHONPATH
    export PYTHONPATH=${TRANSABYSS_PATH}
else
    # Prepend the Trans-ABySS code base directory to PYTHONPATH
    export PYTHONPATH=${TRANSABYSS_PATH}:$PYTHONPATH
fi

# assemble the test dataset
transabyss --se ${reads1} ${reads2} --outdir ${assemblydir} --name ${name} --threads 2 --island 0 -c 1

# check completeness of the assembly
if [ ! -e ${assemblydir}/${name}.FINAL.COMPLETE ] || [ $(cat ${finalassembly} | wc -l) -lt 2 ]
then
    echo 'ERROR: Assembly is incomplete!'
    exit 1
fi

analysesdir=./${name}/analyses
antisensefusionsfile=${analysesdir}/fusion/antisense_fusion.tsv
indelsfile=${analysesdir}/indel/events_filtered.tsv
splicefile=${analysesdir}/splice/events_filtered.tsv

# analyze the assembly for structural variants
transabyss-analyze -a ${finalassembly} -1 ${reads1} -2 ${reads2} -n ${name} --outdir ${analysesdir} --ref hg19

# check completeness of fusion stage
if [ ! -e ${analysesdir}/fusion/${name}.FUSION.COMPLETE ] || [ $(cat ${antisensefusionsfile} | wc -l) -lt 2 ]
then
    echo 'ERROR: No fusion events!'
    #exit 1
else
    echo "fusion:    " $(grep 'antisense_fusion' ${antisensefusionsfile} |cut -f 13)
fi

# check completeness of indel stage
if [ ! -e ${analysesdir}/indel/${name}.INDEL.COMPLETE ] || [ $(cat ${indelsfile} | wc -l) -lt 2 ]
then
    echo 'ERROR: No indel events!'
    #exit 1
else
    echo "insertion: " $(grep 'ins' ${indelsfile} |cut -f 3,4,5)
    echo "deletion:  " $(grep 'del' ${indelsfile} |cut -f 3,4,5)
fi

# check completeness of splice stage
if [ ! -e ${analysesdir}/splice/${name}.SPLICE.COMPLETE ] || [ $(cat ${splicefile} | wc -l) -lt 2 ]
then
    echo 'ERROR: No splice events!'
    #exit 1
else
    echo "skip-exon: " $(grep 'skipped_exon' ${splicefile} |cut -f 8)
fi



#EOF
