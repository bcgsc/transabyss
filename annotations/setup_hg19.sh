#!/bin/bash

# This script is for setting up the human reference genome hg19 fasta, gmap indexes, and annotations.
#
# Setup Requirements:
# * connection to the internet.
# * Python

set -eo pipefail

skipanno=false
skipref=false

while getopts ":arh" flag; do
    case $flag in
        a)
            skipanno=true
            ;;
        r)
            skipref=true
            ;;
        h)
            echo "Usage: $0 [-a|-r|-h]"
            echo "       -a   do not set up annotation files"
            echo "       -r   do not set up reference genome fasta and indexes"
            echo "       -h   print this help and exit"
            exit
            ;;
        *)
            echo "Invalid option: -$OPTARG" >&2
            exit
            ;;
    esac
done


genome='hg19'
scriptdir=$(dirname $0)
annodir=$(readlink -f ${scriptdir})
codebase=$(readlink -f ${scriptdir}/..)

annofileslist=${annodir}/${genome}_anno_files.txt
genomedir=${annodir}/${genome}
refgenomefasta=${genomedir}/${genome}.fa

UCSC_ANNO_URL=ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/database/
UCSC_REF_URL=ftp://hgdownload.cse.ucsc.edu/goldenPath/${genome}/bigZips/
UCSC_EXE_URL=http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
GMAP_INDEX_URL=http://research-pub.gene.com/gmap/index.html


# Set up the environment variables
if [[ -z "$PYTHONPATH" ]]
then
    # Set the Trans-ABySS code base directory as PYTHONPATH
    export PYTHONPATH=${codebase}
else
    # Prepend the Trans-ABySS code base directory to PYTHONPATH
    export PYTHONPATH=${codebase}:$PYTHONPATH
fi
echo PYTHONPATH: $PYTHONPATH
# Prepend the Trans-ABySS bin directory to PATH
export PATH=${codebase}/bin/:$PATH
echo PATH: $PATH

set -x

# Create the output directory
mkdir -p ${genomedir}
cd ${genomedir}


if [ "$skipanno" = true ]; then
    echo "Skipping setup of annotations files ..."
else
    # Download annotation files
    wget --quiet -N -nH -x --cut-dirs=3 -P ${genomedir} -B ${UCSC_ANNO_URL} -i ${annofileslist}

    # KNOWNGENE
    if [ ! -f knownGene_ref.idx ] || [ knownGene_ref.idx -ot knownGene.txt.gz ] || [ knownGene_ref.idx -ot kgXref.txt.gz ] ; then
        join <(zcat knownGene.txt.gz |sort) <(zcat kgXref.txt.gz |sort) -t $'\t' >knownGene_ref.txt
        python ${codebase}/analysis/annotations/knownGene.py knownGene_ref.txt -i knownGene_ref.idx.incomplete -g ${genome}
        mv -f knownGene_ref.idx.incomplete knownGene_ref.idx
    fi

    # ENSEMBL
    if [ ! -f ensGene_ref.idx ] || [ ensGene_ref.idx -ot ensGene.txt.gz ] || [ ensGene_ref.idx -ot ensemblToGeneName.txt.gz ] ; then
        join -1 2 -2 1 <(zcat ensGene.txt.gz |sort -k2) <(zcat ensemblToGeneName.txt.gz |sort) -a 1 >ensGene_ref.txt
        python ${codebase}/analysis/annotations/ensembl.py ensGene_ref.txt -i ensGene_ref.idx.incomplete -g ${genome}
        mv -f ensGene_ref.idx.incomplete ensGene_ref.idx
    fi

    # ACEVIEW
    if [ ! -f acembly_ref.idx ] || [ acembly_ref.idx -ot acembly.txt.gz ]; then
        zcat acembly.txt.gz |awk -v FS="[\.\t]" '{print $0"\t"$2}' >acembly_ref.txt
        python ${codebase}/analysis/annotations/aceview.py acembly_ref.txt -i acembly_ref.idx.incomplete -g ${genome}
        mv -f acembly_ref.idx.incomplete acembly_ref.idx
    fi 

    # REFSEQ
    if [ ! -f refGene.idx ] || [ refGene.idx -ot refGene.txt.gz ]; then
        zcat refGene.txt.gz >refGene.txt
        python ${codebase}/analysis/annotations/refGene.py refGene.txt -i refGene.idx.incomplete -g ${genome}
        mv -f refGene.idx.incomplete refGene.idx
    fi

    # cytoband
    if [ ! -f cytoBand.txt ] || [ cytoBand.txt -ot cytoBand.txt.gz ]; then
        zcat -df cytoBand.txt.gz >cytoBand.txt.incomplete
        mv -f cytoBand.txt.incomplete cytoBand.txt
    fi

    # splice motifs
    cp ${codebase}/annotations/shared/splice_motifs.txt .
    
    echo "Set up of ${genome} annotations completed."
fi


if [ "$skipref" = true ]; then
    echo "Skipping the setup of reference fasta and indexes ..."
else
    # Download the 2bit file
    wget --quiet -N -P ${genomedir} ${UCSC_REF_URL}/${genome}.2bit

    # Convert the 2bit to FASTA
    if [ ! -f ${refgenomefasta} ] || [ ${refgenomefasta} -ot ${genomedir}/${genome}.2bit ]; then
        if [ ! $(command -v twoBitToFa) ]; then
            # download twoBitToFa
            wget --quiet -N -P ${codebase}/bin ${UCSC_EXE_URL}/twoBitToFa
            chmod +x ${codebase}/bin/twoBitToFa
            export PATH=${codebase}/bin/twoBitToFa:$PATH
        fi

        twoBitToFa ${genomedir}/${genome}.2bit ${refgenomefasta}.incomplete
        mv -f ${refgenomefasta}.incomplete ${refgenomefasta}
    fi

    # Build GMAP indexes
    if [ ! $(command -v gmap_build) ]; then
        # Extract the indentifier of GMAP's latest version from the index page
        wget --quiet -N -P ${genomedir} ${GMAP_INDEX_URL}
        gmap_version=$(grep -i -m 1 'latest version' ${genomedir}/index.html |sed -e 's/.*version //g' -e 's/<\/a>//g')
        
        # Download and install GMAP
        wget --quiet -N -P ${codebase}/bin http://research-pub.gene.com/gmap/src/gmap-gsnap-${gmap_version}.tar.gz
        tar -C ${codebase}/bin -zxf ${codebase}/bin/gmap-gsnap-${gmap_version}.tar.gz
        cd ${codebase}/bin/gmap-${gmap_version}/ && ${codebase}/bin/gmap-${gmap_version}/configure --prefix=${codebase}/bin/gmap-${gmap_version} && make --quiet && make --quiet install
        
        ln -sf ${codebase}/bin/gmap-${gmap_version}/bin/{gmap,gsnap,gmap_build} ${codebase}/bin/
        export PATH=${codebase}/bin/gmap_build:$PATH
        rm -f ${genomedir}/index.html
    fi
    gmap_build -d ${genome} ${refgenomefasta} -D ${genomedir}
    
    # Correct the reference path in Trans-ABySS configuration file
    reference=$(echo ${refgenomefasta} |sed -e 's/\//\\\//g')
    sed -i -e "s/^${genome}:.*/${genome}: ${reference}/g" ${codebase}/configs/transcriptome.cfg
    
    echo "Set up of ${genome} reference fasta and indexes completed."
fi

#EOF
