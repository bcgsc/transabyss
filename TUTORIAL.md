The Trans-ABySS package has 2 main applications:  
1. **transabyss** - RNAseq assembler at a single k-mer size  
2. **transabyss-merge** - merge multiple assemblies from (1)  


### Content:  
**PART 1. Assembly with `transabyss`**  
+ [Usage](#11-usage)  
+ [Running transabyss](#12-running-transabyss)  
+ [Input reads](#13-input-reads)  
+ [Output assembly path](#14-output-assembly-path)  
+ [k-mer size](#15-k-mer-size)  
+ [MPI and multi-threading](#16-mpi-and-multi-threading)  
+ [Strand-specific assembly](#17-strand-specific-assembly)  

**PART 2. Merging assemblies with `transabyss-merge`**  
+ [Usage](#21-usage)   
+ [Running transabyss-merge](#22-running-transabyss-merge)  
+ [Minimum and maximum k-mer sizes](#23-minimum-and-maximum-k-mer-sizes)  
+ [Output merged assembly path](#24-output-merged-assembly-path)  
+ [Output sequence prefixes](#25-output-sequence-prefixes)  

Please use our Google Group <trans-abyss@googlegroups.com> for discussions and support. Existing topics can be viewed at <https://groups.google.com/d/forum/trans-abyss>.
  
You may also create issues on our GitHub repository at <https://github.com/bcgsc/transabyss/issues>.

----

#### PART 1. Assembly with `transabyss`

##### [1.1] Usage  

    transabyss --help    

    
##### [1.2] Running transabyss  
  
    transabyss --se reads.fq  

See below for assembling single-end reads and/or pair-end reads.


##### [1.3] Input reads  

Supported formats are compressed (gzip/bzip) FASTQ/FASTA/SAM/QSEQ or BAM files.  
Paired-end reads in FASTQ/A can be interleaved or in separate files.
Paired-end reads in FASTQ/A must have the suffixes `/1` or `/2` in the read name.
  
Use options `--se` and `--pe` to specify the path(s) of single-end reads and paired-end reads, respectively. Example usages:

combination | reads for sequence content | reads for paired-end linkage | assembly
----|----|----|----
`--se r.fq` | r.fq | _none_ | single-end
`--pe r.fq` | r.fq | r.fq | paired-end
`--se SE.fq --pe PE.fq` | SE.fq | PE.fq | paired-end
`--se SE.fq PE.fq --pe PE.fq` | SE.fq, PE.fq | PE.fq | paired-end
                
                                 
##### [1.4] Output assembly path
  
default:  `./transabyss_1.5.3_assembly/transabyss-final.fa`
  
Use options `--name` and `--outdir` to change the output directory and assembly name.
  
  
##### [1.5] k-mer size

default: `32`
  
Use option `--kmer` to adjust the k-mer size.
k=32 has a good trade-off for assembling both rare and common transcripts.
Using larger k-mers improve the assembly quality of common transcripts and transcripts with repetitive regions, but the assembly of rare transcripts may suffer.
  
  
##### [1.6] MPI and multi-threading

default:  no mpi processes; singe-threaded
    
Use option `--threads` to specify the number of threads.
Use option `--mpi` to specify the number of MPI processes.
  
Only the first stage of assembly (de Bruijn graph) could benefit from MPI; all remaining stages may be multi-threaded.


##### [1.7] Strand-specific assembly

Use option `--SS` to indicate that input reads are strand-specific.
Strand-specific reads are expected to have `/1` reads in _reverse_ direction and `/2` reads in _forward_ direction, ie.
```    
                    <-- R/1
 5' =======================> 3' transcript
    F/2 -->
```
----
  
#### PART 2. Merging assemblies with `transabyss-merge`
  
Should you choose to assemble the same dataset with different settings (different k-mer sizes), you can merge the assemblies together into one FASTA file for downstream analyses (with `transsabyss-analyze`). When a sequence is contained in a longer sequence, the longer sequence is kept.

##### [2.1] Usage

    transabyss-merge --help


##### [2.2] Running transabyss-merge
    
Example:
  
    transabyss-merge a.fa b.fa --mink 32 --maxk 64


##### [2.3] Minimum and maximum k-mer sizes
  
`--mink` and `--maxk` are used to specific the smallest and largest k-mer sizes in the input assemblies.
  
  
##### [2.4] Output merged assembly path
  
default:  `./transabyss-merged.fa`
  
Use option `--out` to specify the output path.


##### [2.5] Output sequence prefixes
  
Use option `--prefixes` to specify the prefixes of output sequences.
  
Example:
  
    transabyss-merge a.fa b.fa c.fa --mink 32 --maxk 64 --prefix k32. k48. k64.

file | prefix  
-----|--------  
a.fa | k32.  
b.fa | k48.  
c.fa | k64.  
  
One prefix for sequences from each input assembly FASTA file. This feature helps you keep track of the origin of each seqeunce in the merged assembly.

----
