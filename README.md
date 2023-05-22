# SUNK-based contig scaffolding pipeline

This repository contains a snakemake to assist in the targeted assembly of complex, repetitive regions using singly unique nucleotide k-mers (SUNKs). This version was specifically used for the targeted assembly of centromeres in CHM1 (an effectively haploid genome). You can find a version for diploid assemblies in the branch “diploid”.

# Input requires:
- hifi_assembly: haploid HiFi assembly
- hifi_fofn: file with paths to HiFi fastq files
- ont_fofn: file with paths to ONT fastq files
- contig_bed: bed file with contigs to scaffold

# Output:
PAF of ONT reads that links two contigs that have SUNK placement support (>2 SUNK matches between ONT read and assembly, and alignment length >50 kbp). Output will be available along this path:
```
{sample}/join_reads/{kmer_size}/{region}/{map}-{remap}/joining_reads.paf 
```
# To Run on Eichler Lab Cluster
```
snakesub -j 40
```

Currently operates using the Eichler lab module system.  Currently updating snakemake to use conda envs for off-cluster availability
 
