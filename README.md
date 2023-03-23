# GLAR SUNK Scaffolding Pipeline

This repository contains a snakemake to assist targeted assembly of complex regions using singly unqiue nucleotide k-mers.
This version was specifically used for the targeted assembly of Centromeres in CHM1 (a haploid genome)
You can find a version for diploid assemblies in the branch "diploid"

# Input requires:
- hifi_assembly: haploid HiFi assembly
- hifi_fofn: file with paths to HiFi fastq files
- ont_fofn: file with paths to ONT fastq files
- contig_bed: bed file with contigs to scaffold

# Output:
PAF of ONT reads that link two Contigs that have SUNK placement support (>2 SUNK Matches between ONT read and assembly and alignment length > 50kbp)
Output will be available along this path:
'''
{sample}/join_reads/{kmer_size}/{region}/{map}-{remap}/joining_reads.paf 
'''
