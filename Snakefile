import pandas as pd
import os
import numpy as np
import pysam
import gzip
from Bio import SeqIO
import re

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)

configfile: 'config.yaml'

HIFI_FOFN = config['hifi_fofn']
ONT_FOFN = config['ont_fofn']
KMER_LIST = config['kmer_size']
CONTIG_LIST = config['contig_list']
LEN_FILT=config['length_filter_ONT']
HIFI_ASSEMBLY = config['hifi_assembly']
TYPE_MAP = config['init_map']
REMAP = config['remap']
SAMPLE = config['sample_name']
UPPER_LIM = int(config['upper_lim'])
LOWER_LIM = int(config['lower_lim'])
HIFI_FAI = config['hifi_assembly'] +".fai"
COV=config['hifi_cov']
MONOMERS=config['monomers']
CONTIG_BED=config['contig_bed']

HIFI_FOFN_DF = pd.read_csv(HIFI_FOFN, sep='\t', header=None, names=['file'])
ONT_FOFN_DF = pd.read_csv(ONT_FOFN, sep='\t', header=None, names=['file'])


fofn_dict = {'HiFi' : HIFI_FOFN_DF, 'ONT' : ONT_FOFN_DF}

shell.prefix("source %s/env.cfg; " % (SNAKEMAKE_DIR))

wildcard_constraints:
	part = "[0-9]{4}",
	tech='|'.join(['ONT', 'HiFi']),
	kmer='|'.join(['%d' % size for size in KMER_LIST])

def findHiFiReads(wildcards):
	df = HIFI_FOFN_DF
	return df.at[int(wildcards.read), 'file']

def findONTReads(wildcards):
	df = ONT_FOFN_DF
	return df.at[int(wildcards.read), 'file']

def jellyCounts(wildcards):
	return expand(rules.jellyfish_count.output.tally, sample=wildcards.sample, kmer=wildcards.kmer, read=HIFI_FOFN_DF.index)

def gatherPartsONT(wildcards):
	PARTS = glob_wildcards("{sample}/reads/tmp/ONT/temp_ONT_{read}_{part}".format(sample=wildcards.sample, read=wildcards.read, part='{part}')).part
	return expand(rules.ont_fastq_to_fasta.output.ont_fasta, sample=wildcards.sample, read=wildcards.read, part=PARTS)

def findONTAlign(wildcards):
	PARTS = glob_wildcards("{sample}/reads/tmp/ONT/temp_ONT_{read}_{part}".format(sample=wildcards.sample, read=wildcards.read, part='{part}')).part
	return expand("{sample}/alignments/ONT/{map}/{read}_{part}.bam", sample=wildcards.sample, part=PARTS, map=wildcards.map, read=wildcards.read)

### EDIT THIS FUNCTION ###
def findRegion(wildcards):
	contig_df = pd.read_csv(CONTIG_LIST.get((wildcards.region)), sep ='\t', header=None, names=['file'])
	contig_list = contig_df["file"].tolist()
	string_contig = " ".join(contig_list)
	print(string_contig)
	return (string_contig)

def getContig(wildcards):
	contig_file = CONTIG_BED.get((wildcards.region))
	return(contig_file)

def aggregate_input(wildcards):
	check_file = checkpoints.get_depth_cov.get(sample=wildcards.sample, region=wildcards.region).output[0] 
	if os.stat(check_file).st_size == 0:
		return "{sample}/break_contigs/{region}/contigs.bed"
	else:
		return "{sample}/break_contigs/{region}/final.bed"
def getMonomer(wildcards):
	return MONOMERS.get(wildcards.region)

def get_ONT_per_split(wildcards):
	ONT_READS = glob_wildcards("{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/{ontread}/correction.fa".format(sample=wildcards.sample, split=wildcards.split, kmer=wildcards.kmer, region=wildcards.region, map=wildcards.map, map2=wildcards.map2, ontread='{ontread}')).ontread
	return expand(rules.string_decomposer_ONT.output.tsv, sample=wildcards.sample, map=TYPE_MAP, map2=REMAP, region=wildcards.region, kmer=wildcards.kmer, split=wildcards.split, ontread=ONT_READS)
def gatherAll(wildcards):
	SPLIT = glob_wildcards("{sample}/join_reads/{kmer}/{region}/{map}-{map2}/group_{split}.txt".format(sample=wildcards.sample, map=TYPE_MAP, map2=REMAP, region=wildcards.region, kmer=wildcards.kmer, split='{split}')).split
	return expand(rules.gather_ONT_per_split.output.flag, sample=wildcards.sample, map=TYPE_MAP, map2=REMAP, region=wildcards.region, kmer=wildcards.kmer, split=SPLIT)

localrules: all, sunk_tag_ONT, sunk_tag_HiFi, get_overlapping_reads

###################### SAVE HIFI AND ONT READS ##########################

rule all:
	input: 
		expand('{sample}/join_reads/{kmer}/{region}/{map}-{map2}/joining_reads.paf', sample=SAMPLE, kmer=['%d' % size for size in KMER_LIST], region=['%s' % size for size in CONTIG_LIST], map=TYPE_MAP, map2=REMAP),

rule save_HiFi:
	input:
		reads = findHiFiReads
	output: 
		temp_read = "{sample}/reads/tmp/temp_HiFi_{read}.read",
	resources:
		mem = 8,
		disk = 0,
		hrs = 12
	threads: 8
	run:
		if str(input.reads)[-3:] == ".gz":
			shell_string = 'gunzip -c %s > {output.temp_read}' % (input.reads)
			shell(shell_string)
		else:
			shell_string = 'rsync -av --bwlimit=10000 $(readlink -f {input.reads}) {output.temp_read}'
			shell(shell_string)

checkpoint save_split_ONT:
	input:
		read = findONTReads
	output:
		split = "{sample}/reads/tmp/ONT/temp_ONT_{read}_0000",
		flag = touch("{sample}/reads/tmp/.checkpoint_ONT_{read}")
	priority: 50 
	resources: 
		mem = 8,
		disk = 100,
		hrs = 12
	threads: 8 
	params:
		directory = "{sample}/reads/tmp/ONT/"
	run:

		if str(input.read)[-3:] == ".gz":
			shell_string = 'gunzip -c %s > $TMPDIR/temp_ONT_{wildcards.read}.read' % (input.read)
			shell(shell_string)
		else:
			shell_string = 'rsync -av --bwlimit=10000 $(readlink -f {input.read}) $TMPDIR/temp_ONT_{wildcards.read}.read'
			shell(shell_string)

		shell('/bin/ls $TMPDIR/temp_ONT_{wildcards.read}.read | parallel -j{threads} \'split -a 4 -d -l 800000 {{}} {params.directory}{{/.}}_\'')

####################### BEGIN ONT ONLY RULES ##############################
rule ont_fastq_to_fasta:
	input:
		ont_fastq = "{sample}/reads/tmp/ONT/temp_ONT_{read}_{part}"
	output:
		ont_fasta = temp("{sample}/fasta/tmp/ONT/ONT_{read}_{part}.fasta")
	threads: 1
	resources:
		mem = 18,
		disk = 0,
		hrs = 8
	shell:
		'''
		seqtk seq -A -U -l 60 -L {LEN_FILT}000 {input.ont_fastq} > {output.ont_fasta}
		'''
rule combine_ONT_per_read:
	input:
		reads = gatherPartsONT,
		flag = "{sample}/reads/tmp/.checkpoint_ONT_{read}"
	output:
		combined = temp('{sample}/fasta/ONT/{read}.fasta'),
		flag_2 = touch("{sample}/reads/tmp/.checkpoint_ONT_{read}.done")
	threads: 1
	resources:
		mem = 8,
		disk = 0,
		hrs = 2
	shell:
		'''
		cat {input.reads} > {output.combined}
		'''

rule combine_ONT_reads_all:
	input:
		reads = expand(rules.combine_ONT_per_read.output.combined, sample=SAMPLE, read=ONT_FOFN_DF.index)
	output:
		all_reads = '{sample}/fasta/ONT/all.fasta',
		fai = '{sample}/fasta/ONT/all.fasta.fai'
	threads: 1
	resources:
		mem=12,
		disk =0,
		hrs = 18
	shell:
		'''
		cat {input.reads} > {output.all_reads}
		samtools faidx {output.all_reads}
		'''
rule countKmers:
	input:
		assembly = HIFI_ASSEMBLY
	output:
		dirct = directory("{sample}/alignments/merylDB_{sample}")
	threads: 1
	resources:
		mem = 12,
		disk = 0,
		hrs = 8
	shell:
		'''
		meryl count k=15 output {output.dirct} {input.assembly}
		'''
rule getRepeatKmers:
	input:
		db = rules.countKmers.output.dirct
	output:
		rep = "{sample}/alignments/repetitive_k15_{sample}.txt"
	threads: 1
	resources:
		mem = 8,
		disk = 0,
		hrs = 8
	shell:
		'''
		meryl print greater-than distinct=0.9998 {input.db} > {output.rep}
		'''

rule winnow_map_align:
	input:
		ont_fasta = '{sample}/fasta/tmp/ONT/ONT_{read}_{part}.fasta',
		assembly = HIFI_ASSEMBLY,
		repKmers = rules.getRepeatKmers.output.rep
	output:
		bam = temp('{sample}/alignments/ONT/winnowmap/{read}_{part}.bam')
	threads: 8
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	shell:
		'''
		winnowmap -W {input.repKmers} -t {threads} -I 7G -a -x map-ont {input.assembly} {input.ont_fasta} | samtools sort -o {output.bam} -
		'''

rule minimap_align:
	input:
		ont_fasta = '{sample}/fasta/tmp/ONT/ONT_{read}_{part}.fasta',
		assembly = HIFI_ASSEMBLY
	output:
		bam = temp('{sample}/alignments/ONT/minimap2/{read}_{part}.bam')
	threads: 8
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	shell:
		'''
		minimap2 -a -t {threads} -I 7G -x map-ont {input.assembly} {input.ont_fasta} | samtools sort -o {output.bam} -
		'''


rule combine_ONT_alignments_per_read:
	input:
		align = findONTAlign,
		flags = "{sample}/reads/tmp/.checkpoint_ONT_{read}"
	output:
		combined = temp('{sample}/alignments/ONT/{map}/{read}_all.bam'),
		flag_4 = touch("{sample}/alignments/tmp/.checkpoint_ONT_{map}_{read}.done")
	threads: 4
	resources:
		mem = 8,
		disk = 0,
		hrs = 12
	shell:
		'''
		samtools merge -@ {threads} {output.combined} {input.align}
		'''
rule combine_ONT_alignments_all:
	input:
		align_minimap = expand('{sample}/alignments/ONT/minimap2/{read}_all.bam', sample=SAMPLE, read=ONT_FOFN_DF.index),
		#align_winnowmap = expand('{sample}/alignments/ONT/winnowmap/{read}_all.bam', sample=SAMPLE, read=ONT_FOFN_DF.index)
	output:
		combined_minimap = '{sample}/alignments/ONT/minimap2/all_ont.bam',
		#combined_winnowmap = '{sample}/alignments/ONT/winnowmap/all_ont.bam'
	threads: 4
	resources:
		mem = 8,
		disk = 0,
		hrs = 12
	shell:
		'''
		samtools merge -@ {threads} {output.combined_minimap} {input.align_minimap}
		samtools index {output.combined_minimap}
		'''

rule index_ONT_bam:
	input:
		bam = '{sample}/alignments/ONT/{map}/all_ont.bam'
	output:
		bai = '{sample}/alignments/ONT/{map}/all_ont.bam.bai'
	threads: 1
	resources:
		mem = 24,
		disk = 0,
		hrs = 12
	shell:
		'''
		samtools index {input.bam}
		'''
rule extract_reads:
	input:
		bam = '{sample}/alignments/ONT/{map}/all_ont.bam',
		bai = '{sample}/alignments/ONT/{map}/all_ont.bam.bai'
	output:
		align_filt = '{sample}/alignments/ONT/{map}/{region}/filtered.sam',
		read_list = "{sample}/alignments/ONT/{map}/{region}/read_list.bed"
	threads: 1
	resources:
		mem = 12,
		disk = 0,
		hrs = 24
	params:
		contigs = findRegion
	shell:
		'''
		samtools view {input.bam} {params.contigs} > {output.align_filt}
		cut -f1 {output.align_filt} | sort -u > {output.read_list}
		'''

rule ont_fai_to_bed:
	input: 
		fai =  '{sample}/fasta/ONT/all.fasta.fai',
		read_list = rules.extract_reads.output.read_list
	output:
		bed = temp("{sample}/alignments/ONT/{map}/{region}/reads.bed")
	threads: 1
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	run:
		fai_df = pd.read_csv(input.fai, header=None, names=['contig', 'len', 'offset', 'byte', 'encode'],sep="\t")
		read_df = pd.read_csv(input.read_list, header=None, names=['read_name'], sep="\t")
		read_list=read_df['read_name'].tolist()
                print(read_list[0])
		data = {'chrom':[],'start':[],'stop':[]}
		bed_df =  pd.DataFrame(data)
		bed_df['chrom'] = read_list
		bed_df = bed_df.assign(start='0')
		for read in read_list:
			bed_df.loc[bed_df['chrom'] == read, 'stop'] = str(fai_df.loc[fai_df['contig'] == read, 'len'].iloc[0])
		bed_df.to_csv(output.bed, header=None, index=None, sep='\t')

rule get_filtered_ONT_fasta:
	input:
		fasta = '{sample}/fasta/ONT/all.fasta',
		read_bed = "{sample}/alignments/ONT/{map}/{region}/reads.bed"
	output:
		filtered_fasta = '{sample}/fasta/ONT/{map}/{region}/filtered_fasta.fa',
		filtered_fasta_fai = '{sample}/fasta/ONT/{map}/{region}/filtered_fasta.fa.fai'
	threads: 1
	resources:
		mem = 64,
		disk = 0,
		hrs = 12
	shell:
		'''
		bedtools getfasta -fi {input.fasta} -bed {input.read_bed} -fo {output.filtered_fasta}
		samtools faidx {output.filtered_fasta}
		'''
##################### END OF ONT ONLY RULES ######################################
##################### BEGIN HIFI ONLY RULES ######################################

rule hifi_fastq_to_fasta:
	input: 
		hifi_fastq = "{sample}/reads/tmp/temp_HiFi_{read}.read",
	output:
		hifi_fasta = "{sample}/fasta/HiFi/{read}.fasta.gz"
	threads: 1
	resources:
		mem = 86,
		disk = 0,
		hrs = 8
	shell:
		'''
		seqtk seq -A -U -l 60 {input.hifi_fastq} | bgzip -c > {output.hifi_fasta}
		'''

rule map_hifi_to_hifi:
	input:
		hifi_fasta = "{sample}/fasta/HiFi/{read}.fasta.gz",
		assembly = HIFI_ASSEMBLY
	output:
		bam = temp("{sample}/alignments/HiFi/{read}.bam")
	threads: 8
	resources:
		mem = 12,
		disk = 0,
		hrs = 8
	shell:
		'''
		pbmm2 align --log-level DEBUG --preset CCS --min-length 5000 -j {threads} {input.assembly} {input.hifi_fasta} | samtools sort -@ {threads} -o {output.bam} -
		'''

rule combine_HiFi_alignments_all:
	input:
		align = expand("{sample}/alignments/HiFi/{read}.bam", sample=SAMPLE, read=HIFI_FOFN_DF.index)
	output:
		combined = '{sample}/alignments/HiFi/all_HiFi.bam'
	threads: 4
	resources:
		mem = 8,
		disk = 0,
		hrs = 12
	shell:
		'''
		samtools merge -@ {threads} {output.combined} {input.align}
		samtools index {output.combined}
		'''

rule hifi_fai_to_bed:
	input: 
		fai = HIFI_FAI,
	output:
		bed = "{sample}/break_contigs/{region}/contigs.bed"
	threads: 1
	resources: 
		mem = 12,
		disk = 0,
		hrs = 12
	params:
		contigs = findRegion 
	run:
		fai_df = pd.read_csv(input.fai, header=None, names=['contig', 'len', 'offset', 'byte', 'encode'],sep="\t")
		contig_list = params.contigs.split(" ")
		data = {'chrom':[],'start':[],'stop':[]}
		bed_df =  pd.DataFrame(data)
		bed_df['chrom'] = contig_list
		bed_df = bed_df.assign(start='0')
		for contig in contig_list:
			bed_df.loc[bed_df['chrom'] == contig, 'stop'] = str(fai_df.loc[fai_df['contig'] == contig, 'len'].iloc[0])

		bed_df.to_csv(output.bed, header=None, index=None, sep='\t')

checkpoint get_depth_cov:
	input: 
		bed = rules.hifi_fai_to_bed.output.bed,
		bam = rules.combine_HiFi_alignments_all.output.combined
	output:
		depth = "{sample}/break_contigs/{region}/filtered.bed"
	threads: 1
	resources:
		mem = 24,
		disk = 0,
		hrs = 12
	shell:
		'''
		samtools depth -b {input.bed} -a {input.bam} | awk '{{if ($3 > 10+{COV} || $3==0) printf ("%s\\t%s\\t%s\\n", $1, $2, $2)}}' | bedtools merge -i - > {output.depth} 
		'''
rule nucFreak:
	input:
		depth = rules.get_depth_cov.output.depth,
		bam = rules.combine_HiFi_alignments_all.output.combined
	output:
		bed = "{sample}/break_contigs/{region}/intermediate.bed"
	threads: 4
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	shell:	
		'''
		rustybam nucfreq -t {threads} --bed {input.depth} {input.bam} > {output.bed}                
		'''
rule filter_nucFreak:
	input: 
		bed = rules.nucFreak.output.bed
	output:
		filtered_bed = "{sample}/break_contigs/{region}/intermediate_2.bed"
	threads: 1
	resources:
		mem = 20,
		disk = 0,
		hrs = 12
	run:
		bed_df = pd.read_csv(input.bed, sep="\t", header=0)
		bed_df['second_highest'] = bed_df.iloc[:, ][['A','G','T','C']].apply(lambda row: row.nlargest(2).values[-1],axis=1)
                collapse_df = bed_df.loc[bed_df['second_highest'] > 5]
                collapse_df = collapse_df[['#chr','start','end']]
                missassembly_df = bed_df.loc[(bed_df['A'] == 0) & (bed_df['C'] == 0 ) & (bed_df['G'] == 0) & (bed_df['T'] == 0)]
                missassembly_df = missassembly_df[['#chr','start','end']]
                combined_df = pd.concat([missassembly_df, collapse_df], ignore_index=True)
                combined_df = combined_df.sort_values(by=['#chr', 'start','end'])
                combined_df.to_csv(output.filtered_bed, sep="\t", header=None, index=None)

rule filter_bed:
	input:
		bed = rules.filter_nucFreak.output.filtered_bed
	output:
		filtered_bed = "{sample}/break_contigs/{region}/intermediate_3.bed"
	threads: 1
	resources:
		mem = 20,
		disk = 0,
		hrs = 12
	shell:
		'''
		bedtools merge -i {input.bed} -c 1 -o count -d 5000 | awk '$4 > 2' | bedtools merge -i - -d 15000 | awk '{{printf "%s\\t%s\\t%s\\n", $1, $2-5000, $3+5000}}' > {output.filtered_bed}
		'''

rule find_breaks:
	input:
		bed = rules.filter_bed.output.filtered_bed
	output:
		break_bed = "{sample}/break_contigs/{region}/break.bed"
	threads: 1
        resources:
		mem = 8,
		disk = 0,
		hrs = 1
	run:
		bed_df=pd.read_csv(input.bed, sep="\t", header=None, names=["chr","start","stop"])
		fai_df = pd.read_csv(HIFI_FAI, sep="\t", header=None, names=['contig', 'len', 'offset', 'byte', 'encode'])
		clear_df = bed_df

		for ind in bed_df.index:
			is_cleared=False
			if int(bed_df.at[ind,'start']) < 25000:
				print(clear_df)
				clear_df=clear_df.drop([ind])
				is_cleared=True
			contig = bed_df.loc[ind, 'chr']
			fai_len = int(fai_df.loc[fai_df['contig'] == contig, 'len'])
			bed_stop = int(bed_df.loc[ind,'stop'])
			if bed_stop > fai_len-25000 and not is_cleared:
				print(clear_df)
				clear_df=clear_df.drop([ind])
		clear_df.to_csv(output.break_bed, sep="\t", header=None, index=None)

 
rule bed_subtract:
	input:
		break_bed = rules.find_breaks.output.break_bed,
		contigs_bed = rules.hifi_fai_to_bed.output.bed
	output:
		final_bed = '{sample}/break_contigs/{region}/final.bed'
	threads: 1
	resources:
		mem = 8,
		disk = 0,
		hrs = 1
	shell:
		'''
		bedtools merge -i {input.break_bed} | bedtools subtract -a {input.contigs_bed} -b - > {output.final_bed}
		'''

rule break_contigs:
	input:
		hifi_assembly = HIFI_ASSEMBLY,
		final_bed = aggregate_input
	output:
		broken_contigs = "{sample}/break_contigs/{region}/hifi_assembly_broken.fasta",
		broken_contigs_fai = "{sample}/break_contigs/{region}/hifi_assembly_broken.fasta.fai"
	threads: 1
	resources:
		mem = 24,
		disk = 0,
		hrs = 12
	shell:
		'''
		bedtools getfasta -fi {input.hifi_assembly} -bed {input.final_bed} -fo {output.broken_contigs}
		samtools faidx {output.broken_contigs}
		'''
rule reverse_hifi:
	input:
		bed = getContig,
		fasta = rules.break_contigs.output.broken_contigs,
		fai = rules.break_contigs.output.broken_contigs_fai
	output:
		reverse_comp = '{sample}/break_contigs/{region}/reverse_comp_hifi.fasta'
	resources: 
		mem = 25,
		disk = 0,
		hrs = 12,
	threads: 1
	run:
		bed_df = pd.read_csv(input.bed, header=None, sep="\t", names=['contig', 'start','stop','hifi_contig', 'mapq','strand'])
		fai_df = pd.read_csv(input.fai, header=None, sep="\t", names=['hifi_contig','len','offset','byte','encode'])
		hifi_records = SeqIO.to_dict(SeqIO.parse(input.fasta, "fasta"))

		hifi_contigs=fai_df['hifi_contig'].unique()

		outfile = open(output.reverse_comp, "w+")
		for hifi_contig in hifi_contigs:
			just_contig = hifi_contig.split(":")[0]
			hifi_record = hifi_records[hifi_contig]
			hifi_seq = hifi_record.seq
			contig_df = bed_df.loc[bed_df['hifi_contig'] == just_contig]
			if contig_df['strand'].iloc[0] == "-":
				write_seq = hifi_seq.reverse_complement()
			else:
				write_seq = hifi_seq
			header = hifi_record.description
			outfile.write(">"+str(header)+"\n")
			outfile.write(str(write_seq)+"\n")

		outfile.close()

rule jellyfish_count:
	input:
		fasta = '{sample}/fasta/HiFi/{read}.fasta.gz'
	output:
		tally = temp('{sample}/jellyfish/HiFi/count/{kmer}/{read}.output')
	threads: 16
	resources: 
		mem = 6,
		disk = 50,
		hrs = 24
	shell:
		'''
		zcat {input.fasta} > $TMPDIR/{wildcards.read}.fa
		jellyfish count -m {wildcards.kmer} -C -o $TMPDIR/{wildcards.read}.output -c 3 -s 10000000 --disk -t {threads} $TMPDIR/{wildcards.read}.fa
		rsync -av --bwlimit=10000 $TMPDIR/{wildcards.read}.output {output.tally}
		'''
rule jellyfish_merge:
	input: 
		tally = jellyCounts
	output:
		combined = temp('{sample}/jellyfish/HiFi/count/{kmer}/jellyfish.combined')
	threads: 1
	resources: 
		mem = 16,
		disk = 0,
		hrs = 4
	shell:
		'''
		jellyfish merge -o {output.combined} {input.tally}
		'''

rule jellyfish_hist:
	input:
		tally = rules.jellyfish_merge.output.combined
	output:
		histo = '{sample}/jellyfish/HiFi/hist/{kmer}/jellyfish.hist'
	threads: 4
	resources:
		mem = 8,
		disk = 0,
		hrs = 2
	shell:
		'''
		jellyfish histo -l 1 -h 60 -i 1 -t {threads} -o {output.histo} -v {input.tally}
		'''
rule jellyfish_dump:
	input:
		histo = rules.jellyfish_hist.output.histo,
		counts = rules.jellyfish_merge.output.combined
	output:
		dump = temp('{sample}/jellyfish/HiFi/dump/{kmer}/jellyfish.dump'),
		db = '{sample}/sunks/{kmer}/HiFi/jellyfish.db',
		hist = '{sample}/jellyfish/HiFi/hist/{kmer}/jellyfish.hist.png',
		bounds = temp('{sample}/jellyfish/HiFi/dump/{kmer}/jellyfish.dump.tsv')
	threads: 1
	resources:
		mem = 24,
		disk = 0,
		hrs = 6
	run:
		
		import matplotlib.pyplot as plt
		count_df = pd.read_csv(input.histo, sep='\s+', header=None, names=['freq', 'count'])
		count_df['weights'] = count_df['freq'] * count_df['count']
		mean = int(COV)
		search_lower = np.floor(mean/6)
		search_upper = np.ceil(mean+5*search_lower)
		upper_cut = count_df[count_df['count'] == min(count_df[(count_df['freq'] <= search_upper) & (count_df['freq'] >= mean)]['count'])]['freq']
		lower_cut = count_df[count_df['count'] == min(count_df[(count_df['freq'] <= mean) & (count_df['freq'] >= search_lower)]['count'])]['freq']
		search_df = count_df[(count_df['freq'] <= max(upper_cut)) & (count_df['freq'] >= min(lower_cut))]
		sigma = 0
		for index in search_df.index:
			sigma += search_df.at[index,'count']*(search_df.at[index,'freq']-mean)**2
		st_dev = np.sqrt(sigma/np.sum(search_df['count']))
		upper_lim = int(np.ceil(mean+UPPER_LIM*st_dev))
		lower_lim = int(max(0, np.floor(mean-LOWER_LIM*st_dev)))
		kmer_out = count_df[(count_df['freq'] <= upper_lim) & (count_df['freq'] >= lower_lim)]
		shell('jellyfish dump -c -t -L %d -U %d {input.counts} > {output.dump}' % (lower_lim, upper_lim))
		shell('awk \'{{print $1}}\' {output.dump} > {output.db}')

		with open(output.bounds, 'w') as outFile:
			outFile.write('Mean: %f\nStd: %f\nUpper: %d\nLower: %d\nkmers: %s' % (mean, st_dev, upper_lim, lower_lim, "{:,}".format(np.sum(kmer_out['count']))))

		fig, ax = plt.subplots()
		rects = ax.bar(count_df['freq'], count_df['count'], 1)
		ax.axvline(x=upper_lim, ymin=0, ymax=0.5, color='firebrick')
		ax.axvline(x=lower_lim, ymin=0, ymax=0.5, color='firebrick')
		plt.ylabel('Count')
		plt.title('Kmer Distribution')
		plt.xlabel('Occurence')
		plt.tight_layout()
		plt.savefig(output.hist)



######################### END OF HIFI ONLY RULES ############################
#############################################################################

rule remap_countKmers:
        input:
                ref = rules.reverse_hifi.output.reverse_comp
        output:
                dirct = directory('{sample}/alignments/remap/{map}-winnowmap/{region}/merlyDB_{sample}')
        threads: 1
        resources:
                mem = 12,
                disk = 0,
                hrs = 8
        shell:
                '''
                meryl count k=15 output {output.dirct} {input.ref}
                '''
rule remap_getRepeatKmers:
	input:
		db = rules.remap_countKmers.output.dirct
	output:
		rep = '{sample}/alignments/remap/{map}-winnowmap/{region}/repetitive_k15_{sample}.txt'
	threads: 1
	resources:
		mem = 8,
		disk = 0,
		hrs = 8
	shell:
		'''
		meryl print greater-than distinct=0.9998 {input.db} > {output.rep}
		'''
rule remap_winnow_map_align:
	input:
		ont_fasta = rules.get_filtered_ONT_fasta.output.filtered_fasta,
		assembly = rules.reverse_hifi.output.reverse_comp,
		repKmers = rules.remap_getRepeatKmers.output.rep
	output:
		paf = '{sample}/alignments/remap/{map}-winnowmap/{region}/ont_to_hifi.paf'
	threads: 8
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	shell:
		'''
		winnowmap -W {input.repKmers} -t {threads} -x map-ont {input.assembly} {input.ont_fasta} > {output.paf}
		'''
rule remap_minimap2:
	input: 
		ont_fasta = rules.get_filtered_ONT_fasta.output.filtered_fasta,
		assembly = rules.reverse_hifi.output.reverse_comp,
	output:
		paf = '{sample}/alignments/remap/{map}-minimap2/{region}/ont_to_hifi.paf'
	threads: 8
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	shell:
		'''
		minimap2 -t {threads} -x map-ont {input.assembly} {input.ont_fasta} > {output.paf} 
		'''
rule filter_paf:
	input:
		paf = '{sample}/alignments/remap/{map}-{map2}/{region}/ont_to_hifi.paf'
	output:
		filtered_paf = '{sample}/alignments/remap/{map}-{map2}/{region}/ont_to_hifi_filt.paf',
		read_list = temp('{sample}/alignments/remap/{map}-{map2}/{region}/read_list.txt')
	threads: 1
	resources:
		mem = 64,
		disk = 0,
		hrs = 12
	run:
		paf_df = pd.read_csv(input.paf, header=None, names=['query_name', 'query_len', 'query_start', 'query_end', 'strand', 'target_name', 'target_len', 'target_start', 'target_end', 'residual_matches', 'block_len', 'map_q', 'tp', 'cm', 's1', 's2', '#rl'], sep='\t')

		filtered_paf_df =  paf_df.query('block_len < query_len')

		reads = (filtered_paf_df["query_name"].unique()).tolist()

		for read in reads:
			read_df = filtered_paf_df.loc[filtered_paf_df['query_name'] == read]
			read_contigs_list = read_df["target_name"].unique().tolist()
			if len(read_contigs_list) < 2:
				filtered_paf_df = filtered_paf_df.loc[filtered_paf_df['query_name'] != read]
		reads = (filtered_paf_df["query_name"].unique())
		reads = pd.DataFrame(reads, columns=['query_name'])
		reads.to_csv(output.read_list, sep="\t", header=None, index=None)
		filtered_paf_df.to_csv(output.filtered_paf, sep="\t", header=None, index=None)

rule paf_fai_to_bed:
	input: 
		fai = '{sample}/fasta/ONT/{map}/{region}/filtered_fasta.fa.fai',
		read_list = rules.filter_paf.output.read_list
	output:
		bed = temp('{sample}/alignments/remap/{map}-{map2}/{region}/reads.bed')
	threads: 1
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	run:
		fai_df = pd.read_csv(input.fai, header=None, names=['contig', 'len', 'offset', 'byte', 'encode'],sep="\t")
		read_df = pd.read_csv(input.read_list, header=None, names=['read_name'], sep="\t")
		read_list=read_df['read_name'].tolist()
		data = {'chrom':[],'start':[],'stop':[]}
		bed_df =  pd.DataFrame(data)
		bed_df['chrom'] = read_list
		bed_df = bed_df.assign(start='0')
		for read in read_list:
			bed_df.loc[bed_df['chrom'] == read, 'stop'] = str(fai_df.loc[fai_df['contig'] == read, 'len'].iloc[0])
		bed_df.to_csv(output.bed, header=None, index=None, sep='\t')

rule paf_filtered_fasta:
	input:
		fasta = '{sample}/fasta/ONT/{map}/{region}/filtered_fasta.fa',
		read_bed = '{sample}/alignments/remap/{map}-{map2}/{region}/reads.bed'
	output:
		filtered_fasta = '{sample}/fasta/ONT/{map}-{map2}/{region}/filtered_fasta.fa'
	threads: 1
	resources:
		mem = 64,
		disk = 0,
		hrs = 12
	shell:
		'''
		bedtools getfasta -fi {input.fasta} -bed {input.read_bed} -fo {output.filtered_fasta}
		samtools faidx {output.filtered_fasta}
		'''

rule sunk_tag_ONT:
	input:
		ont_fasta = rules.paf_filtered_fasta.output.filtered_fasta,
		sunks = '{sample}/sunks/{kmer}/HiFi/jellyfish.db'
	output:
		sunk_pos = '{sample}/sunks/{kmer}/ONT/{map}-{map2}/{region}/ONT.sunks_pos.tab'
	threads: 1
	resources: 
		mem = 250,
		disk = 0, 
		hrs = 24
	shell:
		'''
		{SNAKEMAKE_DIR}/scripts/KmersPos2.py -k {input.sunks} -fa {input.ont_fasta} -p {output.sunk_pos}
		'''

rule sunk_tag_HiFi:
	input:
		hifi_fasta = rules.reverse_hifi.output.reverse_comp,
		sunks = '{sample}/sunks/{kmer}/HiFi/jellyfish.db'
	output:
		sunk_pos = '{sample}/sunks/{kmer}/HiFi/{region}/HiFi.sunks_pos.tab'
	threads: 1
	resources: 
		mem = 250,
		disk = 0, 
		hrs = 24
	shell:
		'''
		{SNAKEMAKE_DIR}/scripts/KmersPos2.py -k {input.sunks} -fa {input.hifi_fasta} -p {output.sunk_pos}
		'''
checkpoint get_overlapping_reads:
	input:
		sunk_reads = '{sample}/sunks/{kmer}/ONT/{map}-{map2}/{region}/ONT.sunks_pos.tab',
		sunk_asm = '{sample}/sunks/{kmer}/HiFi/{region}/HiFi.sunks_pos.tab',
		paf = '{sample}/alignments/remap/{map}-{map2}/{region}/ont_to_hifi_filt.paf'
	output:
		filt_paf = '{sample}/join_reads/{kmer}/{region}/{map}-{map2}/joining_reads.paf',
		read_list = '{sample}/join_reads/{kmer}/{region}/{map}-{map2}/group_0.txt',
	threads: 1
	params:
		dirct = "{sample}/join_reads/{kmer}/{region}/{map}-{map2}/"
	resources:
		mem = 60,
		disk = 0,
		hrs = 48
	run:
		sunk_reads = pd.read_csv(input.sunk_reads, sep='\t', header=None, names=['contig', 'sunk', 'pos'])
		sunk_asm = pd.read_csv(input.sunk_asm, sep='\t', names=['contig', 'sunk', 'pos'], header=None)
		paf_filt = pd.read_csv(input.paf, header=None, names=['query_name', 'query_len', 'query_start', 'query_end', 'strand', 'target_name', 'target_len', 'target_start', 'target_end', 'residual_matches', 'block_len', 'map_q', 'tp', 'cm', 's1', 's2', 'rl'], sep='\t')

		sunk_reads['contig']=sunk_reads['contig'].str.split(pat=":").str[:2].str.join(":")
		
		for index in paf_filt.index:
			asm_sub = sunk_asm.loc[((sunk_asm['pos'] >= paf_filt.at[index, 'target_start']) & (sunk_asm['pos'] <= paf_filt.at[index, 'target_end']) & (sunk_asm['contig'] == paf_filt.at[index, 'target_name']))]
			read_sub = sunk_reads.loc[((sunk_reads['pos'] >= paf_filt.at[index, 'query_start']) & (sunk_reads['pos'] <= paf_filt.at[index, 'query_end']) & (sunk_reads['contig'] == paf_filt.at[index, 'query_name']))]
			read_check = read_sub.groupby('sunk').size().reset_index().rename({0: 'occurrence'}, axis=1)
			read_sub = pd.merge(read_sub, read_check.loc[read_check['occurrence'] == 1][['sunk']])
			merged = pd.merge(asm_sub, read_sub, on='sunk', suffixes=['_asm', '_reads'])
			#merged_sorted = merged.sort_values(by=['pos_asm'])
			#for ind in merged_sorted.index:
			#	if ind < len(merged_sorted) -1:
			#		target_dis = abs(int(merged_sorted.at[ind+1, 'pos_asm']) - int(merged_sorted.at[ind, 'pos_asm']))
			#		query_dis = abs(int(merged_sorted.at[ind+1, 'pos_reads']) - int(merged_sorted.at[ind, 'pos_reads']))
			#		if target_dis != query_dis:
			#			merged = merged.loc[merged['sunk'] != merged_sorted.at[ind,'sunk']]
			#			record = merged_sorted.iloc[ind]
			#			record['dis_asm'] = target_dis
			#			record['dis_reads'] = query_dis
						
			paf_filt.at[index, 'sunk_intersect'] = len(merged['sunk'].unique())

		paf_join = paf_filt.loc[(paf_filt['sunk_intersect'] > 2) & ((paf_filt['block_len'] >= 50000) | (paf_filt['block_len']/paf_filt['target_len'] > 0.5))]

		paf_join.groupby('query_name')

		paf_group = paf_join.groupby(['query_name', 'target_name']).size().reset_index()

		paf_group_filt = pd.DataFrame()

		for read_name in paf_group['query_name'].unique():
			if len(paf_group.loc[paf_group['query_name'] == read_name]) > 1:
				paf_group_filt = paf_group_filt.append(pd.DataFrame.from_dict({'query_name' : [read_name]}))

		paf_join = pd.merge(paf_join, paf_group_filt)
		paf_join.to_csv(output.filt_paf, sep="\t", header=None, index=None)
		reads = paf_join['query_name'].unique()
		matching_df = pd.DataFrame(columns=['reads','target'])
		for read in reads:
			read_df = paf_join.loc[paf_join['query_name']==read]
			target_list = (read_df['target_name'].unique()).tolist()
			target_list.sort()
			target_string = ",".join(target_list)
			matching_df = matching_df.append({'read': read, 'target': target_string}, ignore_index=True)

		groups = matching_df.groupby(['target'])
		count = 0
		for (splitno, split) in groups:
			split[['read']].to_csv(params.dirct+"group_"+str(count)+".txt", header=None, index=None, sep="\t")
			count = count+1

#################################END OF FUNCTIONAL PIPELINE############################################################

rule paf_final_to_bed:
	input:
		fai = '{sample}/fasta/ONT/{map}/{region}/filtered_fasta.fa.fai',
		read_list = '{sample}/join_reads/{kmer}/{region}/{map}-{map2}/group_{split}.txt'
	output:
		bed = '{sample}/alignments/remap/{map}-{map2}/{region}/{kmer}/reads_final_{split}.bed'
	threads: 1
	resources:
		mem = 12,
		disk = 0,
		hrs = 12
	run:
		fai_df = pd.read_csv(input.fai, header=None, names=['contig', 'len', 'offset', 'byte', 'encode'],sep="\t")
		read_df = pd.read_csv(input.read_list, header=None, names=['read_name'], sep="\t")
		read_list=read_df['read_name'].tolist()
		data = {'chrom':[],'start':[],'stop':[]}
		bed_df =  pd.DataFrame(data)
		bed_df['chrom'] = read_list
		bed_df = bed_df.assign(start='0')
		for read in read_list:
			bed_df.loc[bed_df['chrom'] == read, 'stop'] = str(fai_df.loc[fai_df['contig'] == read, 'len'].iloc[0])
		bed_df.to_csv(output.bed, header=None, index=None, sep='\t')

rule paf_final_fasta:
	input:
		fasta = '{sample}/fasta/ONT/{map}/{region}/filtered_fasta.fa',
		read_bed = rules.paf_final_to_bed.output.bed
	output:
		filtered_fasta = '{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/filtered_fasta_final_{split}.fa',
		fai = '{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/filtered_fasta_final_{split}.fa.fai'
	threads: 1
	resources:
		mem = 64,
		disk = 0,
		hrs = 12
	shell:
		'''
		bedtools getfasta -fi {input.fasta} -bed {input.read_bed} -fo {output.filtered_fasta}
		samtools faidx {output.filtered_fasta}
		'''

checkpoint generate_medaka:
	input:
		fasta = '{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/filtered_fasta_final_{split}.fa',
		fai = '{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/filtered_fasta_final_{split}.fa.fai',
	output:
		flag = touch("{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/.medaka.done"),
		correction = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/0/correction.fa",
		ref = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/0/ref.fa"
	threads:1
	params:
		dirct = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/"
	resources:
		mem = 8,
		hrs = 12,
		disk = 0
	shell:
		'''
		count=0; for file in $(cut -f1 {input.fai}); do mkdir -p {params.dirct}${{count}}; samtools faidx {input.fasta} ${{file}} > {params.dirct}${{count}}/ref.fa; grep -v ${{file}} {input.fai} | cut -f 1 | xargs -i samtools faidx {input.fasta} {{}} > {params.dirct}${{count}}/correction.fa; samtools faidx {params.dirct}${{count}}/ref.fa; samtools faidx {params.dirct}${{count}}/correction.fa; count=$((count+1)); done
		'''

rule medaka_consensus:
	input:
		correction = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/{ontread}/correction.fa",
		ref = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/{ontread}/ref.fa"
	output:
		 dirct = directory("{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/{ontread}/medaka_consensus/"),
		 consensus = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/{ontread}/medaka_consensus/consensus.fasta",
	threads: 16
	resources:
		mem = 6,
		hrs = 12,
		disk = 0
	shell:
		'''
		module load bcftools/1.12 samtools/1.12 htslib/1.12 medaka/1.2.6
		medaka_consensus -i {input.correction} -d {input.ref} -o {output.dirct} -t {threads} -m r941_min_high_g360
		'''
rule string_decomposer_ONT:
	input:
		monomers=getMonomer,
		consensus = rules.medaka_consensus.output.consensus
	output:
		tsv = "{sample}/string_decomposer/ONT/{map}-{map2}/{region}/{kmer}/{split}/{ontread}/string.tsv"
	threads: 8
	resources:
		mem = 12,
		hrs = 12,
		disk = 0
	shell:
		'''
		run_decomposer.py {input.consensus} {input.monomers} -t {threads} -o {output.tsv}
		'''

rule string_decomposer_bed:
	input:
		paf = '{sample}/join_reads/{kmer}/{region}/{map}-{map2}/joining_reads.paf'

	output:
		bed = "{sample}/string_decomposer/bed/{region}/{kmer}/{map}-{map2}/filtered-region.bed"
	threads: 1
	resources:
		mem = 12,
		hrs = 6,
		disk = 0
	run:
		paf_df = pd.read_csv(input.paf, header=None, names=['query_name', 'query_len', 'query_start', 'query_end', 'strand', 'target_name', 'target_len', 'target_start', 'target_end', 'residual_matches', 'block_len', 'map_q', 'tp', 'cm', 's1', 's2', 'rl', 'num'], sep='\t')
		bed_df = pd.DataFrame(columns=['contig','start','stop'])

		contig_list = (paf_df['target_name'].unique()).tolist()
		region_list = []

		for contig in contig_list:
			read_df = paf_df.loc[paf_df['target_name'] == contig]
			read_df.reset_index(inplace=True)
			least_start = int(read_df.at[0,'target_len'])
			greatest_end = 0
			for ind in read_df.index:
				start_dis = int(paf_df.at[ind, 'target_start'])
				end_dis = int(paf_df.at[ind, 'target_end']) 

				if start_dis < least_start:
					least_start = start_dis
				if end_dis > greatest_end:
					 greatest_end = end_dis


			least_start = least_start - 5000
			if least_start < 0:
				least_start = 0

			greatest_end = greatest_end + 5000
			if greatest_end > int(read_df.at[0,'target_len']):
				greatest_end = int(read_df.at[0,'target_len'])

			bed_df = bed_df.append({'contig' : contig, 'start': str(least_start), 'stop' : str(greatest_end)}, ignore_index=True)
		bed_df.to_csv(output.bed, header=None, index=None, sep="\t")

rule string_decomposer_fasta:
	input:
		fasta = "{sample}/break_contigs/{region}/hifi_assembly_broken.fasta",
		bed = rules.string_decomposer_bed.output.bed 
	output:
		filtered_fasta = "{sample}/string_decomposer/fasta/{region}/{kmer}/{map}-{map2}/filtered_fasta.bed"
	threads: 1
	resources:
		mem = 12,
		hrs = 6,
		disk = 0
	shell:
		'''
		bedtools getfasta -fi {input.fasta} -bed {input.bed} -fo {output.filtered_fasta}
		samtools faidx {output.filtered_fasta}
		'''

rule string_decomposer_HiFi:
	input:
		fasta = rules.string_decomposer_fasta.output.filtered_fasta,
		monomers=getMonomer
	output:
		tsv = "{sample}/string_decomposer/HiFi/{map}-{map2}/{region}/{kmer}/string.tsv"
	threads: 8
	resources:
		mem =12,
		hrs = 12,
		disk = 0
	shell:
		'''
		run_decomposer.py {input.fasta} {input.monomers} -t {threads} -o {output.tsv}
		'''
rule gather_ONT_per_split:
	input:
		parts = get_ONT_per_split,
		flag = "{sample}/fasta/ONT/{map}-{map2}/{region}/{kmer}/medaka/{split}/.medaka.done"
	output:
		flag = touch('{sample}/.{map}_{map2}_{region}_{kmer}_{split}')
	threads: 1
	resources:
		mem =8,
		hrs = 12,
		disk =0

rule trigger_checkpoint:
	input:
		gather=gatherAll,
	output:
		final_output = touch("{sample}/.{kmer}_{region}_{map}-{map2}_finished")
	threads:1
	resources:
		mem =12,
		hrs = 12,
		disk = 0