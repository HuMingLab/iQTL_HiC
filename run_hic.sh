#!/bin/bash


phase1=$1
phase2=$2
s=$3  ### read1 and read2: ${s}_R1.fastq.gz, ${s}_R2.fastq.gz
F1=$4
F2=$5
chrom_size=
ref=
#ref=/home/leeh7/iQTL/fasta/reg_bwa/${phase1}x${phase2}_genome.fa
out_dir=/output_dir/${phase1}x${phase2}_${s}


bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${phase1}x${phase2}_genome.fa ${F1} ${F2} | samtools view -bhS - > ${out_dir}/${s}.bam
##bwa mem -SP -t 16 ${ref} /proj/HuMingLab/60_mice_data/shallow_seq/MH213/${s}_240104_1.fastq.gz /proj/HuMingLab/60_mice_data/shallow_seq/MH213/${s}_240104_2.fastq.gz | samtools view -bhS - > ${out_dir}/${s}.bam
./allelic_hic_mapping_v3.sh -i ${out_dir}/${s}.bam -o ${out_dir} -c ${chrom_size}/${phase1}x${phase2}_chr_sizes.txt -p1 ${phase1} -p2 ${phase2}


## converting .pairs to .hic in order to get counts using juicer dump
./convert_hic.sh $phase1 $phase2 $s

## grabbing the IDs from unmapped reads and get the unmapped reads from the fastq files
gunzip -k ${out_dir}/${s}.phase_mix.dedup.pairs.gz
awk -F $'\t' 'FNR > 90 {print $1}' ${out_dir}/${s}.phase_mix.dedup.pairs > ${out_dir}/id.txt  
seqtk subseq $F1 ${out_dir}/id.txt > ${out_dir}/unmapped_${s}_1.fastq.gz
seqtk subseq $F2 ${out_dir}/id.txt > ${out_dir}/unmapped_${s}_2.fastq.gz

## mapping the unmapped reads to the parental genomes separately

bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${phase1}_genome.fa ${out_dir}/unmapped_${s}_1.fastq.gz ${out_dir}/unmapped_${s}_2.fastq.gz | samtools view -bhS - > ${out_dir}/${phase1}_${s}.bam
bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${phase2}_genome.fa ${out_dir}/unmapped_${s}_1.fastq.gz ${out_dir}/unmapped_${s}_2.fastq.gz | samtools view -bhS - > ${out_dir}/${phase2}_${s}.bam

## converting BAM to .pairs
pairtools parse -c ${chrom_size}/${phase1}_chr_sizes.txt --drop-sam ${out_dir}/${phase1}_${s}.bam -o ${out_dir}/${phase1}_${s}.pairs 
pairtools parse -c ${chrom_size}/${phase2}_chr_sizes.txt --drop-sam ${out_dir}/${phase2}_${s}.bam -o ${out_dir}/${phase2}_${s}.pairs 

## converting .pairs to .hic in order to get counts using juicer dump
./multi_convert_hic.sh $phase1 $phase2 $s


## getting the total counts

python hic_allelic_counts.py ${out_dir}/all_${s}.phase_allele1_counts.txt ${out_dir}/all_${s}.phase_allele2_counts.txt ${out_dir}/all_${phase1}_${s}_counts.txt ${out_dir}/all_${phase2}_${s}_counts.txt /home/leeh7/iQTL/liftover/HiC_loops/${phase1}_liftover_HiC_loops.txt /home/leeh7/iQTL/liftover/HiC_loops/${phase2}_liftover_HiC_loops.txt /home/leeh7/iQTL/liftover/111323_345582_loops_from_59mice.bedpe ${out_dir}/${phase1}x${phase2}_${s}_allelic_total_counts.txt

