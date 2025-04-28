#!/bin/bash


CC1=$1
CC2=$2
s=$3  ### read1 and read2: ${s}_R1.fastq.gz, ${s}_R2.fastq.gz
fastq1=$4
fastq2=$5
chrom_size=
ref=
#ref=${CC1}x${CC2}_genome.fa
out_dir=/output_dir/${CC1}x${CC2}_${s}


bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${CC1}x${CC2}_genome.fa ${fastq1} ${fastq2} | samtools view -bhS - > ${out_dir}/${s}.bam
##bwa mem -SP -t 16 ${ref} /proj/HuMingLab/60_mice_data/shallow_seq/MH213/${s}_240104_1.fastq.gz /proj/HuMingLab/60_mice_data/shallow_seq/MH213/${s}_240104_2.fastq.gz | samtools view -bhS - > ${out_dir}/${s}.bam
./allelic_hic_mapping_v3.sh -i ${out_dir}/${s}.bam -o ${out_dir} -c ${chrom_size}/${CC1}x${CC2}_chr_sizes.txt -p1 ${CC1} -p2 ${CC2}


## converting .pairs to .hic in order to get counts using juicer dump
./convert_hic.sh $CC1 $CC2 $s $out_dir

## grabbing the IDs from unmapped reads and get the unmapped reads from the fastq files

###add by wanying
gunzip -k ${out_dir}/${s}.phase_mix.dedup.pairs.gz
awk -F $'\t' 'FNR > 90 {print $1}' ${out_dir}/${s}.phase_mix.dedup.pairs > ${out_dir}/id.txt  
total_lines=$(wc -l < "${out_dir}/id.txt")
half_size=$((total_lines / 2))

# Randomly shuffle and split the file
shuf "${out_dir}/id.txt" | split -l "$half_size" - "${out_dir}/half_"

seqtk subseq $fastq1 ${out_dir}/half_aa > ${out_dir}/unmapped_${CC1}_1.fastq.gz
seqtk subseq $fastq2 ${out_dir}/half_aa > ${out_dir}/unmapped_${CC1}_2.fastq.gz

seqtk subseq $fastq1 ${out_dir}/half_ab > ${out_dir}/unmapped_${CC2}_1.fastq.gz
seqtk subseq $fastq2 ${out_dir}/half_ab > ${out_dir}/unmapped_${CC2}_2.fastq.gz
## mapping the unmapped reads to the parental genomes separately

bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${CC1}_genome.fa ${out_dir}/unmapped_${CC1}_1.fastq.gz ${out_dir}/unmapped_${CC1}_2.fastq.gz | samtools view -bhS - > ${out_dir}/${CC1}_${s}.bam
bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${CC2}_genome.fa ${out_dir}/unmapped_${CC2}_1.fastq.gz ${out_dir}/unmapped_${CC2}_2.fastq.gz | samtools view -bhS - > ${out_dir}/${CC2}_${s}.bam

## converting BAM to .pairs
pairtools parse -c ${chrom_size}/${CC1}_chr_sizes.txt --drop-sam ${out_dir}/${CC1}_${s}.bam -o ${out_dir}/${CC1}_${s}.pairs 
pairtools parse -c ${chrom_size}/${CC2}_chr_sizes.txt --drop-sam ${out_dir}/${CC2}_${s}.bam -o ${out_dir}/${CC2}_${s}.pairs 

## converting .pairs to .hic in order to get counts using juicer dump
./split_multi_convert_hic.sh $CC1 $CC2 $s $out_dir


## getting the total counts

python new_hic_allelic_counts.py ${out_dir}/all_${s}.phase_allele1_counts.txt ${out_dir}/all_${s}.phase_allele2_counts.txt ${out_dir}/all_${CC1}_${s}_counts.txt ${out_dir}/all_${CC2}_${s}_counts.txt /home/leeh7/iQTL/Yang_HiC/AlleliC/new_liftover/${CC1}_union.hic.bedpe_0407 /home/leeh7/iQTL/Yang_HiC/AlleliC/new_liftover/${CC2}_union.hic.bedpe_0407 /home/leeh7/iQTL/Yang_HiC/AlleliC/new_liftover/mouse_union.hic.bedpe_0407 /home/leeh7/iQTL/Yang_HiC/AlleliC/split_results/${CC1}x${CC2}_${s}_allelic_total_counts.txt


