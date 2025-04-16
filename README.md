# iQTL_HiC
## Allelic HiC mapping for CCGG.



### Requirements
- BWA-mem2 v2.2.1
- Samtools v1.20
- Pairtools v1.0.2
- Seqtk v1.4
- Juicer v.1.22.01
- Python v3.6.8
	* Numpy v1.19.5 
	* Pandas v1.1.5

### Running the pipeline
This pipeline can be simply run by 
```
./allelic_hic_mapping_v3.sh [Mouse Name] [Maternal Genome Name] [Paternal Genome Name] [fastq R1] [fastq R2]
```

For rest of this README, we define each of the main inputs (required) as:
1. Mouse Name = ${NAME}
2. Maternal Genome Name = ${CC1}
3. paternal Genome Name = ${CC2}
4. fastq R1 = ${F1}
5. fastq R2 = ${F2}
6. ${OUT} = output directory
7. ${map_dir} = working directory

Here are the detailed steps to this pipeline:
Step 1: Mapping fastq file to the psuedo-genome
	Input: pseudo-genome fasta, fastq R1/R2 files
	Output: mapped reads in BAM file
```
bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${CC1}x${CC2}_genome.fa $fastq1 $fastq2 | samtools view -bhS - > ${map_dir}/${s}.bam                                              
```
Step 2: Convert the mapped pair-end reads into ligation events
	Input: mapped BAM from step 1, chromosome size file
	Output: unphased .pairs file
```
pairtools parse --min-mapq 0 --add-columns XA,NM,AS,XS --nproc-in 16 --nproc-out 16 --drop-sam --walks-policy all -c ${chrom_size}/${CC1}x${CC2}_chr_sizes.txt -o ${map_dir}/${s}.unphase.pairs.gz ${map_dir}/${s}.bam
```

Step 3: Assign unphased pairs file into the allele 1, allele 2, unmapped (mixed), or bad pairs
	Input: unphased pairs file from step 2
	Output: 4 .pairs files for: allele 1, allele 2, unmapped (mixed), or bad
 ```
python allelic_mphase_v2.py --in ${map_dir}/${s}.unphase.pairs.gz --out ${map_dir}/${s}.phase --phase1 $CC1 --phase2 $CC2
```
Step 4: Deduplicate and grab QC
	Input: 4 .pairs file from step 3
	Output: 4 deduplicated .pairs for: allele 1, allele 2, unmapped (mixed), or bad file and QC files ending in .stats and .pairs.info

```
for f in ${map_dir}/${s}.phase*pairs
do
        fname=`basename $f .pairs`
        bgzip ${f}
        pairtools sort -o ${map_dir}/${fname}.sorted.pairs.gz ${f}.gz
        pairtools dedup --mark-dups --extra-col-pair phase1 phase2 --output-stats ${map_dir}/${fname}.dedup.stats -o ${map_dir}/${fname}.dedup.pairs.gz ${map_dir}/${fname}.sorted.pairs.gz
        zcat ${map_dir}/${fname}.dedup.pairs.gz | grep -v '#' | cut -f17-18 | sort | uniq -c > ${map_dir}/${fname}.dedup.pairs.info
done
```
Step 5: Converting .pairs to .hic in order to get counts using juicer dump
Input: 2 deduplicated allele pairs file, chromosome size file
Output: maternal and paternal uniquely mapped count files
```
./convert_hic.sh $CC1 $CC2 $s
```
In the shellscript, it sorts and modifies the .pairs file in the way that Juicer requires the .pairs file to be. Then, we use juicer pre to convert .pairs to .hic and use Juicer dump to get all the counts in a txt file. 
Step 6: Grabbing the IDs from unmapped reads and get the unmapped reads from the fastq files
	Input: deduplicated unmapped (mix) pairs file
	Output: unmapped fastq R1 and R2 files
```
gunzip -k ${map_dir}/${s}.phase_mix.dedup.pairs.gz
awk -F $'\t' 'FNR > 90 {print $1}' ${map_dir}/${s}.phase_mix.dedup.pairs > ${map_dir}/id.txt  
seqtk subseq $fastq1 id.txt > ${map_dir}/unmapped_${s}_1.fastq.gz
seqtk subseq $fastq2 id.txt > ${map_dir}/unmapped_${s}_2.fastq.gzStep 7: Map the unmapped fastq files to maternal and paternal genomes separately
```	
Step 7: Mapping the unmapped reads to the parental genomes separately
	Input: unmapped fastq files, maternal and paternal reference genome files (fasta)
	Ouput: multi-mapped maternal and paternal reads in BAM file
```
bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${CC1}_genome.fa ${map_dir}/unmapped_${s}_1.fastq.gz ${map_dir}/unmapped_${s}_2.fastq.gz | samtools view -bhS - > ${map_dir}/${CC1}_${s}.bam
bwa-mem2 mem -SP5M -T0 -t16 ${ref}/${CC2}_genome.fa ${map_dir}/unmapped_${s}_1.fastq.gz ${map_dir}/unmapped_${s}_2.fastq.gz | samtools view -bhS - > ${map_dir}/${CC2}_${s}.bam
```
Step 8: Converting BAM to .pairs
	Input: multi-mapped maternal and paternal reads in BAM files from step 7
	Output: multi-mapped maternal and paternal reads in .pairs format
```
pairtools parse -c ${chrom_size}/${CC1}_chr_sizes.txt --drop-sam ${map_dir}/${CC1}_${s}.bam -o ${map_dir}/${CC1}_${s}.pairs 
pairtools parse -c ${chrom_size}/${CC2}_chr_sizes.txt --drop-sam ${map_dir}/${CC2}_${s}.bam -o ${map_dir}/${CC2}_${s}.pairs
```



Step 9: Converting .pairs to .hic in order to get counts using juicer dump
	Input: multi-mapped maternal and paternal reads in .pairs format
	Output: maternal and paternal multi mapped count files
```
./multi_convert_hic.sh $CC1 $CC2 $s
```
The same process is being done as the step 5.

Step 9: Merge the uniquely and multi-mapped counts and calculate the total counts
	Input: the count files for maternal/paternal uniquely and multi-mapped reads, liftover HiC loops
	Output: the total and allelic-counts (Note that the rows that contain all 0s or 0s for CC1 or CC2 positions, these are meant to be NA – including the 0 count respectively) 
``` 
Python hic_allelic_counts.py [maternal uniquely mapped counts] [paternal uniquely mapped counts] [maternal multi-mapped counts] [paternal multi-mapped counts] [maternal liftover HiC loops] [paternal liftover HiC loops] [mm10 based HiC loops] [output file name]
```

	

This pipeline was developed with [Dr. Yang Xie](https://github.com/Xieeeee/AlleliC/). Contact us at Ming Hu (hum@ccf.org) or Lindsay Lee (leeh7@ccf.org).
