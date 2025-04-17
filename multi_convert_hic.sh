#!/bin/bash

CC1=$1
CC2=$2
name=$3
out_dir=$4


CCs=(${CC1} ${CC2}
)

cd $out_dir

for CC in "${CCs[@]}" ; do 
pairtools parse -c /home/leeh7/iQTL/fasta/chr_size/${CC}_chr_sizes.txt --drop-sam ${CC}_${name}.bam -o ${CC}_${name}.pairs
pairtools sort -o test${CC}.sorted.pairs ${CC}_${name}.pairs;  
grep -v '#' test${CC}.sorted.pairs | awk '{if ($2 > $4){print $6"\t"$2"\t"$3"\t""0""\t"$7"\t"$4"\t"$5"\t""1"}else{print $7"\t"$4"\t"$5"\t""0""\t"$6"\t"$2"\t"$3"\t""1"}}' > tmp.pairs;
sed 's/+/0/g' tmp.pairs > test.pairs;
sed -i 's/-/1/g' test.pairs;
sort -k2,2d -k6,6d test.pairs > tmp.pairs;
java -Xmx100g -jar juicer_tools_1.22.01.jar pre -r 10000 tmp.pairs ${CC}_${name}.hic /home/leeh7/iQTL/fasta/chr_size/${CC}_chr_sizes.txt;
rm -rf test${CC}.sorted.pairs test.pairs tmp.pairs;


for((j=1;j<20;j++)); do java -Xmx20g -jar juicer_tools_1.22.01.jar dump observed NONE ${CC}_${name}.hic chr${j}_${CC} chr${j}_${CC} BP 10000 ${CC}_${name}_chr${j}.txt; done;

for((j=1;j<20;j++)); do awk -F $'\t' -v a="chr${j}" '{print a,$0}' OFS=$'\t' ${CC}_${name}_chr${j}.txt > ${CC}_${name}_chr${j}_fixed.txt; done;

cat ${CC}_${name}_chr*_fixed.txt > all_${CC}_${name}_counts.txt; done







#for((i=1;i<20;i++)); do java -Xmx20g -jar /home/leeh7/software/juicer_tools_1.22.01.jar dump observed NONE MH213.phase_allele2.hic chr${i}_CC021 chr${i}_CC021 BP 10000 MH213.phase_allele2_chr${i}.txt; done

#for((i=1;i<20;i++)); do awk -F $'\t' -v a="chr${i}_CC075" '{print a,$0}' OFS=$'\t' MH213.phase_allele1_chr${i}.txt > MH213.phase_allele1_chr${i}_fixed.txt; done
#for((i=1;i<20;i++)); do awk -F $'\t' -v a="chr${i}_CC021" '{print a,$0}' OFS=$'\t' MH213.phase_allele2_chr${i}.txt > MH213.phase_allele2_chr${i}_fixed.txt; done
