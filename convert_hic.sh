#!/bin/bash

CC1=$1
CC2=$2
name=$3


cd ../${CC1}x${CC2}_${name}

for((i=1;i<3;i++)); 
do 
#xgunzip -k ${name}.phase_allele${i}.dedup.pairs.gz
pairtools sort -o test${i}.sorted.pairs ${name}.phase_allele${i}.dedup.pairs;  
grep -v '#' test${i}.sorted.pairs | awk '{if ($2 > $4){print $6"\t"$2"\t"$3"\t""0""\t"$7"\t"$4"\t"$5"\t""1"}else{print $7"\t"$4"\t"$5"\t""0""\t"$6"\t"$2"\t"$3"\t""1"}}' > tmp.pairs;          
sed 's/+/0/g' tmp.pairs > test.pairs;
sed -i 's/-/1/g' test.pairs;
sort -k2,2d -k6,6d test.pairs > tmp.pairs;
java -Xmx100g -jar /home/leeh7/software/juicer_tools_1.22.01.jar pre -r 10000 tmp.pairs ${name}.phase_allele${i}.hic /home/leeh7/iQTL/fasta/chr_size/${CC1}x${CC2}_chr_sizes.txt;
rm -rf test${i}.sorted.pairs test.pairs tmp.pairs; done


for((j=1;j<20;j++)); do java -Xmx20g -jar /home/leeh7/software/juicer_tools_1.22.01.jar dump observed NONE ${name}.phase_allele1.hic chr${j}_${CC1} chr${j}_${CC1} BP 10000 ${name}.phase_allele1_chr${j}.txt; done
for((j=1;j<20;j++)); do java -Xmx20g -jar /home/leeh7/software/juicer_tools_1.22.01.jar dump observed NONE ${name}.phase_allele2.hic chr${j}_${CC2} chr${j}_${CC2} BP 10000 ${name}.phase_allele2_chr${j}.txt; done

for((i=1;i<3;i++)); do
for((j=1;j<20;j++)); do awk -F $'\t' -v a="chr${j}" '{print a,$0}' OFS=$'\t' ${name}.phase_allele${i}_chr${j}.txt > ${name}.phase_allele${i}_chr${j}_fixed.txt; done;

cat ${name}.phase_allele${i}_chr*_fixed.txt > all_${name}.phase_allele${i}_counts.txt; done
cd ..






#for((i=1;i<20;i++)); do java -Xmx20g -jar /home/leeh7/software/juicer_tools_1.22.01.jar dump observed NONE MH213.phase_allele2.hic chr${i}_CC021 chr${i}_CC021 BP 10000 MH213.phase_allele2_chr${i}.txt; done

#for((i=1;i<20;i++)); do awk -F $'\t' -v a="chr${i}_CC075" '{print a,$0}' OFS=$'\t' MH213.phase_allele1_chr${i}.txt > MH213.phase_allele1_chr${i}_fixed.txt; done
#for((i=1;i<20;i++)); do awk -F $'\t' -v a="chr${i}_CC021" '{print a,$0}' OFS=$'\t' MH213.phase_allele2_chr${i}.txt > MH213.phase_allele2_chr${i}_fixed.txt; done
