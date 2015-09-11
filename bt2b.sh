#!/bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -N bt2b
#$ -o bt2bout
#$ -e bt2berr
#$ -l mem_free=10G
#$ -V
# #$ -h
#$ -t 1-3:1

i=$(expr $SGE_TASK_ID - 1)

PATH=$PATH:/raid1/home/bpp/knausb/bin/samtools-0.1.18/bcftools/

FILE=( `cat "samplesb.txt" `)

#for j in "${FILE[@]}"
#do
# echo $j
#done

REF="/raid1/labs/Grunwald_Lab/home/everhars/pram/GBS/bowtie2/ramorum"
#REF='/raid1/labs/Grunwald_Lab/web_data/pinf_broad/mt/pinf_mtiia'

IFS=';' read -a arr <<< "${FILE[$i]}"
#IFS=',' read -a arr <<< "$FILE[0]"

echo ${arr[0]}
#echo "${arr[1]}"

echo -n "Running on: "
hostname
echo "SGE job id: $JOB_ID"
date

CMD="/home/bpp/knausb/bin/bowtie2-2.0.6/bowtie2 -q --no-unal --rg-id ${arr[0]}.sra --rg SM:${arr[0]} --local $REF -U ${arr[1]} -S sams/${arr[0]}.sam"
echo $CMD
$CMD

echo "Bowtie2 done"
date

#cd sams

pwd
#ls

#CMD="samtools view -S -b ./sams/${arr[0]}.sam > ./bams/${arr[0]}.bam"
#CMD="samtools view -bS sams/${arr[0]}.sam > /raid1/labs/Grunwald_Lab/home/knausb/pinf_bowtie2/yoshida/${arr[0]}.bam"
CMD="samtools view -bS -o bams/${arr[0]}.bam sams/${arr[0]}.sam"
echo $CMD
$CMD
CMD="samtools sort bams/${arr[0]}.bam bams/${arr[0]}.sorted"
echo $CMD
$CMD
CMD="samtools index bams/${arr[0]}.sorted.bam"
echo $CMD
$CMD

echo "Samtools done"
date

# EOF.
