#!/bin/bash

#$ -cwd
# #$ -m bes
# #$ -M $email
#$ -S /bin/bash
#$ -N part4
#$ -o part4_outv3
#$ -e part4_errv3
#$ -l mem_free=3G
#$ -V


############  DO NOT EDIT BELOW THIS LINE  ###################################################
#============================================================================================#
mkdir part4
touch part4/avh2blast2.list
        ls -1 part3/avh2blast2*.list >> part4/avh2blast2.list
ablast2=( `cat "part4/avh2blast2.list" `)
#touch part4/avh2blast2rem.list
#        ls -1 part3/avh2removedblast2*.list >> part4/avh2removedblast2.list
#ablastrem=( `cat "part4/avh2removedblast2.list" `)

####  Echo files that are used in this script
for t in "${ablast2[@]}"
        do
                echo $t
        done

#for t in "${ablastrem[@]}"
#	do
#		echo $t
#	done

#######  Summarizing the number of candidate Avh candidates are identified in each blast search ###########
touch part4/summary.out
sum=part4/summary.out
#ls -1 *list >> $sum
	# summarizing for blast to ramorum reference results:
#	echo avh2blast2ref.list >> $sum
#	sort part3/avh2blast2ref.list | wc -l >> $sum
#	sort -u part3/avh2blast2ref.list | wc -l >> $sum
#	echo avh2removedblast2ref.list >> $sum
#	sort part3/avh2removedblast2ref.list | wc -l >> $sum
#	sort -u part3/avh2removedblast2ref.list | wc -l >> $sum
#	# summarizing for blast to avh family results:
#	echo avh2blast2fam.list >> $sum
#	sort part3/avh2blast2fam.list | wc -l >> $sum
#	sort -u part3/avh2blast2fam.list | wc -l >> $sum
#	echo avh2removedblast2fam.list >> $sum
#	sort part3/avh2removedblast2fam.list | wc -l >> $sum
#	sort -u part3/avh2removedblast2fam.list | wc -l >> $sum

for t in "${ablast2[@]}"
	do
		echo $t >> $sum
		sort $t | wc -l >> $sum
		sort -u $t | wc -l >> $sum
	done

#for t in "${ablastrem[@]}"
#	do
#		echo $t >> $sum
#		sort $t | wc -l >> $sum
#		sort -u $t | wc -l >> $sum
#	done


#EOF
