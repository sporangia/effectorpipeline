#!/bin/bash

#$ -cwd
# #$ -m bes
# #$ -M $email
#$ -S /bin/bash
#$ -N NA1-part3
#$ -o part3_outv3
#$ -e part3_errv3
#$ -l mem_free=3G
#$ -V

PATH=$PATH:/raid1/home/bpp/everhars/meme/bin/
echo $host
date
####  MEME MANUAL CHECKPOINT #2   ################################################################
#   Manual checkpoint #2 -- performed in Excel
#   Exact criteria were:  
#   1. Removed both dEER and RxLR candidate if <10AA at either N or C-terminus
#   2. Also removed RxLR with pval<1E-4
#   Created list sites to be removed in Excel, list of sites to be removed named removelist.txt and save it in part2 directory 
#====================================================================================#
#################  DO NOT EDIT BELOW   ###############################################
#===================================================================================#
#orfseq=part1/seq.faa
orfseq=/nfs0/Grunwald_Lab/home/everhars/pram/effectors/orfseq/orf4.faa
#orfref=/nfs0/Grunwald_Lab/home/everhars/pram/effectors/orfseq/orfref.faa
#orffam=/nfs0/Grunwald_Lab/home/everhars/pram/effectors/orfseq/orffam.faa

# all amino acid translations of genomes for blast comparison  -- including reference and family 
mkdir part3
touch part3/orfseq.list
ls /nfs0/Grunwald_Lab/home/everhars/pram/effectors/orfseq/*faa >> part3/orfseq.list
orfs=( 'cat "part3/orfseq.list" ')

####  Echo files that are used in this script
echo "$orfseq"
#echo "$orfref"
#echo "$orffam"
touch part3/to-be-blasted.list
j=1
orfs=( `cat "part3/orfseq.list" `)
	for t in "${orfs[@]}"
        do
		echo "orfseq$j" >> part3/to-be-blasted.list
                echo $t >> part3/to-be-blasted.list
		j=$(($j+1))		 
        done
###########################################################################################
## Adding line to filter out sites from previous list of sites
#touch part3/avhcandlist2_reduced.txt
#grep -v -f part2/removelist.txt part2/avhcandlist2_round2.txt >> part3/avhcandlist2_reduced.txt
#touch part3/sitequickfilt.out
#printf "sitesbefore\n" >> part3/sitequickfilt.out
#wc -l part2/avhcandlist2_round2.txt >> part3/sitequickfilt.out
#printf "sitestoremove\n" >> part3/sitequickfilt.out
#wc -l part2/removelist.txt >> part3/sitequickfilt.out 
#printf "sitesafter\n" >> part3/sitequickfilt.out
#wc -l part3/avhcandlist2_reduced.txt >> part3/sitequickfilt.out

##########  NOTES    ##############################################################################
# Now there are two lists that are going to be used for subsequent blast searches
# Can also run MEME again with those removed sites if you want it for a summary output
###################################################################################################

########  Making blast database of good candidates list  ########################################################
# Create fasta file of reads, using my script, from the list of sites in previous step 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d keepall.txt -a $orfseq -o part3/avh2_forblast.fasta
# Trim sequences to start with M-signal peptide using my perl script 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part3/avh2_forblast.fasta -o part3/msignalseq2_forblast.fasta
# Convert fasta list of avh candidates into a blastdb
makeblastdb -in part3/msignalseq2_forblast.fasta -dbtype 'prot' -out part3/avh2blastdb -title avh2blastdb -parse_seqids

#########   Making blast database of poor candidates list  #######################################################
# Create fasta file of reads, using my script, from the list of sites in previous step 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d removelist.txt -a $orfseq -o part3/removelist_forblast.fasta
# Trim sequences to start with M-signal peptide using my perl script 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part3/removelist_forblast.fasta -o part3/msignal_removelist_forblast.fasta
# Convert fasta list of avh candidates into a blastdb
makeblastdb -in part3/msignal_removelist_forblast.fasta -dbtype 'prot' -out part3/avh2removedblastdb -title avh2removedblastdb -parse_seqidsi

########  Blast searching both good and poor candidates against other assembled genomes, reference genome, and Avh database
touch blast_key.list
j=1
for t in "${orfs[@]}"
	do 
	# good candidates blast
		blastp -query $t -db part3/avh2blastdb -out part3/avh2blast2orf$j.out -outfmt "6 sseqid" -evalue 1e-5
		blastdbcmd -db part3/avh2blastdb -dbtype prot -entry_batch part3/avh2blast2orf$j.out -outfmt %a -out part3/avh2blast2orf$j.list
	# poor candidates blast
		blastp -query $t -db part3/avh2removedblastdb -out part3/avh2removedblast2orf$j.out -outfmt "6 sseqid" -evalue 1e-5
		blastdbcmd -db part3/avh2removedblastdb -dbtype prot -entry_batch part3/avh2removedblast2orf$j.out -outfmt %a -out part3/avh2removedblast2orf$j.list
	# creating a list of names corresponding to blast search		
		echo "orf$j.faa" >> part3/blast_key.list
		echo $t >> part3/blast_key.list
		j=$(($j+1))
	done


## reference search
#		blastp -query $orfref -db part3/avh2blastdb -out part3/avh2blast2ref.out -outfmt "6 sseqid" -evalue 1e-5
#		blastdbcmd -db part3/avh2blastdb -dbtype prot -entry_batch part3/avh2blast2ref.out -outfmt %a -out part3/avh2blast2ref.list
#		blastp -query $orfref -db part3/avh2removedblastdb -out part3/avh2removedblast2ref.out -outfmt "6 sseqid" -evalue 1e-5
#		blastdbcmd -db part3/avh2removedblastdb -dbtype prot -entry_batch part3/avh2removedblast2ref.out -outfmt %a -out part3/avh2removedblast2ref.list
## Avh database search
#		blastp -query $orffam -db part3/avh2blastdb -out part3/avh2blast2fam.out -outfmt "6 sseqid" -evalue 1e-5
#		blastdbcmd -db part3/avh2blastdb -dbtype prot -entry_batch part3/avh2blast2fam.out -outfmt %a -out part3/avh2blast2fam.list
#		blastp -query $orffam -db part3/avh2removedblastdb -out part3/avh2removedblast2fam.out -outfmt "6 sseqid" -evalue 1e-5
#		blastdbcmd -db part3/avh2removedblastdb -dbtype prot -entry_batch part3/avh2removedblast2fam.out -outfmt %a -out part3/avh2removedblast2fam.list

#EOF
