#!/bin/bash

#$ -cwd
# #$ -m bes
# #$ -M $email
#$ -S /bin/bash
#$ -N part2
#$ -o part2_outv1
#$ -e part2_errv1
#$ -l mem_free=3G
#$ -V

PATH=$PATH:/raid1/home/bpp/everhars/meme/bin/
echo $host
date
################  MEME MANUAL CHECKPOINT #1  ######################################################################
#   Manually review the html output from MEME:
#   1. remove if length of residue is <60 AA
#   2. remove if RxLR is <20 residues from C-terminus (end)                                                       #
#   3. remove if dEER is <10 residues from C-terminus                                                             #
#   4. remove if either motif is "too close" to C- or N-terminus (not sure what "too close" is)                   #
#   5. remove if score is <1E-5                                                                                   #
#  
#   Create a cleaned up list by copying msignalseq3.fasta and removing the ones listed above and saving as msignalseq4.fasta in part1 directory
###################################################################################################################
#################  DO NOT EDIT BELOW   ###############################################
#====================================================================================#
mkdir part2
orfseq=part1/seq.faa
orfseq=/nfs0/Grunwald_Lab/home/everhars/pram/effectors/orfseq/orf4.faa
#  amino acid translations of other genomes for blast comparison  
echo "$orfseq"

################ PART 2 ###########################################################################
#  Candidate Avh list is now going to be used perform new HMM search in the original data         #
###################################################################################################
###Create a multiple sequence alignment of the avh candidates in a new folder named part2/ and named msignalseq4.fasta
muscle -in part1/msignalseq4.fasta -out part2/msig4algnd.fasta

###Build HMM model of aligned amino acid sequences   
hmmbuild --amino --seed 854 -o part2/hmmbuild.out part2/Avh_new.hmm part2/msig4algnd.fasta

###Runs the original trimed AA sequence thorugh the HMM software using the new Avh database - have to use new version to avoid error
## don't use this line unless the next line throws an error: /raid1/home/bpp/everhars/bin/hmmer-3.1b1-linux-intel-x86_64/src/hmmsearch -o avh_round2.out /nfs0/Grunwald_Lab/home/everhars/pram/effectors/smEU1_030002/Avh_new.hmm msignalseq1.fasta
/raid1/home/bpp/everhars/bin/hmmer-3.1b1-linux-intel-x86_64/src/hmmsearch -o part2/avh_round2.out part2/Avh_new.hmm part1/msignalseq1.fasta

### Remove candidates with hmm score <10 and those that lack BOTH the RxLR and dEER motif  			******* This is a low-threshold hmm score filter, normally it's <20  **************
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/avhfilterless.pl -i part2/avh_round2.out -o part2/avhfilt_round2.out -l part2/avhfiltlist_round2.out

# Excludes exact site repeats using my perl script
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d part2/avhfiltlist_round2.out -a $orfseq -o part2/avh1_round2.fasta

# Trim sequences to start with M-signal peptide using my perl script 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part2/avh1_round2.fasta -o part2/msignalseq2_round2.fasta

# SignalP used to give signal peptide score:
/raid1/home/bpp/everhars/bin/signalp/signalp-3.0/signalp -t euk -d part2/sigpout_round2.txt -f short -m hmm -q part2/msignalseq2_round2.fasta

## Filter out sequences with signal peptide score <0.00001 using my perl script
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/sigpfilt.pl -i part2/sigpout_round2.txt -o part2/avhcandlist2_round2.txt -v T -f 0.00001

# Sequences with poor (<0.8) are kept with this filter script
#/nfs0/Grunwald_Lab/home/everhars/pram/effectors/sigpfiltpoor.pl -i part2/sigpout_round2.txt -o part2/avhcandlist2_round2_poor.txt -v T -f 0.8

# Create fasta file of reads, using my script, from the list of sites in previous step 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d part2/avhcandlist2_round2.txt -a $orfseq -o part2/avh2_round2.fasta
#/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d part2/avhcandlist2_round2_poor.txt -a $orfseq -o part2/avh2_round2_poor.fasta

# Trim sequences to start with M-signal peptide using my perl script
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part2/avh2_round2.fasta -o part2/msignalseq3_round2.fasta
#/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part2/avh2_round2_poor.fasta -o part2/msignalseq3_round2_poor.fasta

# Use MEME to detect short motifs and weed out (in a manual step) the unlikely candidates
##### NOTE:  in test running this script, meme throws an error of too many candidates, so increased -minsize to 200000 from 100000
/raid1/home/bpp/everhars/meme/bin/meme part2/msignalseq3_round2.fasta -minw 4 -maxw 8 -nmotifs 4 -minsize 200000 -o part2/memeout_round2
#/raid1/home/bpp/everhars/meme/bin/meme part2/msignalseq3_round2_poor.fasta -minw 4 -maxw 8 -nmotifs 10 -o part2/memeout_round2_poor

####  MEME MANUAL CHECKPOINT #2   ###############################################################################################
#   Manual checkpoint #2 -- performed in Excel
#   Exact criteria were:  
#   1. Removed both dEER and RxLR candidate if <10AA at either N or C-terminus
#   2. Also removed RxLR with pval<1E-4
#   Created list sites to be removed in Excel, list of sites to be removed named removelist.txt 
######################################################################################################################


#EOF
