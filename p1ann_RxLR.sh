#!/bin/bash

#$ -cwd
# #$ -m bes
# #$ -M $email
#$ -S /bin/bash
#$ -N part1
#$ -o part1_outv1
#$ -e part1_errv1
#$ -l mem_free=3G
#$ -V

PATH=$PATH:/raid1/home/bpp/everhars/meme/bin/
echo $host
date
#===================================================================================#
#################  SET PARAMETERS HERE  #############################################
# !! User should only change parameters within this section to have plug-&-play     #
#  Specify the fasta file (including path) that you want to analyze:
#seq=/nfs0/Grunwald_Lab/home/everhars/pram/ramorum1.fasta
#################  DO NOT EDIT BELOW   ###############################################
#====================================================================================#
mkdir part1
#getorf -minsize 90 -find 0 -sequence $seq -outseq part1/seq.faa
orfseq=/nfs0/Grunwald_Lab/home/everhars/pram/effectors/orfseq/orf4.faa

################  PART 1  #########################################################################
#  Translate to ORFs, trim to start with Signal-P, and make list of candidate effectors           # 
#  using the HMM model from Danyu Shen (January, 2014)                                            # 
###################################################################################################
#Translate the genomic sequence to six-frame open-reading frame 
# the option -find 0 will perform translation of regions between STOP codons <<<<<<
# by default, the program will find ORF's in the reverse complement of the sequence
# resultant list of AA sequences will be labeled according to whether they were found in the forward or reverse direction
#########################################################################################################

#Trims sequence so that all start with the M-signal peptide
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a $orfseq -o part1/msignalseq1.fasta

#Runs the trimmed AA sequence through the HMM software to generate candidates
###/local/cluster/spatafora/hmmer-3.0/src/hmmsearch -o part1/avh_raw.out /nfs0/Grunwald_Lab/home/everhars/pram/effectors/Avh.hmm part1/msignalseq1.fasta
############# !!!  below is an updated version of hmmer that should be used in future runs
/raid1/home/bpp/everhars/bin/hmmer-3.1b1-linux-intel-x86_64/src/hmmsearch -o part1/avh_raw.out /nfs0/Grunwald_Lab/home/everhars/pram/effectors/Avh.hmm part1/msignalseq1.fasta

# Remove candidates with hmm score <20 and those that lack BOTH the RxLR and dEER motif
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/avhfilter.pl -i part1/avh_raw.out -o part1/avhfilt.out -l part1/avhfiltlist.out

# Excludes exact site repeats using my perl script
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d part1/avhfiltlist.out -a $orfseq -o part1/avh1.fasta

# Trim sequences to start with M-signal peptide using my perl script 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part1/avh1.fasta -o part1/msignalseq2.fasta

# SignalP used to give signal peptide score:
/raid1/home/bpp/everhars/bin/signalp/signalp-3.0/signalp -t euk -d part1/sigpout.txt -f short -m hmm -q part1/msignalseq2.fasta

# Filter out sequences with signal peptide score <0.8 using my perl script
# (verify at MEME checkpoint that each contains RxLR motif; set to 0.9 if you want to avoid this check)    #
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/sigpfilt.pl -i part1/sigpout.txt -o part1/avhcandlist2.txt -v T -f 0.8

# Create fasta file of reads, using my script, from the list of sites in previous step 

# Filter out sequences with signal peptide score <0.8 using my perl script
# (verify at MEME checkpoint that each contains RxLR motif; set to 0.9 if you want to avoid this check)    #
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/sigpfilt.pl -i part1/sigpout.txt -o part1/avhcandlist2.txt -v T -f 0.8

# Create fasta file of reads, using my script, from the list of sites in previous step 
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/filterfasta.pl -d part1/avhcandlist2.txt -a $orfseq -o part1/avh2.fasta

# Trim sequences to start with M-signal peptide using my perl script
/nfs0/Grunwald_Lab/home/everhars/pram/effectors/trimfasta.pl -a part1/avh2.fasta -o part1/msignalseq3.fasta

# Use MEME to detect short motifs and weed out (in a manual step) the unlikely candidates
/raid1/home/bpp/everhars/meme/bin/meme part1/msignalseq3.fasta -minw 4 -maxw 8 -nmotifs 10 -o part1/memeout


################  MEME MANUAL CHECKPOINT #1  ######################################################################
#   Manually review the html output from MEME:
#   1. remove if length of residue is <60 AA
#   2. remove if RxLR is <20 residues from C-terminus (end)                                                       #
#   3. remove if dEER is <10 residues from C-terminus                                                             #
#   4. remove if either motif is "too close" to C- or N-terminus (not sure what "too close" is)                   #
#   5. remove if score is <1E-5                                                                                   #
#  
#   Create a cleaned up list by copying msignalseq3.fasta and removing the ones listed above and saving as msignalseq4.fasta
###################################################################################################################

################ PART 2 ###########################################################################
#  Candidate Avh list is now going to be used perform new HMM search in the original data         #
###################################################################################################

#EOF
