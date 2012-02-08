#!/bin/bash

#-------------------------------------------
# Parameters

# Transcript annotation (BED file)
BED=sample.bed

# Read length
READLEN=75

# Number of reads generated
NREAD=100000

# The mean and std of spans for paired-end reads. 
PAIREDEND="200,20"

# positional bias file
POSBIAS=sampleposbias.txt

# Read error position profile
READERR=samplereaderror.txt

# output FASTA prefix
FASTAFILE=outpair

# reference chromosome
REFERENCE=reference.fa

#-----------------------------------------------
# Commands to randomly assign weights to each transcript

CMD0="sort -k 1,1 -k 2,2n $BED | genexplvprofile.py > explvprofile.txt"

echo "Commands to randomly assign weights to each transcript:"
echo $CMD0

sort -k 1,1 -k 2,2n $BED | genexplvprofile.py > explvprofile.txt

# Commands to simulate reads (output to STDOUT in BED format)  
# If you want single-end reads, don't use the "-p" option.
CMD1="gensimreads.py -e explvprofile.txt  -n $NREAD -b $POSBIAS -l $READLEN -p $PAIREDEND  $BED "

echo "Commands to generate simulated paired-reads in BED format:"
echo $CMD1


# Commands to convert BED file to fasta file
CMD2="getseqfrombed.py -b $READERR -f A -r 0.01 -l $READLEN -  $REFERENCE"

echo "Commands to generate FASTA file from the last command:"
echo $CMD2
echo "Output FASTA prefix: $FASTAFILE"



# Execute two commands simultaneously
# If you want single-end reads, do not use the splitfasta.py. Instead, execute:
# $CMD1 | $CMD2 > $FASTAFILE
$CMD1  | $CMD2 | splitfasta.py -o $FASTAFILE 




