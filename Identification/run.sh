#!/bin/bash

# This script processes a given input file to perform various bioinformatics tasks including TPM cutoff, ORF translation, signal peptide recognition, and mature toxin liberation.
# Usage: run.sh <input_file> <output_directory_name>
# Typically I would take the input file and name output directory name the same.

# TPM cutoff
echo "Running"
input_file="$1"
name="$2"
if [ -z "$input_file" ]; then
  echo "Error: No input file provided."
  echo "Usage: $0 <input_file>"
  exit 1
fi
mkdir $2
cp $1 $2/
cd $2/
read -p "Press enter to continue"
python3 ../tpmCutoff.py $1 100
read -p "Press enter to continue"

# ORF translation
getorf -sequence aboveThreshold.fa -outseq aboveORF.fa
getorf -sequence belowThreshold.fa -outseq belowORF.fa
read -p "Press enter to continue"

# All-against-all ToxProt search for belowThreshold
mmseqs easy-search belowORF.fa ../uniprotToxProt.fasta belowSel.fa tmp --format-output "query,target,evalue,qseq" -e 1.000E-3
read -p "Press enter to continue"

# Performing signal peptide recognition by DTU/DeepTMHMM server
echo "Running signal peptide recognition (by DTU/DeepTMHMM)"
cd ../Deep/DeepTMHMM
python3 predict.py --fasta ../../$2/aboveORF.fa --output-dir $2-aboveSignal
python3 predict.py --fasta ../../$2/belowSelCut.fa --output-dir $2-belowSignal
read -p "Press enter to continue"
cd ../../$2

# Performing mature toxin liberation by ResolveSP script
# This requires that the local install of DeepTMHMM is available (replace file path with your local path)
echo "Runnning ResolveSP to attain mature sequence"
python3 ../ResolveSP.py ../Deep/DeepTMHMM/$2-aboveSignal/predicted_topologies.3line extra.fasta aboveClean.fasta
python3 ../ResolveSP.py ../Deep/DeepTMHMM/$2-belowSignal/predicted_topologies.3line extraBelow.fasta belowClean.fasta
read -p "Press enter to continue"

# Cleaning sequence for 30 length mark
python3 ../LengthCutoff.py aboveClean.fasta 30 aboveClean30.fasta
read -p "Press enter to continue"