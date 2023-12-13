#!/bin/bash

genome=$1
dbSNP=$2
interactions=$3

## minimal example
time differentialBindingAffinity_multipleSNPs examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${genome} necessaryInputFiles/estimatedScalesPerMotif_1.9.txt

## Example 1: Consider only TFs expressed in the cell type or tissue of interest
time differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_expression/ -c 0.001 -n 10 -t examples/RNA-seq_humanLV_hiPSC-CM.txt -e examples/TF_ensemblID_name_human_JASPAR2022_GRCh38p13.txt  -d 0.5  -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${genome} necessaryInputFiles/estimatedScalesPerMotif_1.9.txt

## Example 2: Add open chromatin regions of the cell type of interest

#download encode data 
wget 'https://www.encodeproject.org/files/ENCFF199VHV/@@download/ENCFF199VHV.bed.gz'
gunzip ENCFF199VHV.bed.gz

time differentialBindingAffinity_multipleSNPs  -o examples/SNEEP_output_open_chromatin/ -c 0.001 -n 10 -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt -f ENCFF199VHV.bed examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${genome} necessaryInputFiles/estimatedScalesPerMotif_1.9.txt

## Example 3: Associate the SNPs, which significantly affect the binding behavior of a TF to their target genes
time differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_REM_PRO_HiC/ -c 0.001 -n 10  -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt  -r ${interactions} -g ensemblID_GeneName.txt  examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${genome} necessaryInputFiles/estimatedScalesPerMotif_1.9.txt

## Example 4: Compute a proper random background control and highlight the cell-ype specific TFs
time differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_background_sampling/ -c 0.001 -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt  -n 20 -j 100 -k ${dbSNP} -l 2 -q 0 -r interactionsREM_PRO_HiC.txt -g ensemblID_GeneName.txt  examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${genome} necessaryInputFiles/estimatedScalesPerMotif_1.9.txt

