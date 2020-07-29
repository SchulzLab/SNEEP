# SNEEP
SNp Exploration and Analysis using EPigenomics data

Dependencies:
c++ 
python3.x
bedtools
knitr and R for overview pdf
openmp?!


Compile: 
run make in src folder

Small example:
 ./src/differentialBindingAffinity_multipleSNPs -o test/ -r necessaryInputFiles/REMAnnotationModelScore_1.csv -g necessaryInputFiles/REMsEnsemblIDs_geneName.txt  necessaryInputFiles/JASPAR2020_HUMAN_transfac_P0.txt  examples/someSNPs.txt 
 
 TODO: add REMannotation file and hg38.fa
 


