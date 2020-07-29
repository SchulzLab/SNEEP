# SNEEP
SNp Exploration and Analysis using EPigenomics data

# Dependencies:
- c++ 
- python3.x
- bedtools
- knitr and R for overview pdf
- openmp?!


# Compile: 
To compile the C++ code of SNEEP run the following 
```
cd src/
make
```

# Small example:
For the help information type
```
./src/differentialBindingAffinity_multipleSNPs -h
```

To run a small example, perform the following command
```
 ./src/differentialBindingAffinity_multipleSNPs -o test/ -r necessaryInputFiles/REMAnnotationModelScore_1.csv -g necessaryInputFiles/REMsEnsemblIDs_geneName.txt  necessaryInputFiles/JASPAR2020_HUMAN_transfac_P0.txt  examples/someSNPs.txt 
 ```
# TODOs: 
- add REMannotation file and hg38.fa to necessaryInputFiles
 


