# SNEEP
SNp Exploration and Analysis using EPigenomics data

# Dependencies:
- c++11 
- g++ (9.3.0)
- python3.x
- bedtools (v2.27.1)
- knitr and R for overview pdf (details?!)
- openmp
- R (4.0.4), ggplot2 library


# Compile: 
To compile our SNEEP code, run the following 
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
 


