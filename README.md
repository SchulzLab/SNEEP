# SNEEP: SNP exploration and functional analysis using epigenetic data

SNEEP is a method that helps to prioritize GWAS SNPs to study the impact of genetically induced transcriptional mis-regulation in human diseases and other phenotypes.

In more detail, our SNEEP approach prioritizes SNPs as targets of one or several Transcription Factors (TFs) and infers whether a gene’s expression is influenced by the change in the TF binding behavior. To do so, we 
-	evaluate the impact of a SNP on a potential TF binding site by calculating a probabilistic differential binding score for the difference in TF binding affinity in wild type versus a mutated sequence. 
-	associate SNPs to potential target genes, by making use of different types of information, like linked epigenetic elements using Hi-C data or catalouged elements that associate with gene expression changes from the EpiRegio database. 

SNEEP can easily handle large collections of SNPs provided by the user. Additional, various kinds of user-specific epigenetic data can be incorporated, like 
-	open chromatin data to exclude SNPs located in region not accessible in the cell type of interest 
-	gene expression data to restrict the analysis to TFs which are expressed. 

Further, SNEEP provides a comprehensive report containing different summary statistics, which
-	highlight TFs whose binding affinity is most often affected by the analyzed SNPs, 
-	list TFs associated to a gain, or a loss of TF binding affinity based on the input data,
-	provide information about how many SNPs are linked to a target gene.  

All statistical assessments are compared against proper random controls to highlight only biologically interesting results.

The following overview figure visualizes the various ways how to use our SNEEP approach. 

TODO: create figure

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
 


