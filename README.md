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

# Build our tool

Please make sure that the following software is available on your machine: 

- c++11 
- g++ (9.3.0)
- python3.x
- bedtools (v2.27.1)
- knitr and R for overview pdf (details?!)
- openmp
- R (4.0.4), ggplot2 library

To install SNEEP,  run the following commands: 

```
cd SNEEP/src/
make
```

TODO bash script with knitr packages 
TODO: small test script with all combinations 


# Usage 

The following 3 files are required as minimal input to run SNEEP:

1)	a file containing the TF motifs in TRANSFAC format. (see also …) 
2)	a bed-like SNP file
3)	a reference genome file in fasta format

We provide the human TF motifs from the JASPAR database (version 2020) in the required format in the examples directory (SNEEP/examples/ JASPAR2020_HUMAN_transfac_P0.txt ).  However, ever set of TF motifs can be used instead like from another database (e.g. HOCOMOCO or …) or a different species. 
The bed-like SNP file needs is a tab-separated file containing the following entries: 
-	chr
-	start
-	end
-	var1 (e.g. effector allele or alternative allele) 
-	var2 (e.g. wild type allele)
-	rsID if known, otherwise - 
-	minor allele frequency (MAF) if known, otherwise -1 (explain why we need MAF)

A properly format SNP files looks as following: 

```
chr1    109274569       109274570       G       A       rs7528419       0.2009
chr1    109275907       109275908       C       T       rs646776        0.2384
chr1    154424939       154424940       G       T       rs12118721      1e-07
chr1    154424939       154424940       G       T       -      0.3
chr12   111569951       111569952       G       C       rs653178        -1
```

If you want to consider a SNP which has multiple alternative alleles like for instance rs11206510 (https://www.ncbi.nlm.nih.gov/snp/rs11206510), please add one line per alternative allele in the bed-like SNP file. An example is shown below: 

```
chr1    55030365        55030366        A       T       rs11206510      0.1018
chr1    55030365        55030366        C       T       rs11206510      0.1018
chr1    55030365        55030366        G       T       rs11206510      0.1018
```


In the reference genome file, the different chromosome must be named as chr1, chr2 etc. resulting in the following format: 

```
>chr1
ATCGGGTCA…
>chr2
TTTGAGACCAT…
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
 


