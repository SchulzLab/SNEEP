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

We provide the human TF motifs from the JASPAR database (version 2020) in the required format in the examples directory (SNEEP/examples/ JASPAR2020_HUMAN_transfac_P0.txt ).  However, every set of TF motifs can be used instead, for instance from another database (e.g. HOCOMOCO or Kellis) or a different species. 
The bed-like SNP file needs is a tab-separated file containing the following entries: 
-	chr
-	start position (0-based)
-	end position (0-based)
-	var1 (e.g. effector allele or alternative allele) 
-	var2 (e.g. wild type allele)
-	rsID if known, otherwise - 
-	minor allele frequency (MAF) if known, otherwise -1 (explain why we need MAF)

A properly formated SNP files looks as following: 

```
chr1    109274569       109274570       G       A       rs7528419       0.2009
chr1    109275907       109275908       C       T       rs646776        0.2384
chr1    154424939       154424940       G       T       rs12118721      1e-07
chr1    154424939       154424940       G       T       -      0.3
chr12   111569951       111569952       G       C       rs653178        -1
```

If you want to consider a SNP, which has multiple alternative alleles like for instance [rs11206510](https://www.ncbi.nlm.nih.gov/snp/rs11206510) (T -> A,C,G) , please add one line per alternative allele in the bed-like SNP file. An example is shown below: 

```
chr1    55030365        55030366        A       T       rs11206510      0.1018
chr1    55030365        55030366        C       T       rs11206510      0.1018
chr1    55030365        55030366        G       T       rs11206510      0.1018
```

SNEEP can only handle mutations effecting a single base pair (no deletions or insertions).

In the reference genome file, the different chromosome must be named as chr1, chr2 etc. resulting in the following format: 

```
>chr1
ATCGGGTCA…
>chr2
TTTGAGACCAT…
```

## Minimal example: 

To try SNEEP with the minimal required input, make sure you are in the SNEEP folder and run: 

```
./src/differentialBindingAffinity_multipleSNPs  examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path-to-genome-file> 
```

Per default the result is stored in the directory ‘SNEEP_output’. For more details about the result files, see Section … 

## Optional input parameters

To display all optional input parameters, type 

```
./src/differentialBindingAffinity_multipleSNPs  -h
```

which prints

```
Call program with ./src/differentialBindingAffinity_multipleSNPs

optinal parameters:
-o outputDir (default SNEEP_output/, if you want to specific it, it must be done as first argument)
-n number threads (default 1)
-p pvalue for motif hits (default 0.05)
-c pvalue differential binding (default 0.01)
-b base frequency for PFMs -> PWMs (default; /necessaryInputFiles/frequency.txt)
-a if flag is set,  all computed differential bindinding affinities are stored in <outputDir>/AllDiffBindAffinity.txt
-f additional footprint/open chromatin region file in bed file format
-m if flag is set, the  maximal differential binding affinity per SNP is stored in <outputDir>/MaxDiffBindingAffinity.txt
-t file where expression values of TFs are stored (e.g RNA-seq in a tab-seperated format e.g. ensemblID	expression-value)
-d threshold TF activity (must be given if -t is given)
-e tab-seperated file containing ensemblID to gene name mapping of the TFs (must be given if -t is given)
-r bed-like file with epigenetic interactions
-g path to file containing ensemblID to gene name mapping, must be given if -r is given (,-seperated)(mapping for all genes within EpiRegio)
-j rounds sampled background (default 0)
-k path to sorted dbSNP file (if our provided file is used only SNPs in coding regions are considered)
-i path to the source GitHub dir (default .)
-l start seed (default 1)
-q minimal TF count which needs to be exceeded to be considered in random sampling (default 0)
-h help

transfac PFM file,  bed-like SNP file and path to genome file (fasta format)  must be given
help function end

```

# Use SNEEP on a realistic example


For an realistic example we consider SNPs associated to myocardial infarction (downloaded from the [GWAS catalog](https://www.ebi.ac.uk/gwas/efotraits/EFO_0000612) and the corresponding proxy SNPs (determined with [SNiPA](https://snipa.helmholtz-muenchen.de/snipa3/index.php?task=proxy_search), R2 value >= 0.8). The following section provides example runs with different combination of optional input parameters. The example data is located in the directory SNEEP/example/. 

## Example 1: only consider TFs expressed in the cell type or tissue of interest

To do so, we need to set the optional parameters -t, -d and -e. For our example the file containing the expression values (flag -t) are provided in the example directory and derived from cardio myocytes. Additionally, we provide the file containing the mapping between ensembl ID (used in the expression value file) and the names of the TFs specified in the motif file (…). As flag d, we use a rather less stringent expression value threshold of 0.5. In general, you can choose any value which is most suitable for you.

So, the resulting command is: 

```
./src/differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_expression/ -t examples/RNA-seq_humanLV_hiPSC-CM.txt -e examples/ensemblID_geneName_TFs.txt -d 0.5 examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed < path-to-genome-file>
```

Note, that we specified the output directory with the -o flag as examples/SNEEP_output_expression/. 

## Example 2: open chromatin regions

If open chromatin data of your cell type of interest is available, it is possible to integrate this data and automatically exclude SNPs from the analysis in closed, most likely inactive chromatin regions. 
Therefore, a bed-file holding the open chromatin regions can be specified using the flag -f. 

For our example, we want to use an ATAC-seq on human heart right ventricle from ENCODE. 

To download the data run: 

```
wget ‘https://www.encodeproject.org/files/ENCFF199VHV/@@download/ENCFF199VHV.bed.gz’
```

Next unzip the file via gunzip.

The resulting SNEEP call is 

```
./src/differentialBindingAffinity_multipleSNPs  -o examples/SNEEP_output_open_chromatin/ -f <path-to-ENCODE-data> examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path-to-genome-file>
```

## Example 3: Associate the SNPs, which significantly affect the binding behavior of a TF to their target genes

To associate the target genes, we need to provide a file holding information about epigenetic interactions (flag -e). We provide this data via a Zenodo repository which contains three different epigenetic interaction files (for more detail explanation see …). For our example the most suitable one is the file interactionsREM_PRO_HiC.txt. The HiC data is retrieved from whole human heart, so we can benefit from the interactions for our example analysis. Additionally, the file ensemblID_GeneName.txt containing the ensembl ID to gene name mapping for all genes listed in the epigenetic interaction file is required (flag -g).
 

```
./src/differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_REM_PRO_HiC/ -r interactionsREM_PRO_HiC.txt -g ensemblID_GeneName.txt  examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed /home/nbaumgarten/hg38.fa
```

