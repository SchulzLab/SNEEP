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

# Download Zenodo repository
Additionally, to the GitHub repository, which contains the source code and small example files, we created a Zenodo repository. The repository allows the easy download of larger data files required to run SNEEP and contains 4 files, explained in more detail in the following. 

## The dbSNP database: 

To identify TFs effected more often in the given input SNP set than expected, SNEEP can perform a statistical assessment to compare the result against proper random controls. To do so, the pipeline randomly samples SNPs from the dbSNP database and rerun the analysis on these SNPs. 
In order to sample the SNPs in a fast and efficient manner, we provide a file containing the SNPs of the dbSNP database.  The file is a slightly modified version of the [public available one]( https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/) (file GCF_000001405.38). In detail, we 

-	removed all SNPs overlapping with a protein-coding region (annotation of the [human genome (GRCh38), version 36 (Ensembl 102)]( https://www.gencodegenes.org/human/release_36.html)),
-	removed all information not important for SNEEP,
-	removed mutations longer than 1 bp,
-	and sorted SNPs according to their MAF distribution in ascending order. 


## Epigenetic interactions: 

We provide three files containing epigenetic interactions associated to target genes:

-	interactionsREMs.txt provides regulatory elements (REMs) linked to their target genes. The data was derived with the STITCHIT algorithm, which is a peak-calling free approach to identify gene-specific REMs by analyzing epigenetic signal of diverse human cell-types with regard to gene expression of a certain gene. For more information, you can also have a look at our public [EpiRegio database](https://epiregio.de) holding all REMs stored in the interactionsREMs.txt file. 
-	interactionsREM_PRO.txt: Additional to the REMs the promoters (+/- 500 bp around TSS) of the genes are included as regions linked to their target genes. 
-	interactionsREMs_PRO_HiC.txt: We further added HiC regions linked to target genes via the ABC algorithm on human heart data from a [recent published paper from Anene-Nzelu *et al.*]( https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.046040?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed).


Please download the zenodo repository to the SNEEP (SNEEP/) directory and gunzip the files. 


It is also possible to use your own epigenetic interactions file or extend on of ours with for instance cell type specific data. Please stick to our tab-separated format: 

-	chr of the linked region
-	start of the linked region (0-based)
-	end of the linked region (0-based)
-	target gene (ensembl ID, may you also need to add the ensembl ID and the corresponding gene name to the file ensemblID_GeneName.txt)
-	unique identifier of the interaction region not longer than 10 letters/digits (e.g., PRO0000001, HiC0000234, … ), 
-	7 tab-separated dots (or additional information which you wish to keep -> displayed in the result.txt file but not in the summary pdf). 

# Usage 

The following 3 files are required as minimal input to run SNEEP:

1)	a file containing the TF motifs in TRANSFAC format,
2)	a bed-like SNP file,
3)	a reference genome file in fasta format.

We provide the human TF motifs from the JASPAR database (version 2020) in the required format in the examples directory (SNEEP/examples/ JASPAR2020_HUMAN_transfac_P0.txt).  However, every set of TF motifs can be used instead, for instance from another database (e.g. HOCOMOCO or Kellis) or a different species. 
The required bed-like SNP file is a tab-separated file containing the following entries: 
-	chr,
-	start position (0-based),
-	end position (0-based),
-	var1 (e.g. effector allele or alternative allele) ,
-	var2 (e.g. wild type allele),
-	rsID if known, otherwise - ,
-	minor allele frequency (MAF) if known, otherwise -1. The minor allele frequency is important to provide if you want to assess the result of SNEEP against random controls. Then, SNEEP samples SNPs based on the MAF distribution of the input SNPs. 

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

SNEEP can only handle mutations effecting a single base pair (no deletions or insertions). Deletions and insertions are identified by the pipeline and ignored. Also duplicated entries are only considered once.

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

Per default the result is stored in the directory ‘SNEEP_output’. For more details about the result files, see Section [Detailed explanation of the output files](#Detailed-explanation-of-the-output-files)

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


For an realistic example we consider SNPs associated to myocardial infarction (downloaded from the [GWAS catalog](https://www.ebi.ac.uk/gwas/efotraits/EFO_0000612)) and the corresponding proxy SNPs (determined with [SNiPA](https://snipa.helmholtz-muenchen.de/snipa3/index.php?task=proxy_search), R2 value >= 0.8). The following section provides example runs with different combination of optional input parameters. The example data is located in the directory SNEEP/example/. The default parameters (SNP-file, motif file and human genome file) are the once we already used in the minimal example. Make sure you are located in the SNEEP main folder (SNEEP/).

Notice, that the optional parameters can be combined in any way, not only in the ways shown in the following examples.

## Example 1: Consider only TFs expressed in the cell type or tissue of interest

To do so, we need to set the optional parameters -t, -d and -e. For our example the file containing the expression values (flag -t) are provided in the example directory and derived from cardio myocytes. Additionally, we provide the file containing the mapping between ensembl ID (used in the expression value file) and the names of the TFs specified in the motif file. As flag d, we use a rather less stringent expression value threshold of 0.5. In general, you can choose any value which is most suitable for you.

So, the resulting command is: 

```
./src/differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_expression/ -t examples/RNA-seq_humanLV_hiPSC-CM.txt -e examples/ensemblID_geneName_TFs.txt -d 0.5 examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path-to-genome-file>
```

Note, that we specified the output directory with the -o flag as examples/SNEEP_output_expression/. 

## Example 2: Add open chromatin regions of the cell type of interest

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
./src/differentialBindingAffinity_multipleSNPs  -o examples/SNEEP_output_open_chromatin/ -f <path-to-ENCODE-data/ENCFF199VHV.bed> examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path-to-genome-file>
```

## Example 3: Associate the SNPs, which significantly affect the binding behavior of a TF to their target genes

To associate the target genes, we need to specify a file holding information about epigenetic interactions (flag -e). We provide this data via a Zenodo repository, which contains three different epigenetic interaction files (for more detail explanation see …). For our example the most suitable one is the file interactionsREM_PRO_HiC.txt. The HiC data is retrieved from whole human heart, so we can benefit from the interactions for our example analysis. Additionally, the file ensemblID_GeneName.txt containing the ensembl ID to gene name mapping for all genes listed in the epigenetic interaction file is required (flag -g).
 

```
./src/differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_REM_PRO_HiC/ -r interactionsREM_PRO_HiC.txt -g ensemblID_GeneName.txt  examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path-to-genome-file>
```

## Example 4: Compute a proper random background control and highlight the biologically interesting results 

To perform a random background sampling the optional parameters -j, -k, -l and -q need to be specified. We recommend to sample 100 background rounds, meaning set -j to 100. The random SNPs are sampled from the dbSNP database. We provide the corresponding file in the Zenodo repository (for more information, see …), which is used to specify the flag -k. To allow reproducible results, we ask the user to set a random seed via the -l flag. Please use varying random seed for runs with different input SNPs. The flag -q is used to speed up the background sampling by exclude TFs, which did not have or did have less significant differential binding affinities on the input SNPs. Per default -q is set not 0, meaning only TFs with at least 1 significant change in the binding affinity are considered in the background sampling. 
Further we recommend running SNEEP in the parallel mode by specifying the number of threads via the -n flag. 

A possible SNEEP run with background sampling can look as following: 

```
./src/differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_background_sampling/ -n 20 -j 100 -k dbSNPs_sorted.txt -l 2 -q 0 -r interactionsREM_PRO_HiC.txt -g ensemblID_GeneName.txt  examples/JASPAR2020_HUMAN_transfac_P0.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path-to-genome-file>
```

Note, that we also associated the SNPs to their potential target genes (as shown in example 3).

To visualize the biologically interesting results, we generate a summary file via RMarkdown. 
In this summary information about TFs, which are associated to a change in the binding affinity caused by a SNP more often than expected, are presented. 
To generate the summary file, we provide the Rscript builtPDF.R. As input the 

-	absolut path to the SNEEP_output directory,
-	number of sampled background rounds,
-	number of significant changes in the binding affinity a TF must exceed to be considered in the summary file,
-	and the path the source folder (SNEEP/src/) 

must be given.
So for our current example, we run the following:

```
Rscript  src/builtPDF.R <absolute-path>/SNEEP/examples/SNEEP_output_background_sampling/ 100 4  src/
```

The resulting summary pdf is called summaryReport.pdf, and can be found in the SNEEP output directory, hence in our case in the  examples/SNEEP_output_background_sampling/ directory.

# Detailed explanation of the output files 
The output files of a SNEEP run can either be found in the default output directory (SNEEP_output/) or in the user-defined one. In the following the most important output files are explained in more detail: 

## Main result file (result.txt) 

The file lists all SNPs which cause a significant change of the binding affinity for a TF. A single SNP can affect multiple TFs, which results in multiple lines in the result file. The first 14 entries are the following: 

-	*SNP_position* chr:start-end (0-based),
-	*var1* e.g., effector allele or alternative allele,
-	*var2* e.g., wild type allele,
-	*rsID* 
-	*MAF* 
-	*peakPosition* if additional footprint/open chromatin region file is used (flag -f) the position of the overlapping region is given, otherwise . (meaning the option was not used).
-	*TF* name of the affected TF 
-	*TF-binding_position* genomic region, which is bound by the affected TF 
-	*strand* strand to which the TF binds (r) -> reverse strand, (f) -> forward
-	*effectedPositionInMotif* position within the TF binding motif which is affected by the SNP
-	*pvalue_BindAff_var1* p-value of the binding affinity for var1 
-	*pvalue_BindAff_var2* p-value of the binding affinity for var2
-	*log_pvalueBindAffVar1_pvalueBindAffVar2* logarithmic ratio of the p-value of the binding affinity of var1 and the p-value of the binding affinity of var2
-	*pvalue_DiffBindAff* corresponding p-value for the previous entry (log_pvalueBindAffVar1_pvalueBindAffVar2)

If the current SNPs overlaps with an epigenetic interaction, the following entries, provide more information about the interaction regions. Otherwise, the remaining entries are filled with ‘.’.

-	*REM_positions* genomic region of the epigentic interaction 
-	*ensemblIDs* ensemble ID of the associated gene
-	*geneNames* gene name of the associated gene
-	*REMIds* unique identifier of the epigenetic interactions (if it is an identifier for a REM it can be used to search for the region in our EpiRegio database) 

The last 8 entries contain information specific to EpiRegio. For more information, please have a look at epiregio.de.  


## Info file (info.txt) 

The file info.txt holds the input parameters used for the SNEEP run. 

## Summary report (summaryReport.pdf)

The file summaryReport.pdf presents our summary analysis and is only generated if the background sampling was performed. All information shown in the file are taken from result.txt, info.txt and TF_count.txt file. The TF_count.txt file contains for each TF how often this TF was significantly affected by a SNP for the input data and each random sampled round. All the tables and figures illustrated in the summary report are explained in the file itself. 

## Result of the random background sampling

The results for the random background sampling can be found in the directory sampling and contains per round the randomly sampled SNPs (random_SNPs_<round>.txt) and the randomResult.txt file (similar to result.txt). 
