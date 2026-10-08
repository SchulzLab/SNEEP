
=======================================
Optional parameters
=======================================

In the following, we outline all optional parameters that can be used to run SNEEP in more detail. 

Flag -o: Specify an output folder
===================================
  
As a default, the result of our pipeline is stored within the folder SNEEP_output/.  Using the flag -o, a user-defined path for the output folder can be given (at any position among the optional parameters, with or without a final /). 

The output folder must either be empty (or not exist yet) or contain the output of a former SNEEP run (recognized by the file info.txt). In the latter case, the former output is deleted. If the folder contains other files, SNEEP stops with an error message and does not change the folder. This prevents deleting other data by mistake, e.g., with -o . or a mistyped path.

Flag -n: Number of threads
==========================
  
To speed up the p-value computation of the binding affinity and the random background analyses, the parameter n can be specified. Default: -n 1. 

Flag -p: p-value threshold for the TF binding score
===================================================
  
This flag specifies the p-value threshold for the TF binding score. For a TF, the binding score is computed for all possible shifts that overlap with the SNP. If a shift exceeded the p-value threshold, an absolute maximal differential TF binding score was computed. We recommend using a moderate p-value threshold of 0.5. Default: -p 0.5.
  
Flag -b and flag -x: base frequency for TF binding score computation
=========================================================================
The approach to derive the p-value for the TF binding score allows us to include the base frequency of the bases and the transition frequency between to bases of the considered sequences. We computed both frequencies for the human genome and provided them within our GitHub repository (necessaryInputFiles/frequency.txt and necessaryInputFile/transition_matrix.txt). If these frequencies are not specified, we assume a default of 0.25.

Flag -c: p-value threshold for the absolute maximal differential TF binding score
===============================================================================
The p-value threshold for D\ :sub: `max` is set to 0.01 by default. In our benchmarking analyses outlined in our `paper <sneep paper>`_, the best results were observed for p-value thresholds between 0.01 and 0.001.

Flag -k: dbSNP database (dbSNPs_sorted.txt.gz)
=============================================== 
To identify TFs that are more often affected by the given data than one would expect from random data, SNEEP can perform a statistical assessment to compare the results against proper random controls. To do so, the pipeline randomly samples SNPs from the `dbSNP database <https://www.ncbi.nlm.nih.gov/snp/>`_ and rerun the analysis on these SNPs. The random SNPs are matched to the minor allele frequency (MAF) distribution of the input SNPs (bins of width 0.01) and, optionally, to their GC content (flag -s).
To sample the SNPs in a fast and efficient manner, we provided a file (in our `Zenodo repository <https://zenodo.org/record/4892591>`_ containing the SNPs of the dbSNP database.  The file is a slightly modified version of the `publicly available one <https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/>`_ (file GCF_000001405.38). In detail, 

-	all information not important for SNEEP were removed,
-	mutations longer than 1 bp were removed,
-	and we sorted SNPs according to their MAF distribution in ascending order. 

In older versions of the file (dbSNP build 154), all SNPs overlapping with a protein-coding region were removed (annotation of the `human genome (GRCh38), version 36 (Ensembl 102) <https://www.gencodegenes.org/human/release_36.html>`_). The newer version (dbSNP build 157) keeps these SNPs and additionally contains the GC content around each SNP, which is needed for flag -s.

.. TODO: adapt this section when the file of dbSNP build 157 (with GC content) is available on Zenodo.

The file is tab-separated, sorted by the MAF (first column), and contains one line per SNP:

.. code-block:: console

  MAF  chr  start  end  ref  alt  rsID  MAF  GC

where MAF is -1 if no allele frequency is given in dbSNP, alt can hold several alleles separated by commas, and GC (only in the newer version) is the GC content in a window of +- 30 bp around the SNP: (#C + #G) / (#A + #C + #G + #T), lower case bases are counted and N is excluded (-1 if the window contains no A, C, G or T). The file is created with src/getSNPInfo.cpp from the dbSNP VCF file and the genome (getSNPInfo <dbSNP VCF> <outputDir> <genome.fa> [numThreads]); since this takes very long, we recommend to use the provided file.

Flag -r and -g: Epigenetic interactions
=============================================== 
We provide three files (in our `Zenodo repository <https://zenodo.org/record/4892591>`_) containing epigenetic interactions associated to target genes:

-	interactionsREMs.txt provides regulatory elements (REMs) linked to their target genes. The data were derived with the STITCHIT algorithm, which is a peak-calling free approach for idenitifying gene-specific REMs by analyzing the epigenetic signals of diverse human cell types with regard to the gene expression of a certain gene. For more information, you can also have a look at our public `EpiRegio database <https://epiregio.de>`_ holding all REMs stored in the interactionsREMs.txt file. 
-	interactionsREM_PRO.txt: Additional in addition to the REMs, the promoters (+/- 500 bp around the TSS) of the genes are included as regions linked to their target genes. 
-	interactionsREMs_PRO_HiC.txt: This file further includes enhancer-gene links predicted with the ABC algorithm on human heart data from a `published paper from Anene-Nzelu *et al.* <https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.046040?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed>`_.

It is also possible to use your own epigenetic interactions file (for instance generated with STARE's gABC score computation) or extend one of ours with for instance cell type specific data. Please adhere to our tab-separated format: 
  
-	chr of the linked region
-	start of the linked region (0-based)
-	end of the linked region (0-based)
-	target gene (ensembl ID, you may also need to add the ensembl ID and the corresponding gene name to the file ensemblID_GeneName.txt)
-	unique identifier of the interaction region no longer than 10 letters/digits (e.g., PRO0000001, HiC0000234, … ), 
-	7 tab-separated dots (or additional information which you wish to keep -> displayed in the result.txt file but not in the summary pdf). 

Furthermore, a file which provides a mapping between ensemblID to gene name must be given. This file comes along with our GitHub repository. 

  
Flag -a: Store D\ :sub:  `max`  values for all considered shifts
=====================================================
If this flag is set, for all shifts that exceed the TF binding score p-value threshold, the resulting D-max value and the corresponding p-value are stored in <outputDir>/AllDiffBindAffinity.txt

Flag -f: Include open chromatin data
======================================

To consider only the SNPs that overlap with  cell type-specific open chromatin data, a peak file in bed-format can be specified with this flag.

Flag -m: Get all D\ :sub:  `max`  values
===============================

If this flag is set all absolute maximal differential TF binding scores are printed (to the console) even if they do not exceed the specified p-value threshold. This flag is useful for estimating the scale parameter

Flag -t, -d and -e: Active TFs of the cell type of interest
=============================================================
To consider only the TFs that are expressed in the analyzed cell type or tissue, our computational approach requires three pieces of information. A file tthat contains the expression value per TF (-t),  a threshold for deciding which TFs are active and a mapping between the ensemblID and the TF name. The last file is provided in our GitHub repository for the TF set used within our analyses. 

Flag -j: Number of sampled background SNP sets
=================================================

With this flag, the number of background rounds can be specified. Default: -j 0.

Flag -i: Use already sampled random SNPs
==========================================
Instead of sampling the random SNPs from the dbSNP file (-k), already sampled ones can be given, e.g., the directory sampling/ of a former SNEEP run. The directory must contain the files randomSNPs_0.txt, ..., randomSNPs_<j-1>.txt for the number of rounds given with -j; -k is then not needed. The directory is only read: all files of the background rounds (e.g., randomResult_<round>.txt) are written to <outputDir>/sampling/ as for sampled SNPs.

Flag -l: Reproducible results for random background analysis
==============================================================
To reproduce the results of the random background analysis, we recommend the use of a specific seed variable. Default: -l 1. 

Note that the same seed reproduces the same random SNPs only with the same compiler and C++ standard library (e.g., the random SNPs differ between macOS and Linux). The random number generator itself is the same everywhere, but the conversion of its numbers into positions (std::uniform_int_distribution) is implemented differently by different standard libraries and versions. The results for the input SNPs are not affected.

Flag -q:  TF count
=====================
This flag allows us to exclude TFs from the background sampling that do not exceed a TF count. Default: -q 0

Flag -s: Match the GC content in the background sampling
=========================================================
If -s is set to true (also accepted: 1), the random SNPs are matched not only to the MAF but also to the GC content of the input SNPs. Default: -s false.

- The GC content of an input SNP is computed from its sequence (<outputDir>/snpRegions.fa) in a window of +- 30 bp around the SNP, exactly as for the dbSNP file (see flag -k): (#C + #G) / (#A + #C + #G + #T), lower case bases are counted, N is excluded, and the reference base at the SNP position is included.
- The random SNPs are sampled per combination of MAF bin (width 0.01) and GC bin (width 0.05).
- The dbSNP file (flag -k) must contain the GC content (column 9, newer version of the file); otherwise SNEEP stops with an error message.
- The flag has no effect if the random SNPs are given (-i) or no background sampling is performed (-j 0); SNEEP then prints a warning.

Each background round always contains exactly as many SNPs as the input. If a bin cannot be filled as intended, SNEEP uses a fallback and writes a warning to the console and to info.txt (once per bin, not per round):

- more than half of the GC window of an input SNP is N: warning only, the GC content might not be meaningful,
- the GC window of an input SNP contains no A, C, G or T: this SNP is matched by MAF only,
- a MAF x GC bin contains fewer dbSNP SNPs than needed: each dbSNP SNP is taken once, the remaining ones are sampled with replacement,
- a GC bin is empty: the SNPs are sampled from the nearest non-empty GC bin within the same MAF bin (random choice if two bins are equally near); if the MAF bin has no SNP with GC content, by MAF only,
- a MAF bin of the input does not exist in the dbSNP file (e.g., MAF > 0.5): the nearest MAF bin is used.

With -s false, the random SNPs are identical to those of former SNEEP versions (same seed -l), except for the rare case of an empty MAF bin in the dbSNP file, where former versions assigned the following SNPs to the wrong bin.
