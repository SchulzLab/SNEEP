
=======================================
Optional parameters
=======================================

In the following, we outline all optional parameters that can be used to run SNEEP in more detail. 

Flag -o: Specify an output folder
===================================
  
As a default, the result of our pipeline is stored within the folder SNEEP_output/.  Using the flag -o, a user-defined path for the output folder can be given. Note that if you set this flag, it must be given as first argument. The (potential) content of the output folder is automatically overwritten when running SNEEP.

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
To identify TFs which are more often affected by the given data than one would expected on random data, SNEEP can perform a statistical assessment to compare the result against proper random controls. To do so, the pipeline randomly samples SNPs from the `dbSNP database <https://www.ncbi.nlm.nih.gov/snp/>`_ and rerun the analysis on these SNPs. 
In order to sample the SNPs in a fast and efficient manner, we provide a file (in our `Zenodo repository <https://zenodo.org/record/4892591>`_ containing the SNPs of the dbSNP database.  The file is a slightly modified version of the `public available one <ttps://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/>`_ (file GCF_000001405.38). In detail, we 

-	removed all SNPs overlapping with a protein-coding region (annotation of the `human genome (GRCh38), version 36 (Ensembl 102) <https://www.gencodegenes.org/human/release_36.html>`_), (TODO: remove this sentence when zenodo dir is updated!)
-	removed all information not important for SNEEP,
-	removed mutations longer than 1 bp,
-	and sorted SNPs according to their MAF distribution in ascending order. 

Flag -r and -g: Epigenetic interactions
=============================================== 
We provide three files (in our `Zenodo repository <https://zenodo.org/record/4892591>`_) containing epigenetic interactions associated to target genes:

-	interactionsREMs.txt provides regulatory elements (REMs) linked to their target genes. The data was derived with the STITCHIT algorithm, which is a peak-calling free approach to identify gene-specific REMs by analyzing epigenetic signal of diverse human cell types with regard to gene expression of a certain gene. For more information, you can also have a look at our public `EpiRegio database <https://epiregio.de>`_ holding all REMs stored in the interactionsREMs.txt file. 
-	interactionsREM_PRO.txt: Additional to the REMs the promoters (+/- 500 bp around TSS) of the genes are included as regions linked to their target genes. 
-	interactionsREMs_PRO_HiC.txt: This file further includes enhancer-gene links predicted with the ABC algorithm on human heart data from a `published paper from Anene-Nzelu *et al.* <https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.046040?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed>`_.

It is also possible to use your own epigenetic interactions file (for instance generated with STARE's gABC score computation) or extend one of ours with for instance cell type specific data. Please stick to our tab-separated format: 
  
-	chr of the linked region
-	start of the linked region (0-based)
-	end of the linked region (0-based)
-	target gene (ensembl ID, you may also need to add the ensembl ID and the corresponding gene name to the file ensemblID_GeneName.txt)
-	unique identifier of the interaction region not longer than 10 letters/digits (e.g., PRO0000001, HiC0000234, … ), 
-	7 tab-separated dots (or additional information which you wish to keep -> displayed in the result.txt file but not in the summary pdf). 

Further, a file which provides a mapping between ensemblID to gene name must be given. This file comes along with our GitHub repository. 

  
Flag -a: Store D\ :sub:  `max`  values for all considered shifts
=====================================================
If this flag is set, for all shifts that exceed the TF binding affinity p-value threshold the resulting D-max value and the corresponding p-value is stored in <outputDir>/AllDiffBindAffinity.txt

Flag -f: Include open chromatin data
======================================

To consider only the SNPs which overlap with  cell type specific open chromatin data, a peak file in bed-format can be specified with this flag.

Flag -m: Get all D\ :sub:  `max`  values
===============================

If this flag is set all absolute maximal differential TF binding scores are printed (to the console) even if they do not exceed the specified p-value threshold. This flag is useful for estimating the scale parameter

Flag -t, -d and -e: Active TFs of the cell type of interest
=============================================================
In order to only consider the TFs which are expressed in your analysed cell type or tissue our computational approach needs three additional information. A file containing the expression value per TF (-t),  a threshold to decide which TFs are active and a mapping between the ensemblID and the TF name, The last file is provided in our GitHub repository for the TF set used within our analyses. 

Flag -j: Number of sampled background SNP sets
=================================================

With this flag the number of background rounds can be specified. Default: -j 0.

Flag -l: Reproducible results for random background analysis
==============================================================
In order to reproduce the result of the random background analysis we recommend to specific a seed variable. Default: -l 1. 

Flag -q:  TF count
=====================
This flags allows to exclude TFs from the baclground sampling which do not exceed a TF count. Default: -q 0
