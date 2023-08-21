
=======================================
Optional parameters
=======================================

In the following we outline all optional parameters which can be used to run SNEEP in more detail. 

Flag -o: Specify an output folder
===================================
  
As default the result of our pipeline is stored within the folder SNEEP_output/.  Using the flag -o a user-defined path for the output folder can be given. Notice that if you set this flag, it must be given as first argument. The (potential) content of the output folder is automatic whenn running SNEEP.

Flag -n: Number of threads
==========================
  
To speed up the p-value computation of the binding affinity and the random background analyses the parameter n can be specified. As default the number of threads is set to 1. 

Flag -p: p-value threshold for TF binding affinity
===================================================
  
This flag specifies the p-value threshold for the TF binding affinity. For a TF the binding affinity is computed for all possible shifts that overlap with the SNP. If a shift exceeds the p-value threshold, a absolute maximal differential TF binding score is computed. Per default this parameter is set to 0.05.
  
Flag -b and flag -x: Base frequency for TF binding affinity computation
=========================================================================
The approach to derive the p-value for the TF binding score allows to include the base frequency of the bases and the transition frequency between to bases of the considered sequences. We computed both frequencies for the human genome and provide them within our GitHub repository (necessaryInputFiles/frequency.txt and necessaryInputFile/transition_matrix.txt). If these frequencies are not specified, we assume as default 0.25.


Flag -c: p-value threshold for absolute maximal differential TF binding score
===============================================================================
The p-value threshold for D-max is per default set to 0.01. In our benchmarking analyses outlined in our `paper <sneep paper>`_, best results were observed for p-value thresholds between 0.01 and 0.001.

Flag -k: dbSNP database (dbSNPs_sorted.txt.gz)
=============================================== 
To identify TFs which are more often affected than expected by chance in the given input SNP set, SNEEP can perform a statistical assessment to compare the result against proper random controls. To do so, the pipeline randomly samples SNPs from the `dbSNP database <??>`_ and rerun the analysis on these SNPs. 
In order to sample the SNPs in a fast and efficient manner, we provide a file (in our `Zenodo repository <https://zenodo.org/record/4892591>`_ containing the SNPs of the dbSNP database.  The file is a slightly modified version of the `public available one <ttps://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/>`_ (file GCF_000001405.38). In detail, we 

-	removed all SNPs overlapping with a protein-coding region (annotation of the `human genome (GRCh38), version 36 (Ensembl 102) <https://www.gencodegenes.org/human/release_36.html>`_), (TODO: remove this sentence when zenodo dir is updated!)
-	removed all information not important for SNEEP,
-	removed mutations longer than 1 bp,
-	and sorted SNPs according to their MAF distribution in ascending order. 


Flag -r and -g: Epigenetic interactions
=============================================== 
We provide three files (in our `Zenodo repository <??>`_) containing epigenetic interactions associated to target genes:

-	interactionsREMs.txt provides regulatory elements (REMs) linked to their target genes. The data was derived with the STITCHIT algorithm, which is a peak-calling free approach to identify gene-specific REMs by analyzing epigenetic signal of diverse human cell types with regard to gene expression of a certain gene. For more information, you can also have a look at our public `EpiRegio database <https://epiregio.de>`_ holding all REMs stored in the interactionsREMs.txt file. 
-	interactionsREM_PRO.txt: Additional to the REMs the promoters (+/- 500 bp around TSS) of the genes are included as regions linked to their target genes. 
-	interactionsREMs_PRO_HiC.txt: This file further includes enhancer-gene links predicted with the ABC algorithm on human heart data from a `published paper from Anene-Nzelu *et al.* <https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.046040?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed>`_.

It is also possible to use your own epigenetic interactions file or extend on of ours with for instance cell type specific data. Please stick to our tab-separated format: 
  
-	chr of the linked region
-	start of the linked region (0-based)
-	end of the linked region (0-based)
-	target gene (ensembl ID, you may also need to add the ensembl ID and the corresponding gene name to the file ensemblID_GeneName.txt)
-	unique identifier of the interaction region not longer than 10 letters/digits (e.g., PRO0000001, HiC0000234, … ), 
-	7 tab-separated dots (or additional information which you wish to keep -> displayed in the result.txt file but not in the summary pdf). 

Further, a file which provides a mapping between ensemblID to gene name must be given. This file comes along with our GitHub repository. 
  
Flag -s: Estimated scale parameters for the TFs used
=====================================================

Our modified Laplace distribution is dependent on two parameters: n, which is two times the length of the TF model and b, which needs to be estimated. 
For the TF set we provide within our GitHub repository, we also estimated the scale parameter listed in XX. 
In case a customized TF motif set is used, one also needs to estimate the scale parameter for each TF. Therefore we provide a script XXX (TODO: provide more details here).
  
Flag -a: Store Dmax values for all considered shifts
=====================================================
If this flag is set, for all shifts that exceed the TF binding affinity p-value threshold the resulting D-max value and the corresponding p-value is stored in <outputDir>/AllDiffBindAffinity.txt

Flag -f: Include open chromatin data
======================================

To consider only the SNPs which overlap with  cell type specific open chromatin data, a peak file in bed-format can be specified with this flag.

Flag -m: Get all Dmax values
===============================

If this flag is set all absolute maximal differential TF binding scores are printed (to the console) even if they do not exceed the specified p-value threshold. This flag is useful for estimating the scale parameter

Flag -t, -d and -e: Active TFs of the cell type of interest
=============================================================
In order to only consider the TFs which are expressed in your analysed cell type or tissue our computational approach needs three additional information. A file containing the expression value per TF (-t),  a threshold to decide which TFs are active and a mapping between the ensemblID and the TF name, The last file is provided in our GitHub repository for the TF set used within our analyses. 

Flag -j: Number of sampled background SNP sets
=================================================

With this flag the number of background rounds can be specified, default 0.

Flag -l: Reproducible results for random background analysis
==============================================================
In order to reproduce the result of the random background analysis we recommend to specific a seed variable. Default is 1. 

Flag -q:  TF count
=====================
This flags allows to exclude TFs from the baclground sampling which do not exceed a TF count (default 0).
