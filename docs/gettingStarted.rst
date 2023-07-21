===============
Getting started
===============

SNEEP is fast method to identify regulatory non-coding SNPs (rSNPs) that modify the binding sites of Transcription Factors (TFs) for large collections of SNPs provided by the user. SNEEP combines our statistical approach introduced in our preprint `A statistical approach to identify regulatory DNA variants <https://www.biorxiv.org/content/10.1101/2023.01.31.526404v1>`_ with epigenetic data and additional functionalities.

A graphical summary of SNEEP is shown below:

TODO: add figure

Installation 
==============
We provide a bioconda package to install the main functionality of our approach. Therefore an installation of  Bioconda `here <https://bioconda.github.io/>`_ is requiered. 

TODO: add command

If you want to directly install SNEEP, please clone our `GitHub repository <https://github.com/SchulzLab/SNEEP/>`_ and make sure that the following software is available on your machine: 

- c++11 
- g++ (9.3.0)
- python3.x
- bedtools (v2.27.1)
- openmp

To build SNEEP, run the following commands: 

.. code-block:: console

  cd SNEEP/src/
  make


Please add the path to our software (SNEEP/src) to our PATH environment (otherwise internally called scrips might not be found)

We tested the code and the Makefile only on a linux machine. 

Testing your installation 
==========================

We provide a test script to verify if your installation worked. To download the test data and scripts, please clone the latest version of our GitHub repository

.. code-block:: console

  clone git@github.com:SchulzLab/SNEEP.git

, download our `Zenodo repository <https://doi.org/10.5281/zenodo.4892591>`_ and unzip the files. 

Additional a reference genome in fasta format is requiered. The different chromosomes within the file must be named as chr1, chr2. An example is shown below:

.. example::
  >chr1
  ATCGGGTCA…
  >chr2
  TTTGAGACCAT…

For the provided examples in the following, please use genome version hg38.

To run our tests, please redirect into the SNEEP folder downloaded from GitHub and perform 

.. code-block:: console

  bash runTests.sh <pathToGenome> <pathTodbSNP>  <pathToInteractions>

where <pathTodbSNP> is the path to the dbSNPs_sorted.txt downloaded from Zenodo repository and  <pathToInteractions> needs to be the path to one of interactiosn file e.g. interactionsREM_PRO_HiC.txt

Basic usage
============

The following 3 files are required as minimal input to run SNEEP:

1)	a file containing the TF motifs in TRANSFAC format, 
2)	a bed-like SNP file,
3)	a reference genome file in fasta format.

Minimal example
---------------

To try SNEEP with the minimal required input, make sure you are in the SNEEP folder and run: 

.. code-block:: console

  ./src/differentialBindingAffinity_multipleSNPs examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed  <path-to-genome-file> 

Per default the result is stored in the directory ‘SNEEP_output’. The file ‘result.txt’ in the SNEEP output directory contains the predicted rSNPs. For more details about the result files, see Section [Detailed explanation of the output files](Detailed-explanation-of-the-output-files). The run takes about 2 to 3 minutes. 


Detailed description of the requiered input files
----------------------------------------------------

(1) We provide human TF motifs from the JASPAR database (version 2022), HOCOMOCO and  Kellis ENCODE database in the required format in the examples directory (SNEEP/examples/ combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt).  

(2) The required bed-like SNP file is a tab-separated file containing the following entries: 
-	chr,
-	start position (0-based),
-	end position (0-based),
-	var1 (e.g. effector allele or alternative allele) ,
-	var2 (e.g. wild type allele),
-	rsID if known, otherwise - ,
-	minor allele frequency (MAF) if known, otherwise -1. The minor allele frequency is important to provide if you want to assess the result of SNEEP against random controls. Then, SNEEP samples SNPs based on the MAF distribution of the input SNPs. 

An example of a properly formated SNP files can be found below: 

.. example::
  chr1    109274569       109274570       G       A       rs7528419       0.2009
  chr1    109275907       109275908       C       T       rs646776        0.2384
  chr1    154424939       154424940       G       T       rs12118721      1e-07
  chr1    154424939       154424940       G       T       -      0.3
  chr12   111569951       111569952       G       C       rs653178        -1


If you want to consider a SNP, which has multiple alternative alleles like for instance [rs11206510](https://www.ncbi.nlm.nih.gov/snp/rs11206510) (T -> A,C,G) , please add one line per alternative allele in the bed-like SNP file. An example is shown below: 

.. example::
  chr1    55030365        55030366        A       T       rs11206510      0.1018
  chr1    55030365        55030366        C       T       rs11206510      0.1018
  chr1    55030365        55030366        G       T       rs11206510      0.1018


SNEEP can only handle mutations effecting a single base pair (no deletions or insertions). Deletions and insertions are identified by the pipeline and ignored. Also duplicated entries are only considered once.

(3) In the reference genome file, the different chromosome must be named as chr1, chr2 etc. resulting in the following format: 

.. example::
  >chr1
  ATCGGGTCA…
  >chr2
  TTTGAGACCAT…

For the provided examples in the following, please use genome version hg38.






