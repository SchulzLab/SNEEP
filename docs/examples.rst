=====================
Examples of how to run SNEEP
=====================

In this section, we provide additional details on how to use our SNEEP pipeline for different applications. 

Optional input parameters
=========================

To display all possible input parameters, type 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs  -h

which results in 

.. code-block:: console

  Call program with ./src/differentialBindingAffinity_multipleSNPs
  optinal parameters:
  -o outputDir (default SNEEP_output/, if you want to specify it, it must be done as first argument)
  -n number threads (default 1)
  -p pvalue for motif hits (default 0.5)
  -c pvalue indicating differential binding (default 0.01)
  -b base frequency for PFMs -> PWMs (/necessaryInputFiles/frequency.txt)
  -a if flag is set,  all computed differential binding scores are stored in <outputDir>/AllDiffBindAffinity.txt
  -f additional footprint/open chromatin region file in bed file format
  -m if flag is set, the  maximal differential binding score per SNP is printed
  -t file where expression values of TFs are stored (e.g RNA-seq in a tab-seperated format e.g. ensemblID	expression-value)
  -d threshold TF activity (must be given if -t is given)
  -e tab-seperated file containing ensemblID to gene name mapping of the TFs (must be given if -t is given)
  -r bed-like file with epigenetic interactions
  -g path to file containing ensemblID to gene name mapping, must be given if -r is given (,-seperated)(mapping for all genes within EpiRegio)
  -j rounds sampled background (default 0)
  -k path to sorted dbSNP file (if our provided file is used only SNPs in coding regions are considered)
  -l start seed (default 1)
  -q minimal TF count that needs to be exceeded to be considered in random sampling (default 0)
  -u gene background analysis is performed (defaul false), -j must be set 
  -v perform TF enrichment analysis (default  false), -j must be set
  -x transition matrix for binding affinity p-value, (default all transitions are equally likely) (necessaryInputFiles/transitionMatrix.txt)-h help
  transfac PFM file,  bed-like SNP file, path to genome file (fasta format) and scale file (see necessaryInputFiles/estimatedScalesPerMotif_1.9.txt for human data) must be given
  help function end


Examples of realistic applications
===================================

For a realistic example, we considered SNPs associated with myocardial infarction (downloaded from the `GWAS catalog <https://www.ebi.ac.uk/gwas/efotraits/EFO_0000612>`_) and the corresponding proxy SNPs (determined with `SNiPA <https://snipa.helmholtz-muenchen.de/snipa3/index.php?task=proxy_search>`_, R2 value >= 0.8). The following section provides example runs with different combinations of the optional input parameters. The example data are located in the directory SNEEP/example/. The default parameters (SNP-file, motif file, human genome file and scale file) are the once we already used in the minimal example. Make sure you are located in the SNEEP main folder (SNEEP/).

Note, that the optional parameters can be combined in many additional ways, as presented in the following examples.

Example 1: Consider only TFs expressed in the cell type or tissue of interest
------------------------------------------------------------------------------

To do so, we need to set the optional parameters -t, -d and -e. For our example, the file containing the expression values (flag -t) are provided in the example directory and derived from cardiomyocytes. Additionally, we provided the a file containing the mapping between the Ensembl ID (used in the expression value file) and the names of the TFs specified in the motif file (flag -e). For flag -d, we used a rather less stringent expression value threshold of 0.5. In general, you can choose any value, that is most suitable for you (also depending on the normalization that you did for the gene expression data).
Additionally, we specify a p-value threshold for the TF binding score (-p), a p-value threshold D\ :sub:`max` as 0.001 (-c), the number of threads to use as 10 (-n) and two files to refine for the binding affinity p-value computation (-b and -x).

Therefore, the resulting command is: 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_expression/ -t examples/RNA-seq_humanLV_hiPSC-CM.txt -e examples/TF_ensemblID_name_human_JASPAR2022_GRCh38p13.txt -d 0.5 -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt  -p 0.5 -c 0.001 -n 10 examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <pathToGenome> necessaryInputFiles/estimatedScalesPerMotif_1.9.txt
 
Note, that we specified the output directory with the -o flag as examples/SNEEP_output_expression/. 

Example 2: Consider only SNPs in the open chromatin of the cell type of interest
---------------------------------------------------------------------------------------

If open chromatin data of your cell type of interest are available, it is possible to integrate these data and automatically exclude SNPs from the analysis in closed, most likely inactive chromatin regions. 
Therefore, a bed-file holding the open chromatin regions can be specified using the flag -f. 

For our example, we used an ATAC-seq data from the  human heart right ventricle from ENCODE. 

To download the data run: 

.. code-block:: console

  wget 'https://www.encodeproject.org/files/ENCFF199VHV/@@download/ENCFF199VHV.bed.gz'

Next, unzip the file via gzip.

The resulting SNEEP call is 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs  -o examples/SNEEP_output_open_chromatin/  -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt -f ENCFF199VHV.bed  -p 0.5 -c 0.001 -n 10  examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt examples/SNPs_EFO_0000612_myocardial_infarction.bed <pathToGenome> necessaryInputFiles/estimatedScalesPerMotif_1.9.txt
  
Example 3: Associating regulatory SNPs with their target genes
------------------------------------------------------------------------------------------------------------

To associate the target genes, we need to specify a file that contains enhancer-gene interactions (flag -r). We provide these data via our Zenodo repository, that contains three different epigenetic interaction files (for more detail explanation click `here <https://sneep.readthedocs.io/en/latest/detailsInput.html#flag-r-and-g-epigenetic-interactions>`_). For our example, the most suitable interactions are stored in the file interactionsREM_PRO_HiC.txt. The HiC data was retrieved from the whole human heart, so we can benefit from these interactions in our example analysis. Please specify the path to this file in the following command. Additionally, the file ensemblID_GeneName.txt containing the Ensembl ID to gene name mapping for all genes listed in the epigenetic interaction file is needed (flag -g).
 
.. code-block:: console

  differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_REM_PRO_HiC/   -r <pathToInteractions> -g examples/ensemblID_GeneName.txt  -p 0.5 -c 0.001 -n 10 -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path_to_genome> necessaryInputFiles/estimatedScalesPerMotif_1.9.txt 

Example 4: Compute a proper random background control and highlight cell type-specific TFs
---------------------------------------------------------------------------------------------

To perform a random background sampling, the optional parameters -j, -k, -l and -q need to be specified. We recommend to sample at least 100 background rounds, meaning that we set -j to 100. However, in our applications, we usually set -j to 500 or 1.000. The random SNPs were sampled from the dbSNP database. We provide the corresponding file in the Zenodo repository (unzipped file: dbSNPs_sorted.txt; additional information is found in the Section *Optional parameters*), which is used to specify the flag -k. To allow reproducible results, we ask the user to set a random seed via the -l flag. Please use varying random seeds for runs with different input SNPs. The flag -q is used to speed up the background sampling by excluding TFs, that did not have or did have less significant differential binding affinities on the input SNPs. Per default -q is set not 0, meaning that only TFs with at least 1 significant change in the binding score were considered in the background sampling. 
Furthermore, we recommend running SNEEP in the parallel mode by specifying the number of threads via the -n flag. 

A possible SNEEP run with background sampling can look as follows: 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_background_sampling/ -p 0.5 -c 0.001 -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt  -n 20 -j 100 -k <pathTodbSNP> -l 2 -q 0 -r <pathToInteractions> -g ensemblID_GeneName.txt  examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <path_to_genome> necessaryInputFiles/estimatedScalesPerMotif_1.9.txt 
