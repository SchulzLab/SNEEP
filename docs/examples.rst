=====================
Examples how to run SNEEP
=====================

In this section we provide more details how to use our SNEEP pipeline for different applications. 

Optional input parameters
=========================

To display all possible input parameters, type 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs  -h

which results in 

.. code-block:: console

  Call program with ./src/differentialBindingAffinity_multipleSNPs
  optinal parameters:
  -o outputDir (default SNEEP_output/, if you want to specific it, it must be done as first argument)
  -n number threads (default 1)
  -p pvalue for motif hits (default 0.05)
  -c pvalue differential binding (default 0.01)
  -b base frequency for PFMs -> PWMs ( /necessaryInputFiles/frequency.txt)
  -s file where the computed scales per motif are stored ( necessaryInputFiles/estimatedScalesPerMotif_1.9.txt) 
  -a if flag is set,  all computed differential bindinding affinities are stored in <outputDir>/AllDiffBindAffinity.txt
  -f additional footprint/open chromatin region file in bed file format
  -m if flag is set, the  maximal differential binding affinity per SNP is printed
  -t file where expression values of TFs are stored (e.g RNA-seq in a tab-seperated format e.g. ensemblID	expression-value)
  -d threshold TF activity (must be given if -t is given)
  -e tab-seperated file containing ensemblID to gene name mapping of the TFs (must be given if -t is given)
  -r bed-like file with epigenetic interactions
  -g path to file containing ensemblID to gene name mapping, must be given if -r is given (,-seperated)(mapping for all genes within EpiRegio)
  -j rounds sampled background (default 0)
  -k path to sorted dbSNP file (if our provided file is used only SNPs in coding regions are considered)
  -l start seed (default 1)
  -q minimal TF count which needs to be exceeded to be considered in random sampling (default 0)
  -u gene background analysis is performed (defaul false), -j must be set 
  -v perform TF enrichment  analysis (default  false), -j must be set
  -x transition matrix for binding affinity p-value, (default all transitions are equally likely) (necessaryInputFiles/transitionMatrix.txt)-h help
  transfac PFM file,  bed-like SNP file and path to genome file (fasta format)  must be given
  help function end


Examples of realistic applications
===================================

For a realistic example we consider SNPs associated to myocardial infarction (downloaded from the `GWAS catalog <https://www.ebi.ac.uk/gwas/efotraits/EFO_0000612>`_ and the corresponding proxy SNPs (determined with `SNiPA <https://snipa.helmholtz-muenchen.de/snipa3/index.php?task=proxy_search>`_, R2 value >= 0.8). The following section provides example runs with different combination of the optional input parameters. The example data is located in the directory SNEEP/example/. The default parameters (SNP-file, motif file and human genome file) are the once we already used in the minimal example. Make sure you are located in the SNEEP main folder (SNEEP/).

Notice, that the optional parameters can be combined in many more ways than presented in the following examples. We recommand to always set the -c flag to provide for each TF a specific scale parameter (see XX).

Example 1: Consider only TFs expressed in the cell type or tissue of interest
------------------------------------------------------------------------------

To do so, we need to set the optional parameters -t, -d and -e. For our example the file containing the expression values (flag -t) are provided in the example directory and derived from cardiomyocytes. Additionally, we provide the file containing the mapping between ensembl ID (used in the expression value file) and the names of the TFs specified in the motif file. As flag d, we use a rather less stringent expression value threshold of 0.5. In general, you can choose any value which is most suitable for you.
Additional, we specify the file holding the scale value per TF (-s), a p-value threshold of D_{max} as 0.001 (-c), the number of threads to use as 10 (-n) and two files to refinefor the binding affinity p-value computation (-b and -x).

So, the resulting command is: 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_expression/ -t examples/RNA-seq_humanLV_hiPSC-CM.txt -e examples/TF_ensemblID_name_human_JASPAR2022_GRCh38p13.txt -d 0.5  -s necessaryInputFiles/estimatedScalesPerMotif_1.9.txt -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt -c 0.001 -n 10 examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed <pathToGenome>
 
Note, that we specified the output directory with the -o flag as examples/SNEEP_output_expression/. 

Example 2: Consider only SNPs in open chromatin of the cell type ir tissue of interest
---------------------------------------------------------------------------------------

If open chromatin data of your cell type of interest is available, it is possible to integrate this data and automatically exclude SNPs from the analysis in closed, most likely inactive chromatin regions. 
Therefore, a bed-file holding the open chromatin regions can be specified using the flag -f. 

For our example, we want to use an ATAC-seq on human heart right ventricle from ENCODE. 

To download the data run: 

.. code-block:: console

  wget 'https://www.encodeproject.org/files/ENCFF199VHV/@@download/ENCFF199VHV.bed.gz'


Next unzip the file via gunzip.

The resulting SNEEP call is 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs  -o examples/SNEEP_output_open_chromatin/  -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_m    atrix.txt -f ENCFF199VHV.bed  -c 0.001 -n 10 -s necessaryInputFiles/estimatedScalesPerMotif_1.9.txt examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt examples/SNPs_EFO_0000612_myocardial_infarction.bed <pathToGenome>
  
Example 3: Associate regulatory SNPs to their target genes
------------------------------------------------------------------------------------------------------------

To associate the target genes, we need to specify a file that holds enhancer-gene interactions (flag -r). We provide this data via a Zenodo repository, which contains three different epigenetic interaction files (for more detail explanation see XX). For our example the most suitable one is the file interactionsREM_PRO_HiC.txt. The HiC data is retrieved from whole human heart, so we can benefit from the interactions for our example analysis. Please specify the path to this file in the following command. Additionally, the file ensemblID_GeneName.txt containing the ensembl ID to gene name mapping for all genes listed in the epigenetic interaction file is required (flag -g).
 
.. code-block:: console

  differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_REM_PRO_HiC/   -r <pathToInteractions> -g ensemblID_GeneName.txt -c 0.001 -n 10 necessaryInputFiles/estimatedScalesPerMotif_1.9.txt -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${genome}

Example 4: Compute a proper random background control and highlight the cell type specific TFs
---------------------------------------------------------------------------------------------

To perform a random background sampling the optional parameters -j, -k, -l and -q need to be specified. We recommend to sample at least 100 background rounds, meaning set -j to 100. The random SNPs are sampled from the dbSNP database. We provide the corresponding file in the Zenodo repository (unzipped file: dbSNPs_sorted.txt) which is used to specify the flag -k. To allow reproducible results, we ask the user to set a random seed via the -l flag. Please use varying random seed for runs with different input SNPs. The flag -q is used to speed up the background sampling by exclude TFs, which did not have or did have less significant differential binding affinities on the input SNPs. Per default -q is set not 0, meaning only TFs with at least 1 significant change in the binding affinity are considered in the background sampling. 
Further we recommend running SNEEP in the parallel mode by specifying the number of threads via the -n flag. 

A possible SNEEP run with background sampling can look as following: 

.. code-block:: console

  differentialBindingAffinity_multipleSNPs -o examples/SNEEP_output_background_sampling/ -c 0.001 -s necessaryInputFiles/estimatedScalesPerMotif_1.9.txt -b necessaryInputFiles/frequency.txt -x necessaryInputFiles/transition_matrix.txt  -n 20 -j 100 -k <pathTodbSNP> -l 2 -q 0 -r <pathToInteractions> -g ensemblID_GeneName.txt  examples/combined_Jaspar2022_Hocomoco_Kellis_human_transfac.txt  examples/SNPs_EFO_0000612_myocardial_infarction.bed ${geno    me}
