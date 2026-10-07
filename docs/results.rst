====================
SNEEP result files
====================

The resulting files are either provided in the user-specified directory or as defaults in SNEEP_output/.

Main result file (result.txt)
==============================

The file lists all the SNPs that cause a significant change in the binding score for a TF (sorted in reserve order, the smallest p-value is listed as the last entry). A single SNP can affect multiple TFs, which results in multiple lines in this file. The first 14 entries are the following: 

-	*SNP_position* chr:start-end (0-based),
-	*var1* wildtype allele (usually the reference allele),
-	*var2* alternative allele (e.g., the effect allele),
-	*rsID* 
-	*MAF* 
-	*peakPosition* if an additional footprint/open chromatin region file is used (flag -f), the position of the overlapping region is given; otherwise . (meaning the option was not used).
-	*TF* name of the affected TF 
-	*TF-binding_position* genomic region, which is bound by the affected TF 
-	*strand* strand to which the TF binds (r) -> reverse strand, (f) -> forward
-	*effectedPositionInMotif* position within the TF binding motif, which is affected by the SNP
-	*pvalue_BindAff_var1* p-value of the binding affinity for var1 
-	*pvalue_BindAff_var2* p-value of the binding affinity for var2
-	*log_pvalueBindAffVar1_pvalueBindAffVar2* logarithmic ratio of the p-value for the binding score of var1 to the p-value for the binding score of var2, i.e., the differential TF binding score D = log(p(var1)/p(var2)) of our publication. A positive value means that var2 binds the TF better than var1 (smaller p-value), i.e., the exchange from the wildtype allele (var1) to the alternative allele (var2) may lead to a gain of a TF binding site; a negative value indicates a possible loss
-	*pvalue_DiffBindAff* corresponding p-value for the previous entry (log_pvalueBindAffVar1_pvalueBindAffVar2)

If the current SNPs overlap with an epigenetic interaction, the following entries, provide more information about the interaction regions. Otherwise, the remaining entries are filled with ‘.’.

-	*REM_positions* genomic region of the epigentic interaction 
-	*ensemblIDs* ensemble ID of the associated gene
-	*geneNames* gene name of the associated gene
-	*REMIds* unique identifier of the epigenetic interactions (if it is an identifier for a REM, it can be used to search for the region in our EpiRegio database) 

The last 8 entries contain information specific to EpiRegio. For more information, please have a look at `epiregio.de <https://epiregio.de>`_.


Info file (info.txt) 
=====================

The file info.txt holds the input parameters used for the SNEEP run. It also contains the number of SNPs for which neither var1 nor var2 matches the reference genome (numSNPsAllelesNotMatchingReference) together with these SNPs (lines starting with "WARNING allele"), warnings about motifs that can never be significant because they are missing in the scale file or have no positive scale (lines starting with "WARNING scale") and warnings of the background sampling (lines starting with "WARNING background sampling"), e.g., if a MAF or GC bin could not be filled as intended (see `flag -s <https://sneep.readthedocs.io/en/latest/detailsInput.html#flag-s-match-the-gc-content-in-the-background-sampling>`_).

Results of the random background sampling
=========================================

The results for the random background sampling can be found in the directory sampling and include the randomly sampled SNPs (randomSNPs_<round>.txt) and the results per round (randomResult_<round>.txt, same columns as result.txt). The number of significant hits per TF in the input data (row realData) and in each round is summarized in TF_count.txt in the output folder.



