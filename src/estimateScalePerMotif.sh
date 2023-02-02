## command to call : 

#06.12.22
##bash estimateScalePerMotif.sh 200000 combined_Jaspar_Hocomoco_Kellis_human_transfac.txt scalesPerMotif/ motifNames_combined_Jaspar_Hocomoco_Kellis_human.txt 1.9

numberSNPs=$1
pathToMotifs=$2
outputDir=$3
listOfMotifs=$4
thresholdEntropy=$5
pathToDbSNP=$6

mkdir ${outputDir}
mkdir ${outputDir}/maximalScores_${thresholdEntropy}/
mkdir ${outputDir}/distributionPerMotif_${thresholdEntropy}/

##step1 sample randomly snps from dbSNP catalog 
shuf -n ${numberSNPs} ${pathToDbSNP}  > ${outputDir}/randomlySelecteddbSNPs_${numberSNPs}.txt 

cut -f2,3,4,5,6,7,8 ${outputDir}/randomlySelecteddbSNPs_${numberSNPs}.txt >${outputDir}/randomlySelecteddbSNPs_${numberSNPs}_.txt 
mv ${outputDir}/randomlySelecteddbSNPs_${numberSNPs}_.txt ${outputDir}/randomlySelecteddbSNPs_${numberSNPs}.txt 

## run SNEEP for the snps
## set -m flag to get D_max printed to the terminal
time ./src/differentialBindingAffinity_multipleSNPs -o ${outputDir}/sneep_${thresholdEntropy}/ -m -n 1 -p 1 -c 1.0  -j 0  ${pathToMotifs} ${outputDir}/randomlySelecteddbSNPs_${numberSNPs}.txt  /home/nbaumgarten/hg38.fa >${outputDir}/maxDiffBindingScores_randomlySelectedSNPs_${numberSNPs}_${thresholdEntropy}.txt

echo -e "motif\testimatedScale\tMSE\toptimizedScale\toptimizedMSE\tmotifLength\tnumberOfKmers" > ${outputDir}/estimatedScalesPerMotif_${thresholdEntropy}.txt
while read -r line; 
do 
	echo "$line"
#	##grep per motif the corrsponding maximal differential binding scores
	echo -e "diffBind\tlength\tmotif" > ${outputDir}/maximalScores_${thresholdEntropy}/${line}_maximalScores.txt 
	grep -P "\t${line}$" ${outputDir}/maxDiffBindingScores_randomlySelectedSNPs_${numberSNPs}_${thresholdEntropy}.txt >> ${outputDir}/maximalScores_${thresholdEntropy}/${line}_maximalScores.txt

	wc -l ${outputDir}/maximalScores_${thresholdEntropy}/${line}_maximalScores.txt 
	## determine scale per motif and write to file
	python src/findZerosEachMotif.py  ${outputDir}/maximalScores_${thresholdEntropy}/${line}_maximalScores.txt  ${outputDir}/estimatedScalesPerMotif_${thresholdEntropy}.txt ${outputDir}/distributionPerMotif_${thresholdEntropy}/ sneep_randomlySelectedDbSNPs_${thresholdEntropy}/motifInfo.txt

done < ${listOfMotifs}

