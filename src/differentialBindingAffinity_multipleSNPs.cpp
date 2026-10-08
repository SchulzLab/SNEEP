#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <bitset>
#include <ctime> //for time
#include <chrono> //for time
#include <iomanip>  // for precision

//own classes
#include "pvalue_copy.hpp"
#include "Matrix_new.hpp"
#include "callBashCommand.hpp"
#include "HandleInOutput.hpp"
#include "sampleRandomRsIDs2.hpp"
#include "stringUtils.hpp"
#include "parallelError.hpp"

//neccessary for function strcmp
#include <stdio.h>
#include <string.h>

//important for reading a directory
#include <dirent.h>
#include <sys/types.h>
#include <getopt.h> //parse the command line arguments

//for log
#include <math.h>

//for parallelization 
#include <omp.h>

using namespace std;

//------------------------------
//Global variables default values 
const int COMPLEMENT[] = {0,4,0,3,1,0,0,2,0,0,0,0,0,0,5};// considers also N
const int POSITION[] = {0,1,0,2,4,0,0,3,0,0,0,0,0,0,5}; // considers also N 
const double EPSILON = 0.001; // accuracy of the PWMs, EPSILON which is added to zeor PWM entries

//functions
vector<string> readDirectory(const char *path); //stores all files of a directory
vector<double> readFrequence(string frequence); //read the base pair frequence of the genome or 0.25 per base
vector<Matrix<double>> PFMsToPWMs( vector<string>& PWM_files, vector<double>& freq, string path_to_pwms); //determines PWMs based on the input count matrices
int probProSeq(Matrix<double>& PWM, const string& line, const vector<double>& pvalues, double pvalue, int pos_snp, vector<double>& prob_sequences);//determines the binding affinity pro sequence
double probKmer(Matrix<double>& PWM,string::const_iterator start, string::const_iterator end); // determines binding affinity for each kmer
double probKmerComplement(Matrix<double>& PWM,string::const_iterator start, string::const_iterator end); // same for the complement
int getMutatedPos(string mut_pos); //calculates which position in the input sequence is the mutated position 
string getChr(string mut_pos); //return the chr of  chr:start-end 
string createMutatedSeq(string header, string& seq, vector<string>& splittedHeader); //determines the mutated sequence based on the given SNP
//double differentialBindingAffinity(double first_elem, double second_elem, unordered_map<double, double>& pre_log,  double& log_, double& scale, int numberKmers, int& round);
double differentialBindingAffinity(double first_elem, double second_elem,  double& log_, double& scale, int numberKmers);
string forwardOrReverse(int distance_); //determines on which strand the motif bindes
vector<string> parseHeader(string header, string delim); //parse the header of the fasta file
void writeHeadersOutputFiles(bool writeOutput,  string& currentOutput, string seq,string header, string mutSeq, vector<string>& splittedHeader, int& pos, string& chr);
void determineMAFsForSNPs(string SNPsFile,vector<double>& MAF);
void determineMAFsAndGCForSNPs(string fastaFile, vector<double>& MAF, vector<double>& GC, vector<string>& warnings);
bool sortbyth(const tuple<double, string>& a, const tuple<double, string>& b);
bool sortby(const tuple<double, string , string, string>& a, const tuple<double,string,  string, string>& b);
double cdf_laplace_abs_max(double scale, double numberKmers, double value);

// result of scoring one motif for one SNV: the k-mer with the maximal differential binding score (smallest p-value)
struct MotifHit{
	bool sigHit = false; // at least one k-mer (wildtype or mutated) with a binding p-value <= threshold (-p)
	double pvalue = 1.0; // p-value of the maximal differential binding score
	int maxPosSigHit = 0; // genomic position of the last base of the best k-mer
	string orientation = ""; // (f) or (r)
	double val1 = 0.0, val2 = 0.0, log = 0.0; // binding p-values wildtype and mutated sequence, log ratio
	int lenMotif = 0;
};
MotifHit scoreMotif(Matrix<double>& PWM, const vector<double>& pvalues, double scale, int lenMotif, const string& wildtypeSeq, const string& mutSeq, double pvalueThreshold, int pos, const string& motif, const string& chr, string* allOutput);
string formatHit(const MotifHit& hit, const string& chr, int pos);
string fastaName(const string& headerLine);
string scientific(double value);

int runSNEEP(int argc, char *argv[]);

// errors (e.g. a failed shell command) stop SNEEP with a message and exit code 1
// (exceptions inside OpenMP parallel regions cannot be caught here, the message is printed before)
int main(int argc, char *argv[]){
	try{
		return runSNEEP(argc, argv);
	} catch (exception& e){
		cerr << "ERROR: " << e.what() << endl;
		return 1;
	}
}

int runSNEEP(int argc, char *argv[]){

	cout << "HELLO github sneep version" << endl;
	
	//to output doubles whith a total of 17 digits
	typedef numeric_limits< double > dbl;
	std::cout.precision(dbl::max_digits10 - 1);

	//handle input 
	InOutput io; //create object
	try{
		io.parseInputPara(argc, argv); //set all paths according to input 
	} catch (exception& e){
		if (string(e.what()) == "help function end"){ // -h
			return 0;
		}
		cout << "ERROR: " << e.what() << endl;
		return 1;
	}

	//create instance bashCommand -> moved to handleInOutput.hpp -> parseInputPara
	BashCommand bc(io.getGenome()); //constructor bashcommand class
	//make outputDir; an existing dir is only cleared if it contains a former SNEEP output (info.txt),
	//otherwise SNEEP stops to avoid deleting other data (e.g. with -o . or a mistyped path)
	bc.mkdir(io.getOutputDir(), "-p", false);
	vector<string> existingFiles;
	for (auto& file : readDirectory(io.getOutputDir().c_str())){
		if (file[0] != '.'){ // hidden files (e.g. .DS_Store) are ignored and not deleted
			existingFiles.push_back(file);
		}
	}
	if (!existingFiles.empty()){
		if (find(existingFiles.begin(), existingFiles.end(), "info.txt") == existingFiles.end()){
			cout << "ERROR: output directory " << io.getOutputDir() << " is not empty and contains no former SNEEP output (info.txt), please choose another directory (-o)" << endl;
			return 1;
		}
		bc.anyCommand("rm -r -f " + io.getOutputDir() + "*"); // remove the former SNEEP output
	}
	ofstream info = io.openFile(io.getInfoFile(), false); //open info file
	info << io << endl; // write settings 
	info.close();

	// check if snps file is bed-like format or CVF
	io.fileFormatVCF(); // if file is CVF format parse to bed-like format and store new file in outputDir as formatedSNPs.txt

	io.checkIfSNPsAreUnique(); // removes SNPs from inputSNP list which are not unique  and stores them in info file
	//cout << "checked if SNPs are unique" << endl;
	//overlap with REMs
	string SNPFile = io.getSNPs();
	int entriesSNPFile = io.CountEntriesFirstLine(SNPFile, '\t'); //count entries in the SNPFile

	//open result file and write header
	ofstream resultFile = io.openFile(io.getResultFile(), false);	

	//check if there is a REM file and if so determine overlap with the given SNPs
	string REMs = io.getREMs();
	// header of result.txt and of the background results (sampling/randomResult_<round>.txt), same format
	string resultHeader = (REMs != "") ? "SNP_position\tvar1\tvar2\trsID\tMAF\tpeakPosition\tTF\tTF-binding_position\tstrand\teffectedPositionInMotif\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tREM_positions\tensemblIDs\tgeneNames\tREMIds\tcoefficients\tpvalues_REM\tnormModelScore\tmeanDNase1Signal\tstdDNase1Signal\tconsortium\n" : "SNP_position\tvar1\tvar2\trsID\tMAF\tpeakPosition\tTF\tTF-binding_position\tstrand\teffectedPositionInMotif\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\n";
	resultFile << resultHeader;
	//unordered_map<string, vector<string>> SNPsToOverlappingREMs;
	if (REMs != ""){
		string output = io.getOverlappingREMs();
		bc.intersect(REMs, SNPFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
	}
	//determine overlapping peaks
	string footprintFile = io.getFootprints();
	int numberSNPs = 0;
	if (footprintFile != ""){ //if footprint or region file is given
		string output = io.getOverlappingFootprints();
		//determine intersection
		bc.intersect(SNPFile,footprintFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
		//parse bedfile 
		io.parseSNPsBedfile(output, entriesSNPFile); //skip insertions that are longer than 1 and remember additionl number of entries per line (the first 5 entries are not counted, since they are necessary), result sored in OutputDir + snpRegions.bed, use getSNPBedFile() 

	}else{
		//parse bedfile 
		io.parseSNPsBedfile(SNPFile, entriesSNPFile);
	}
	//check again if file is unique (if a snp overlaps with more than one REM we might include them publicated here
	io.checkUniqAgain();

	//call getFasta
	bc.getFasta(io.getSNPBedFile(), io.getSNPfastaFile(), "-name");

	// determine number of SNPs from fastaa file
	numberSNPs = io.getNumberSNPs(io.getSNPfastaFile());
	cout << "number SNPs: " << numberSNPs <<endl;

	//check which TFs are active and write this TFs seperated in files  (done with a python script) also determine frequence matrix from the count matrices
	string activeTFs = io.getActiveTFs();
	bc.mkdir(io.getPFMsDir(), "-p", false); //create dir
	bc.rm(io.getPFMsDir()); // remove PWMs if there are any
	if (activeTFs != ""){
		bc.callPythonScriptCheckActiveMotifs(io.getActiveTFs(), io.getPFMs(), io.getPFMsDir(), io.getEnsembleIDGeneName(), io.getActivityThreshold(), io.getOutputDir()); 
	}else{
		bc.callPythonScriptSplitPFMs(io.getPFMs(), io.getPFMsDir(), io.getOutputDir()); 
		//for snp selex data (see also callBashCommands for more details) 
		//bc.callPythonScriptSplitPFMsSELEX(io.getPFMs(), io.getPFMsDir(), io.getOutputDir()); 
		
	}
	unordered_map<string, double> scales; 
	if (io.getScaleFile() != ""){
		io.readScaleValues(io.getScaleFile(),scales);
	}

//	for(auto& i : scales){
//		cout << i.first << " " << i.second << endl;
//	}

	//convert PFMs internally to PWMs
	vector<string> PWM_files = readDirectory(io.getPFMsDir().c_str()); //stores names of all pwm files
	vector<double> freq = readFrequence(io.getFrequence()); //TODO determine from the data? is there a bedtools function?
	vector<Matrix<double>> PWMs = PFMsToPWMs(PWM_files,freq, io.getPFMsDir());

	//determine pvalues for PFMs 
	string motif = "";
	Matrix<double> transition_matrix(4,4,0.25); // default to definde transition matrix with 0.25 
	//cout << transition_matrix << endl;
	// if transition matrix is specified use this file otherwise assume all transitions are equally likely
	if (io.getTransitionMatrix() != ""){
		ifstream transition_file(io.getTransitionMatrix());
       		transition_file>> transition_matrix;
	}
	unordered_map<string, vector<double>> all_pvalues;
	ParallelError pvalueError; // see parallelError.hpp
	#pragma omp parallel for private(motif) num_threads(io.getNumberThreads())
	for(int j = 0; j < PWM_files.size(); ++j){ // iterate over all given motifs
		if (pvalueError.failed()) continue; // an error occurred in another iteration
		try{
			motif = PWM_files[j].substr(0 , PWM_files[j].size() -4); // set motif
			//cout << "motif name: " << motif << endl;
			pvalue  pvalue_obj(PWMs[j].ncol(), EPSILON); // EPSILON entspricht accuracy !!! 
			vector<double> pvalues = pvalue_obj.calculatePvalues(PWMs[j], freq, transition_matrix);
			#pragma omp critical (storePvalues)
			{  
				all_pvalues[motif] = pvalues;
			}
		}catch (const exception& e){ // exceptions must not leave the parallel region
			pvalueError.set(e.what());
		}
	}
	pvalueError.rethrow();
	
	//stores all sequences and the according header
	string line = "";
	string helper1 (1000, 'a'); //initialize as string of length 1000 with only a's
	string helper2 (2 * SNV_FLANK + 1, 'a'); //sequences: SNV_FLANK bp + SNV + SNV_FLANK bp
	vector<string> headers (numberSNPs, helper1); //to avoid reallocation
	vector<string> sequences (numberSNPs, helper2);
	//unordered_map<double, double> preLog; //stores log of pvalues, avoid to recalculate them

	//vector<unordered_map<double, double>> vec_overall_preLog;


	ifstream pro_seq(io.getSNPfastaFile()); //open fasta file
	for(int n = 0; n < numberSNPs; ++n){
		getline(pro_seq, line, '\n'); //header
		headers[n] = fastaName(line);
		getline(pro_seq, line, '\n'); //header
		sequences[n] = line;
	}
	pro_seq.close(); // close fasta file

	// check the alleles against the reference genome; mismatches are only reported, since SNEEP cannot tell their cause apart:
	// wrong coordinates, another genome build or alleles on the minus strand (for a non-palindromic SNV the two alleles and
	// their complements cover all four bases, so a minus strand cannot be distinguished from a wrong position)
	int numMismatches = 0;
	ofstream infoAlleles = io.openFile(io.getInfoFile(), true);
	for(int n = 0; n < numberSNPs; ++n){
		vector<string> fields = parseHeader(headers[n], ";"); //chr:start-end, var1, var2, rsID, ...
		char ref = ((int)sequences[n].size() > SNV_FLANK) ? toupper(sequences[n][SNV_FLANK]) : 'N';
		char var1 = toupper(fields[1][0]), var2 = toupper(fields[2][0]);
		if (ref != var1 and ref != var2){
			numMismatches++;
			infoAlleles << "WARNING allele: SNV " << fields[0] << ";" << fields[1] << ";" << fields[2] << ";" << fields[3] << ": reference base " << ref << " matches neither var1 nor var2\n";
		}
	}
	infoAlleles << "!\tnumSNPsAllelesNotMatchingReference: " << numMismatches << '\n';
	if (numMismatches > 0){
		ostringstream percent;
		percent << fixed << setprecision(1) << 100.0 * numMismatches / numberSNPs;
		string warning = "for " + to_string(numMismatches) + " of " + to_string(numberSNPs) + " SNVs (" + percent.str() + "%) neither var1 nor var2 matches the reference base (listed in info.txt)";
		if (100.0 * numMismatches / numberSNPs > 5.0){
			warning += "; this is a high fraction, please check that the coordinates are 0-based, the SNVs and the genome file have the same genome build, and the alleles are given for the plus strand";
		}
		cout << "WARNING: " << warning << endl;
		infoAlleles << "WARNING allele: " << warning << '\n';
	}
	infoAlleles.close();

	bool writeOutput = false;
	bool outputMax =  io.getMaxOutput();
	ofstream output;
	if (io.getOutputAll() != ""){
		writeOutput = true;
		output = io.openFile(io.getOutputAll(), false);
	}

	int numMotifs = PWM_files.size();
	cout << "numMotifs: " << numMotifs << endl;
	vector<string> motifNames; 
	vector<int> lenMotifs;
	for(int i = 0; i < numMotifs; ++i){ //store motif names only once
		motifNames.push_back(PWM_files[i].substr(0 , PWM_files[i].size() -4)); //set motif
		lenMotifs.push_back(PWMs[i].ncol());
	}	
	if (numMotifs == 0){
		cout << "ERROR: no motifs found in " << io.getPFMsDir() << " (splitting the motif file " << io.getPFMs() << " failed?)" << endl;
		return 1;
	}
	vector<double> motifScales; // scale per motif (0 if the motif is not in the scale file, then the p-value is 1)
	ofstream infoScales = io.openFile(io.getInfoFile(), true);
	for (auto& elem : motifNames){
		if (scales.count(elem) == 0 or scales[elem] <= 0.0){ // p-value of the differential binding score is always 1
			string reason = (scales.count(elem) == 0) ? "is not in the scale file" : "has no positive scale in the scale file";
			cout << "WARNING: motif " << elem << " " << reason << " (" << io.getScaleFile() << "), it is never significant (p-value 1)" << endl;
			infoScales << "WARNING scale: motif " << elem << " " << reason << " (" << io.getScaleFile() << "), it is never significant (p-value 1)\n";
		}
		motifScales.push_back(scales[elem]);
	}
	infoScales.close();
	//determine average counts for TFs, genes and REMs
	unordered_map<string, double> realData_TFs;
	for (auto& elem : motifNames){ //initialize the TF map
		realData_TFs[elem] = 0;
	}

	vector<tuple<double,string, string, string>>  resultAllSNPs; // pvalue , motifName, outputPart1 outputPart2
	int rounds = io.getRounds();

	ParallelError snvError; // see parallelError.hpp
	#pragma omp parallel for private(motif) num_threads(io.getNumberThreads())
	for (int i = 0; i <numberSNPs; ++i){ //iterates over all sequences which contain the SNP(each sequence: 50bp + SNP + 50bp)
		if (snvError.failed()) continue; // an error occurred in another iteration
		try{
			//define variables
			vector<tuple<double, string>> overallResultMax; //stores max result per motif
			unordered_map<string, string> helperOverallResult;

			//create mutated sequence
			vector<string> splittedHeader; //snp_pos var1 var2 peak_pos REM_posensemblId genName REMId coefficient pvalue consortium version reverse
			string wildtypeSeq = sequences[i];
			string mutSeq = createMutatedSeq(headers[i], wildtypeSeq, splittedHeader); // sets var1 in the wildtype sequence and returns the sequence with var2
			int pos = getMutatedPos(splittedHeader[0]); // position of the snp within the genome
			string chr = getChr(splittedHeader[0]);

			// store parts of the output that is same per snp in the following to avoid multiple access to splittedHeader
			string part1 = splittedHeader[0] + '\t' + splittedHeader[1] + '\t' + splittedHeader[2] + '\t' + splittedHeader[3] + '\t' + splittedHeader[4] + '\t' + splittedHeader[5];
			string part2 = "";
			if (REMs == ""){
				part2 = "\n";
			}else{
				part2 = '\t' + splittedHeader[6] + '\t' +  splittedHeader[7] + '\t'  + splittedHeader[8] + '\t' + splittedHeader[9] + '\t' + splittedHeader[10] + '\t' +  splittedHeader[11] + '\t' +  splittedHeader[12] + '\t' + splittedHeader[13] + '\t' + splittedHeader[14] +  '\t'  +  splittedHeader[15] +  '\n';
			}
		
			//write parts of the output
			string currentOutput = "", currentMaxOutput = "", currentResult = "";
			if (!mutSeq.empty()){//if mutSeq is not a fit, the stringis empty
				writeHeadersOutputFiles(writeOutput, currentOutput, headers[i], wildtypeSeq, mutSeq, splittedHeader, pos, chr);
			}
			for(int j = 0; j < numMotifs; ++j){ // iterate over all given motifs
				motif = motifNames[j]; //set motif
				MotifHit hit = scoreMotif(PWMs[j], all_pvalues.at(motif), motifScales[j], lenMotifs[j], wildtypeSeq, mutSeq, io.getPvalue(), pos, motif, chr, writeOutput ? &currentOutput : nullptr);
				if (hit.sigHit){ //check if there is a sig hit in one of the kmers otherwise skip diffBind score and pvalue correction
					if (outputMax == true){
						cout << hit.log << "\t" << lenMotifs[j] << "\t" << motif << endl;
					}
					// determined maxLog pro motifs (pro SNP)
					overallResultMax.push_back(make_tuple(hit.pvalue, motif)); //add  pvalue to overallResultMAx
					helperOverallResult[motif] = formatHit(hit, chr, pos) + '\t';
				}else{
					if (io.getPvalueDiff() == 1){ //for ASB/non-ASB testing (we need all snps in the output file) 
						#pragma omp critical 
						resultFile << part1 + '\t' + motif + "\t-\t-\t-\t0.0\t0.0\t0.0\t1.0" + (REMs != "" ? "\t.\t.\t.\t.\t.\t.\t.\t.\t.\t." : "") + "\n"; 
					}
				}
			}
	
			//sort overall_result per snp and all TFs
			string m = "";
			double p = 0.0;
			sort(overallResultMax.begin(), overallResultMax.end(), sortbyth);
			//cout << overallResultMax.size() << endl;
	
			for (int k = 0; k < overallResultMax.size(); ++k){ // k is the number of motifs -1 (starts from zero)
				p = get<0>(overallResultMax[k]); //pvalue as double
				m = get<1>(overallResultMax[k]); //motif name
				
				#pragma omp critical 
				//resultFile << part1 << '\t' << motif << '\t' << value   << std::scientific << maxDiffBinding << "\tneedToBeComputeds\t" << part2;
				//resultAllSNPs.push_back(make_tuple( p,m,  part1 + '\t' + m + '\t' +  helperOverallResult[m],  part2 ));
				if (p <= io.getPvalueDiff()){ // and p <= io.getPvalueDiff()){ // cutoff based on not fdr corrected pvalue for Jayas data (also for background sampling)
					resultFile << part1 << '\t' << m << '\t' <<  helperOverallResult[m] <<  std::scientific << p << part2;
					realData_TFs[m]+=1;// count number of TF hits seen in original data
				}
			}
		
			if (writeOutput){ // and (firstSeq[l] <= io.getPvalue() or secondSeq[l] <= io.getPvalue())){
				#pragma omp critical 
				output << currentOutput;
			}
		}catch (const exception& e){ // exceptions must not leave the parallel region
			snvError.set(e.what());
		}
	}
	snvError.rethrow();


	//sort resultAllSNPs based on p-value 
	resultFile.close();
	///outputMax.close();
	output.close();

	//open info file
	info = io.openFile(io.getInfoFile(), true); //open info file

	//write number motifs to output file 
	info << "!\tnumMotifs: " << numMotifs << endl; // number of used motisf

	//write TFs and co to files
	string outputDir = io.getOutputDir();
	ofstream outputTFs;
	outputTFs.open(outputDir + "/TF_count.txt");
	outputTFs << "."; //write header (all TF names)
	for (auto& elem : motifNames){
		outputTFs << '\t' << elem;
	}
	outputTFs << "\nrealData";

	for (auto& elem: motifNames){
		outputTFs << '\t' << realData_TFs[elem];
	}
	outputTFs << '\n';
	//Randomly sample SNPs
	if (rounds > 0){
	
		cout << "start random sampling" << endl;

		//update used motifs	
		vector<string> randomSampling_motifNames; 
		vector<int> randomSampling_lenMotifs;
		vector<Matrix<double>> randomSampling_PWMs;
		vector<double> randomSampling_scales;
		string c_motif = "";
		int randomSampling_numMotifs = 0;
		// if we want to compute a gene bacground we need to consider all TFs and cannot reduce the set of motifs (no speed up)
		for(int i = 0; i < numMotifs; ++i){ //store motif names only once
			
			c_motif = motifNames[i];
			if  (io.getGeneBackground() == true|| realData_TFs[c_motif] > io.getMinTFCount()){
				randomSampling_numMotifs++;
				randomSampling_motifNames.push_back(c_motif); //set motif
				randomSampling_lenMotifs.push_back(lenMotifs[i]);
				randomSampling_PWMs.push_back(PWMs[i]);
				randomSampling_scales.push_back(motifScales[i]);
			}
		}	
		cout << "number TFs: " << numMotifs << endl; 
		cout << "number TFs considered background sampling: " << randomSampling_numMotifs << endl; 

		// all files of the background rounds are written to <outputDir>/sampling (also with -i, where only the
		// random SNPs are read from the given directory)
		string randomDir = outputDir + "sampling";
		bc.mkdir(randomDir, "-p" , false);
		vector<string> SNP_filenames (rounds, ""); // initalize vector holding the random files 

		//initialize variables
		vector<double> MAF;
		vector<double> GC;
		double pvalue = io.getPvalue();

		if (io.getRandomSNPs() == ""){	
			cout << "sample snps case" << endl;
			vector<string> warnings;
			if (io.getGCMatching()){
				determineMAFsAndGCForSNPs(io.getSNPfastaFile(), MAF, GC, warnings); //read MAF and GC content (from the sequence) of the input SNPs
			}else{
				determineMAFsForSNPs(io.getSNPBedFile(), MAF); //read MAF distribution from input SNP file
			}

			cout << "before sampling" << endl;
			if (io.getGCMatching()){
				rsIDsampler s(0.01, 0.05, io.getdbSNPs(), MAF, GC); //initialize snp sampler, MAF x GC bins
				SNP_filenames = s.determineRandomSNPs(randomDir, rounds, io.getSeed()); //ddetermine random SNPs for number of rounds based on dbSNP file
				vector<string> samplingWarnings = s.getWarnings();
				warnings.insert(warnings.end(), samplingWarnings.begin(), samplingWarnings.end());
			}else{
				rsIDsampler s(0.01, io.getdbSNPs(), MAF); //initialize snp sampler, MAF bins
				SNP_filenames = s.determineRandomSNPs(randomDir, rounds, io.getSeed()); //ddetermine random SNPs for number of rounds based on dbSNP file
				vector<string> samplingWarnings = s.getWarnings();
				warnings.insert(warnings.end(), samplingWarnings.begin(), samplingWarnings.end());
			}
			//print warnings and store them in the info file
			ofstream info_ = io.openFile(io.getInfoFile(), true);
			for (auto& w : warnings){
				cout << "WARNING: " << w << endl;
				info_ << "WARNING background sampling: " << w << '\n';
			}
			info_.close();
			cout << "after sampling" << endl;
		}else{
			cout << "randomly sampled SNPs are provided" << endl; 
			for (int r = 0; r < rounds; r++){ // only read from the given directory
				SNP_filenames[r] = io.getRandomSNPs() + "/randomSNPs_" + to_string(r) + ".txt";
				if (!ifstream(SNP_filenames[r])){
					cout << "ERROR: random SNPs of round " << r << " not found: " << SNP_filenames[r] << " (-i must contain randomSNPs_0.txt ... randomSNPs_" << rounds - 1 << ".txt for -j " << rounds << ")" << endl;
					return 1;
				}
			}
		}

		string SNP_file = "", SNPs_overlappingREMs = "", bedFile = "", fastaFile = "", currentRound = "";
		ofstream randomResult;
		ifstream fasta;
		ParallelError roundError; // see parallelError.hpp
		#pragma omp parallel for private (currentRound, SNP_file, SNPs_overlappingREMs, bedFile, fastaFile, randomResult) num_threads(io.getNumberThreads())
		for (int r = 0; r < rounds; r++){
			if (roundError.failed()) continue; // an error occurred in another iteration
			try{

				currentRound = to_string(r); // store i as string
				SNP_file = SNP_filenames[r];
				SNPs_overlappingREMs = randomDir + "/randomSNPsOverlapingREMs_" + currentRound + ".bed";
				bedFile = randomDir + "/snpsRegions_" + currentRound + ".bed"; 
				fastaFile = randomDir + "/snpsRegions_" + currentRound + ".fa"; 
				randomResult.open(randomDir + "/randomResult_" + currentRound + ".txt"); //open result file
				if (REMs != ""){
			//		cout << "REMs" << endl;
					//write header output file
					randomResult << resultHeader;
					//intersect SNPs with REMs
					bc.intersect(REMs, SNP_file, SNPs_overlappingREMs, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
					//parse bedfile 
					io.parseRandomSNPs(SNP_file, SNPs_overlappingREMs, bedFile, r);
				}else{
					//write header output file
					randomResult << resultHeader;
					//parse bedfile 
					io.parseRandomSNPs(SNP_file, "", bedFile , r);
				}
				randomResult.close();
				//call getFasta
			//	cout << "getFATSA" << endl;
				bc.getFasta(bedFile, fastaFile, "-name");
				//read current fasta file
			}catch (const exception& e){ // exceptions must not leave the parallel region
				roundError.set(e.what());
			}
		}
		roundError.rethrow();
		for (int r = 0; r < rounds; r++){
		//for (int r = 628; r < rounds; r++){
			cout << "round: " << r << endl;

			unordered_map<string, double> TF_counts; //TF counts need to be count for every round seperatly
			for (auto& elem : motifNames){
				TF_counts[elem] = 0;
			}
			currentRound = to_string(r); // store i as string
			fastaFile = randomDir + "/snpsRegions_" + currentRound + ".fa"; 
			randomResult.open(randomDir + "/randomResult_" + currentRound + ".txt", std::ios_base::app); //open result file

			fasta.open(fastaFile);
			vector<string> current_headers (numberSNPs, helper1); //to avoid reallocation
			vector<string> current_sequences (numberSNPs, helper2);
			string current_line = "";
			for(int n = 0; n < numberSNPs; ++n){
				getline(fasta, current_line, '\n'); //header
				current_headers[n] = fastaName(current_line);
				getline(fasta, current_line, '\n'); //header
				current_sequences[n] = current_line;
			}
			fasta.close(); // close fasta file
	//		cout<< "after read fasta" << endl;
			
			vector<tuple<double,string, string, string>>  currentResultAllSNPs; 
			ParallelError backgroundError; // see parallelError.hpp
			#pragma omp parallel for  num_threads(io.getNumberThreads())
			for (int i = 0; i < numberSNPs; ++i){ //iterates over all sequences which contain the SNP(each sequence: 50bp + SNP + 50bp)
				if (backgroundError.failed()) continue; // an error occurred in another iteration
				try{
					//define variables
					vector<tuple<double, string>> current_overallResultMax; //stores max result per motif
					unordered_map<string, string> current_helperOverallResult;

					//create mutated sequence
					vector<string> current_splittedHeader; //snp_pos var1 var2 peak_pos REM_posensemblId genName REMId coefficient pvalue consortium version reverse
					string current_wildtypeSeq = current_sequences[i];
					string current_mutSeq = createMutatedSeq(current_headers[i], current_wildtypeSeq, current_splittedHeader); // sets var1 in the wildtype sequence and returns the sequence with var2
					int current_pos = getMutatedPos(current_splittedHeader[0]); // position of the snp within the genome
					string current_chr = getChr(current_splittedHeader[0]);
					for(int j = 0; j < randomSampling_numMotifs; ++j){ // iterate over all considered motifs
						string current_motif = randomSampling_motifNames[j]; //set motif
						MotifHit hit = scoreMotif(randomSampling_PWMs[j], all_pvalues.at(current_motif), randomSampling_scales[j], randomSampling_lenMotifs[j], current_wildtypeSeq, current_mutSeq, pvalue, current_pos, current_motif, current_chr, nullptr);
						if (hit.sigHit){
							//ouput maximal binding affinity for the current seq
							current_overallResultMax.push_back(make_tuple(hit.pvalue, current_motif));
							current_helperOverallResult[current_motif] = formatHit(hit, current_chr, current_pos) + '\t';
						}
					}
					//sort overall_result per seq and all TFs and store maximal one
					string m = "";
					double p = 0.0;

					//sort(current_overallResultMax.begin(), current_overallResultMax.end());
					sort(current_overallResultMax.begin(), current_overallResultMax.end(), sortbyth);
					for (int k = 0; k< current_overallResultMax.size(); ++k){
						p = get<0>(current_overallResultMax[k]); //pvalue
						m = get<1>(current_overallResultMax[k]); //motif name

						if (REMs == ""){
							//allows only one thread to write in the output 
							#pragma omp critical
							currentResultAllSNPs.push_back(make_tuple(p, m,current_splittedHeader[0] + '\t' + current_splittedHeader[1] + '\t' + current_splittedHeader[2] + '\t' + current_splittedHeader[3] + '\t' + current_splittedHeader[4] + '\t' + current_splittedHeader[5] + '\t' + m + '\t' + current_helperOverallResult[m], "\n"));
						}else{
							//allows only one thread to write in the output 
							#pragma omp critical
							currentResultAllSNPs.push_back(make_tuple(p, m, current_splittedHeader[0] + '\t' + current_splittedHeader[1] + '\t' + current_splittedHeader[2] + '\t' + current_splittedHeader[3] + '\t'  + current_splittedHeader[4] + '\t' + current_splittedHeader[5] + '\t' + m + '\t' +  current_helperOverallResult[m], '\t' + current_splittedHeader[6] + '\t' +  current_splittedHeader[7] + '\t'  + current_splittedHeader[8] + '\t' + current_splittedHeader[9] + '\t' + current_splittedHeader[10] + '\t' +  current_splittedHeader[11] + '\t' + current_splittedHeader[12] + '\t' + current_splittedHeader[13] + '\t' + current_splittedHeader[14] +  '\t'  + current_splittedHeader[15] + '\n'));
						}
					}
				}catch (const exception& e){ // exceptions must not leave the parallel region
					backgroundError.set(e.what());
				}
			}
			backgroundError.rethrow();

			//sort resultAllSNPs based on p-value 
			sort(currentResultAllSNPs.begin(), currentResultAllSNPs.end(),sortby);
			for (auto& elem : currentResultAllSNPs){ // get vector per entry
				double p = get<0>(elem);
				string m = get<1>(elem);
				if (p <= io.getPvalueDiff()){ // cutoff based on the not corrected p-value (as for the input SNPs)
					TF_counts[m] += 1; //count per round number of sig. TF hits
					randomResult << get<2>(elem) << std::scientific << p << get<3>(elem); // same format as result.txt
				}
			}
			randomResult.close();	

			//write output TF_counts for the current round
			outputTFs << r;
			for (auto& elem : motifNames){
				outputTFs << '\t' << TF_counts[elem];	
			}
			outputTFs << '\n';

		}
		outputTFs.close(); //close TF_counts file
	}
	return 0;
}

/*
*input: path to a directory
*output: vector that contains the names of all files of the directory
*WATCH OUT: it contains two files that are not neccessary for further calcultions: . and .. (they are excluded)
*/
vector<string> readDirectory(const char *path){

	vector <string> result;
	//stores file names
  	dirent* de;

  	DIR* dp;

  	dp = opendir(path);
	if (dp == NULL){
		throw invalid_argument("sorry can not open directory");
 	}
	de = readdir(dp); 

    	while (de != NULL ){
		if(strcmp(de->d_name, ".") and strcmp(de->d_name, "..")){
      			result.push_back(de->d_name);
		}
      		de = readdir(dp);
      	}	
    	closedir(dp);
	return result;
}


bool sortbyth(const tuple<double, string>& a, const tuple<double, string>& b){
    //return (get<2>(a) > get<2>(b));
	return (get<0>(a) > get<0>(b));
}

bool sortby(const tuple<double, string, string, string>& a, const tuple<double, string, string, string>& b){
    //return (get<2>(a) > get<2>(b));
	return (get<0>(a) > get<0>(b));
}
/*
/ read or determine MAF distribution of the input SNPs
*/
//TODO: what to do if MAF is not given? 
void determineMAFsForSNPs(string bedFile, vector<double>& MAF){
	string line = "";
	double maf = 0.0;
	ifstream bedfile(bedFile); //open sequence file
	while (getline(bedfile, line, '\n')){ // for each line extract overlapping SNPs
		getToken(line, ';'); //skip chr\tstart\tend\tchr:start-end
		getToken(line, ';'); //skip wildtype 
		getToken(line, ';'); //skip mutatant
		getToken(line, ';');//skip rsID
		try{ //throws an error when MAF is smaller than double precisoins allows -> set after		
			maf = stod(getToken(line, ';'));
		}catch (const std::out_of_range& oor){
			maf = 0.0;	
		}catch (const std::invalid_argument& ia){ // MAF not given as a number (e.g. "-") -> not given
			maf = -1.0;
		}
		MAF.push_back(maf);
	}
	bedfile.close();
	return;
}

/*
/ read MAF (from the header) and GC content (from the sequence) of the input SNPs from the fasta file (snpRegions.fa)
/ the SNP is located at position SNV_FLANK of the sequence (SNV_FLANK bp + SNP + SNV_FLANK bp), the GC content is determined in a window of +- GC_FLANK bp
/ SNPs without GC content (no A, C, G, T in the window or unexpected sequence length) get GC = -1 and are sampled by MAF only
*/
void determineMAFsAndGCForSNPs(string fastaFile, vector<double>& MAF, vector<double>& GC, vector<string>& warnings){
	string header = "", seq = "", snp = "";
	double maf = 0.0, gc = 0.0;
	int numN = 0, windowLength = 0;
	ifstream fasta(fastaFile); //open fasta file
	while (getline(fasta, header, '\n')){
		getline(fasta, seq, '\n');
		header = fastaName(header) + '\n'; // chr:start-end;var1;var2;rsID;MAF;...
		snp = getToken(header, ';'); //chr:start-end
		snp += ";" + getToken(header, ';'); //wildtype
		snp += ";" + getToken(header, ';'); //mutant
		snp += ";" + getToken(header, ';'); //rsID
		try{ //throws an error when MAF is smaller than double precisoins allows -> set to 0
			maf = stod(getToken(header, ';'));
		}catch (const std::out_of_range& oor){
			maf = 0.0;
		}catch (const std::invalid_argument& ia){ // MAF not given as a number (e.g. "-") -> not given
			maf = -1.0;
		}
		MAF.push_back(maf);

		if ((int)seq.size() != 2 * SNV_FLANK + 1){
			gc = -1;
			warnings.push_back("SNV " + snp + ": sequence length " + to_string(seq.size()) + " instead of " + to_string(2 * SNV_FLANK + 1) + ", no GC content, sampled by MAF only");
		}else{
			gc = gcContent(seq, SNV_FLANK, GC_FLANK, numN, windowLength);
			if (gc == -1){
				warnings.push_back("SNV " + snp + ": GC window contains only N, no GC content, sampled by MAF only");
			}else if (numN > windowLength / 2){
				warnings.push_back("SNV " + snp + ": " + to_string(numN) + " of " + to_string(windowLength) + " bases of the GC window are N, GC content might not be meaningful");
			}
		}
		GC.push_back(gc);
	}
	fasta.close();
	return;
}

/*
/ read or determine MAF distribution of the input SNPs
*/

vector<double> readFrequence(string frequence){

	vector<double> freq;
	if(frequence == ""){
		freq.push_back(0.25);
		freq.push_back(0.25);
		freq.push_back(0.25);
		freq.push_back(0.25);
		//throw invalid_argument("path to frequence.txt is not set");
	}else{
		ifstream file(frequence);
		if (not file.is_open())
			throw invalid_argument("Cannot open frequence.txt");
		string word = "";
		while (file >> word){
			freq.push_back(stod(word));
		}
	}
	return freq;	
}

vector<Matrix<double>> PFMsToPWMs( vector<string>& PWM_files, vector<double>& freq, string path_to_pwms){
	
	if(path_to_pwms == "")
		throw invalid_argument("no path to pwm files is set!");	
	Matrix<double> PWM_matrix;
	vector<Matrix<double>> PWMs(PWM_files.size(), PWM_matrix);

	double rounder = 1/EPSILON;
	double value = 0;

	double freq_max = *(max_element(freq.begin(), freq.end()));

	for (int k = 0; k < PWM_files.size(); ++k){// for all PWMs
		//cout << PWM_files[k] << endl;	
		ifstream PWM(path_to_pwms + PWM_files[k]);//open actual PWM file
		PWM >> PWM_matrix; //read PWM file as matrix
		// determine PWM, round matrix and shift it in such a way that only values > 0 are included in the matrix -> important for pvalue calculation
		for (int i = 1; i<= PWM_matrix.ncol(); i++){
			for (int j = 1; j <= 4; j++){
				value = log10(PWM_matrix(j,i)/freq[j-1]) - log10(EPSILON/freq_max) + EPSILON;
				//value = log(PWM_matrix(j,i)/freq[j-1]) - log(EPSILON/freq_max) + EPSILON;
				PWM_matrix(j,i) =  (round(value * rounder ) / rounder);
			}
		}
		//cout << PWM_matrix << endl; 
		PWMs[k] = PWM_matrix; // store matrix in map
	}
	return PWMs;
}

int probProSeq(Matrix<double>& PWM, const string& line, const vector<double>& pvalues, double pvalue, int pos_snp, vector<double>& pvalue_seq){

	//cout << line << endl;
	int min_ = pos_snp - SNV_FLANK;
	int max_ = pos_snp + SNV_FLANK;

	int length_motif = PWM.ncol();
	int start_seq = max(pos_snp - length_motif, min_-1);
	int end_seq = min(pos_snp, max_-1)-1;
	//cout <<"lengthMotif: " << length_motif << " start_seq: " << start_seq << " end seq: " << end_seq << endl;

        double prob = 0.0; // stores prob for each k_mer
        double prob_comp = 0.0; // stores prob for each k_mer of the complement
	int result = 0; //1 if a sig occurs, 0 if not

	string::const_iterator start = line.begin() - length_motif + 1;
	string::const_iterator end = line.begin();

	double rounder = 1/EPSILON;
	int size_pvalues = pvalues.size() -1;

	//store pos for vector that contains the pvalues
	int pos = 0; 
	int pos_comp = 0;
	int check = 0; // checks if there is a N in the sequence
	int counter = min_ - length_motif;
	
	double pvalue_prob = 0.0;
	double pvalue_prob_comp = 0.0;
        for(int i = 0; i < line.size() ; ++i){ // index starts at 0			
		if (line[i] == 'N' || line[i] == 'n'){
			check = length_motif; 
		}
		if (line[i] != 'N' && line[i] != 'n' && check == 0 && start >= line.begin() && counter >= start_seq && counter <= end_seq){
			prob = probKmer(PWM, start, end); //calculate probability of k_mer

			prob_comp = probKmerComplement(PWM,start, end); // determine prob for complement

			pos = round(size_pvalues -(rounder*prob)); //position of the p-value of prob in the vector
			pos_comp = round(size_pvalues -(rounder*prob_comp)); //position of the p-value of prob in the vector

			pvalue_prob = pvalues[pos];
			pvalue_seq.push_back(pvalue_prob);
			if (pvalue_prob <= pvalue){ 
				result = 1;
			}

			pvalue_prob_comp = pvalues[pos_comp];
			pvalue_seq.push_back(pvalue_prob_comp);
			if(pvalue_prob_comp <= pvalue){
				result = 1;
			}
			//cout << "prob: " << prob << " pos : " << pos << " pvalue: " << pvalue_prob << endl;
			//cout << "probComp: " << prob_comp << " posComp : " << pos_comp << " pvalueComp: " << pvalue_prob_comp << endl;
		}

		if (counter ==  end_seq)
			break;

		end++;// update last letter of k_mer
		start++;
		if (check > 0)
			check--;

		//counter for current position
		counter++;
	}
	return(result);
}

/*
*Input: PWM, iterators that point to the first letter of the k_mer and the last
*Output: probability of the actual k_mer
*/
double probKmer(Matrix<double>& PWM,string::const_iterator start, string::const_iterator end){
	//cout << "forward" << endl;
	
        double prob = 0.0; //sum of the PWM entries, starts at 0.0
	int counter = 1;
	for(auto i = start; i != end+1; i++){ //iterators over k_mer
//		cout << POSITION[((*i) & 0x0f)] << " " ; 
	//	cout << *i;
		prob+= PWM(POSITION[((*i) & 0x0f)], counter); //determine which letter we do consider and look up the accodirng prob in the PWM
		counter ++;
	}
	//cout << endl;
	//cout << prob << endl;
        return prob;
}
//same as for probKMer, only different iterates in reverse order and determne complement for actual letter
double probKmerComplement(Matrix<double>& PWM,string::const_iterator start, string::const_iterator end){

	//cout << "reverse" << endl;
        double prob = 0.0; //sum of the PWM entries, starts at 0.0
	int counter = 1;
	for(auto i = end; i != start-1; i--){// iterator in reverse order
//		cout << COMPLEMENT[((*i) & 0x0f)] << " ";
	//	cout << *i;
		prob+= PWM(COMPLEMENT[((*i) & 0x0f)], counter); //determine which letter we do consider and look up the accodirng prob in the PWM
		counter ++;
	}
	//cout << endl;
	//cout << prob << endl;
        return prob;
}

int getMutatedPos(string mut_pos){

	int del = mut_pos.find(":");
	int pos = stoi(mut_pos.substr(del+1));
//	cout << "pos: " << pos << endl;
	return pos;
}

string getChr(string mut_pos){

	int del = mut_pos.find(":");
	string chr = mut_pos.substr(0, del);
//	cout << "chr: " << chr << endl;
	return chr;
}

string createMutatedSeq(string header, string&  seq, vector<string>& splittedHeader){

	splittedHeader = parseHeader(header, ";"); //chr:start-end, var1, var2, rsID, MAF, ...
	char var1 = toupper(splittedHeader[1][0]);
	char var2 = toupper(splittedHeader[2][0]);
	// the SNV is at position SNV_FLANK of the plus-strand sequence; var1 (wildtype allele) is set in seq, var2 (alternative
	// allele) in the returned sequence; alleles not matching the reference are reported after reading the sequences
	string result = seq;
	seq[SNV_FLANK] = var1;
	result[SNV_FLANK] = var2;
	return result;
}

double cdf_laplace_abs_max(double scale, double numberKmers, double value){

	//p-value of the absolute maximum of numberKmers Laplace(0, scale) values: 1 - (1 - exp(-|value|/scale))^numberKmers,
	//computed as -expm1(numberKmers * log1p(-q)) to stay accurate for very small p-values (1 - (1 - q)^n is 0 below ~1e-16)
	double result = 0;
	if (scale > 0.0){ 
		double q = exp(-(abs(value))/ scale);
		result = -expm1(numberKmers * log1p(-q));
	}else{ // if scale is not defind for this motif length 
		result = 1.0;
	}
	return result;
}


//double differentialBindingAffinity(double first_elem, double second_elem, unordered_map<double, double>& pre_log, unordered_map<double, double>& helper_preLog, double& log_, double& scale, int numberKmers){
//double differentialBindingAffinity(double first_elem, double second_elem, unordered_map<double, double>& pre_log,  double& log_, double& scale, int numberKmers, int& round){
double differentialBindingAffinity(double first_elem, double second_elem,  double& log_, double& scale, int numberKmers){

	log_ = 0.0;
	double pvalue = 0.0;

	double log_first = 0.0, log_second = 0.0;
	
	log_first = log(first_elem);
	log_second = log(second_elem);

	log_ = log_first-log_second;

	//determine pvalue
	pvalue = cdf_laplace_abs_max(scale, numberKmers, log_);
	//cout << "pvalue: " << pvalue << endl;
	return pvalue;	
}

/*
* scores one motif for one SNV: binding p-values of all k-mers overlapping the SNV in the wildtype and the mutated sequence,
* differential binding score and its p-value for all k-mers with a binding p-value <= pvalueThreshold in one of the sequences,
* returns the k-mer with the maximal differential binding score (smallest p-value)
* allOutput: if given, all computed differential binding scores are appended (flag -a)
*/
MotifHit scoreMotif(Matrix<double>& PWM, const vector<double>& pvalues, double scale, int lenMotif, const string& wildtypeSeq, const string& mutSeq, double pvalueThreshold, int pos, const string& motif, const string& chr, string* allOutput){

	MotifHit hit;
	vector<double> firstSeq; //binding p-values wildtype seq (forward and reverse strand alternating)
	vector<double> secondSeq; //binding p-values mutated seq
	if (probProSeq(PWM, wildtypeSeq, pvalues, pvalueThreshold, SNV_FLANK, firstSeq) == 1){
		hit.sigHit = true;
	}
	if (probProSeq(PWM, mutSeq, pvalues, pvalueThreshold, SNV_FLANK, secondSeq) == 1){
		hit.sigHit = true;
	}
	if (!hit.sigHit){
		return hit;
	}
	int posSigHit = 0;
	double diffBinding = 0.0, log_ = 0.0;
	for (int l = 0; l < firstSeq.size(); ++l){ 
		if (firstSeq[l] <= pvalueThreshold or secondSeq[l] <= pvalueThreshold){
			posSigHit = pos + floor(l/2); //deterime position of the hit	
			diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l], log_, scale, lenMotif*2);
			if (diffBinding <= hit.pvalue){ // store the maximal differential binding score
				hit.pvalue = diffBinding;
				hit.maxPosSigHit = posSigHit;
				hit.orientation = forwardOrReverse(l);
				hit.val1 = firstSeq[l];
				hit.val2 = secondSeq[l];
				hit.log = log_;
				hit.lenMotif = lenMotif;
			}
			if (allOutput != nullptr){
				allOutput->append(motif + '\t' + chr + ":" + to_string(posSigHit) + forwardOrReverse(l) + '\t' + scientific(firstSeq[l]) +  "\t" +  scientific(secondSeq[l]) +  "\t" + to_string(log_) + '\t' +  scientific(diffBinding)  + '\n');
			}
		}
	}
	return hit;
}

/*
* output columns of a hit: TF-binding_position, strand, effectedPositionInMotif, pvalue_BindAff_var1, pvalue_BindAff_var2, log ratio
*/
string formatHit(const MotifHit& hit, const string& chr, int pos){
	int start = hit.maxPosSigHit - hit.lenMotif + 1;
	string posInMotif = (hit.orientation == "(f)") ? to_string(pos - start + 1) : to_string(hit.maxPosSigHit + 1 - pos);
	return chr + ':' + to_string(start) + "-" + to_string(hit.maxPosSigHit + 1) + '\t' + hit.orientation + '\t' + posInMotif + '\t' + scientific(hit.val1) + '\t' + scientific(hit.val2) + '\t' + to_string(hit.log);
}

/*
* name of a FASTA record from its header line (">name"); newer bedtools versions (getfasta -name) append the
* coordinates as "name::chr:start-end", which would otherwise end up in the last field of the name (e.g. peakPosition)
*/
string fastaName(const string& headerLine){
	string name = headerLine.substr(1);
	size_t pos = name.find("::");
	if (pos != string::npos){
		name = name.substr(0, pos);
	}
	return name;
}

/*
* value in scientific notation with 6 decimals (e.g. 2.960000e-04), as the p-values in result.txt
*/
string scientific(double value){
	ostringstream result;
	result << std::scientific << value;
	return result.str();
}

string forwardOrReverse(int distance_){

	if ((distance_ % 2) == 0)
		return("(f)");
	else
		return( "(r)");
}
//extract one element from the header
vector<string> parseHeader(string header, string delim){

	vector<string> result;	
	int pos = header.find(delim);
	string elem = "";
	while ((pos = header.find(delim)) != std::string::npos) {
		elem = header.substr(0, pos);
		result.push_back(elem);
		header.erase(0, pos + 1);
	}
	result.push_back(header);	
	return result;
}

void writeHeadersOutputFiles(bool writeOutput, string& currentOutput,  string header, string seq, string mutSeq, vector<string>& splittedHeader, int& pos, string& chr){

	if (writeOutput){
		currentOutput.append("\n#\tinfo: " + header + "\n#\twildtyp: " + seq + "\n" + "#\tmutated: " +  mutSeq + "\n\n.\tpos\tBindAff_Var1\tBindAff_Var2\tlog(div)\tpvalue_diffBindAff\n");
	}
}
