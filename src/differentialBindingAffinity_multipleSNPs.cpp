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

//neccessary for function strcmp
#include <stdio.h>
#include <string.h>

//important for reading a directory
#include <dirent.h>
#include <sys/types.h>
#include <getopt.h> //einlesen der argumente

//for log
#include <math.h>

//for parallelization 
#include <omp.h>

using namespace std;


//------------------------------
//Global variables default values 
int COMPLEMENT[] = {0,4,0,3,1,0,0,2,0,0,0,0,0,0,5};// considers also N
int POSITION[] = {0,1,0,2,4,0,0,3,0,0,0,0,0,0,5}; // considers also N 
double EPSILON = 0.001; // accuracy of the PWMs, EPSILON which is added to zeor PWM entries
vector<double> SCALES = {0,0,0,0,0,0, 0.5099464874518722, 0.4667031501892002 , 0.433783753309124 , 0.40715127687298763, 0.3826005840545816 , 0.3787599566935292, 0.34774610222605595 , 0.34514695643150634,  0.334825561698675, 0.31438987031152266  , 0.3011626326230248, 0.2917197382561486, 0.2889961985763163, 0.2834156088976374, 0.2657082044936282, 0.28252991969868096, 0.2327379244822218, 0.2565442603258892 , 0.24312782127390908 , 0,0,0.22659238366923132 , 0,0,0,0,0,0,0.2116550921968819, 0.21663464972735508}; // from motif length 6 to 35 

//functions
vector<string> readDirectory(const char *path); //stores all files of a directory
vector<double> readFrequence(string frequence); //read the base pair frequence of the genome or 0.25 per base
vector<Matrix<double>> PFMsToPWMs( vector<string>& PWM_files, vector<double>& freq, string path_to_pwms); //determines PWMs based on the input count matrices
int probProSeq(Matrix<double>& PWM, string line, vector<double>& pvalues, double pvalue, int pos_snp, vector<double>& prob_sequences);//determines the binding affinity pro sequence
double probKmer(Matrix<double>& PWM,string::iterator start, string::iterator end); // determines binding affinity for each kmer
double probKmerComplement(Matrix<double>& PWM,string::iterator start, string::iterator end); // same for the complement
int getMutatedPos(string mut_pos); //calculates which position in the input sequence is the mutated position 
string getChr(string mut_pos); //return the chr of  chr:start-end 
int getSeqStart(string header); //same but return the start
int getSeqEnd(string header); //return the end
string createMutatedSeq(string header, string& seq, vector<string>& splittedHeader); //determines the mutated sequence based on the given SNP
double cdf_laplace(double my, double beta, double value); //determine laplace distribution
double differentialBindingAffinityBackground(double first_elem, double second_elem, unordered_map<double, double>& pre_log, unordered_map<double, double>& helper_preLog, double& log_, double& scale, int numberKmers);
//double differentialBindingAffinity(double first_elem, double second_elem, unordered_map<double, double>& pre_log,  double& log_, double& scale, int numberKmers, int& round);
double differentialBindingAffinity(double first_elem, double second_elem,  double& log_, double& scale, int numberKmers);
string forwardOrReverse(int distance_); //determines on which strand the motif bindes
vector<string> parseHeader(string header, string delim); //parse the header of the fasta file
void writeHeadersOutputFiles(bool writeOutput,  string& currentOutput, string seq,string header, string mutSeq, vector<string>& splittedHeader, int& pos, string& chr);
void determineMAFsForSNPs(string SNPsFile,vector<double>& MAF);
string getToken(string& line, char delim);
bool sortbyth(const tuple<double, string>& a, const tuple<double, string>& b);
bool sortby(const tuple<double, string , string, string>& a, const tuple<double,string,  string, string>& b);
double cdf_laplace_abs_max(double scale, double numberKmers, double value);

int main(int argc, char *argv[]){

	
	//to output doubles whith a total of 17 digits
	typedef numeric_limits< double > dbl;
	std::cout.precision(dbl::max_digits10 - 1);


	//handle input 
	InOutput io; //create object
	try{
		io.parseInputPara(argc, argv); //set all paths according to input 
	} catch (exception& e){
		cout << e.what() << endl;
		return (0);
	}

	//create instance bashCommand
	BashCommand bc(io.getGenome()); //constructor bashcommand class
	bc.mkdir(io.getOutputDir(), "-p", true); //make outputDir
	ofstream info = io.openFile(io.getInfoFile(), false); //open info file
	info << io << endl; // write settings 
	info.close();

	io.checkIfSNPsAreUnique(); // removes SNPs from inputSNP list which are not unique  and stores them in info file
	//cout << "checked if SNPs are unique" << endl;
	//overlap with REMs
	string SNPFile = io.getSNPs();
	int entriesSNPFile = io.CountEntriesFirstLine(SNPFile, '\t'); //count entries in the SNPFile

	//open result file and write header
	ofstream resultFile = io.openFile(io.getResultFile(), false);	

	//check if there is a REM file and if so determine overlap with the given SNPs
	string REMs = io.getREMs();
	unordered_map<string, vector<string>> SNPsToOverlappingREMs;
	if (REMs != ""){
		//write header output file
		resultFile << "SNP_position\tvar1\tvar2\trsID\tMAF\tpeakPosition\tTF\tTF-binding_position\tstrand\teffectedPositionInMotif\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tfdr_corrected_pvalue\tREM_positions\tensemblIDs\tgeneNames\tREMIds\tcoefficients\tpvalues_REM\tnormModelScore\tmeanDNase1Signal\tstdDNase1Signal\tconsortium\n";
		string output = io.getOverlappingREMs();
		bc.intersect(REMs, SNPFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
	}else{
		
		//write header output file
		resultFile << "SNP_position\tvar1\tvar2\trsID\tMAF\tpeakPosition\tTF-binding_position\tstrand\teffectedPositionInMotif\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tfdr_corrected_pvalue\n";
	}
	//determine overlapping peaks
	string footprintFile = io.getFootprints();
	int numberSNPs = 0;
	if (footprintFile != ""){ //if footprint or region file is given
		string output = io.getOverlappingFootprints();
		//determine intersection
		bc.intersect(SNPFile,footprintFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
		//parse bedfile 
		numberSNPs = io.parseSNPsBedfile(output, entriesSNPFile); //skip insertions that are longer than 1 and remember additionl number of entries per line (the first 5 entries are not counted, since they are necessary), result sored in OutputDir + snpRegions.bed, use getSNPBedFile() 

	}else{
		//parse bedfile 
		numberSNPs = io.parseSNPsBedfile(SNPFile, entriesSNPFile);
	}
	cout << "number SNPs: " << numberSNPs <<endl;
	//call getFasta
	bc.getFasta(io.getSNPBedFile(), io.getSNPfastaFile(), "-name");
	//cout << "after getFastaa" << endl;

	//check which TFs are active and write this TFs seperated in files  (done with a python script) also determine frequence matrix from the count matrices
	string activeTFs = io.getActiveTFs();
	bc.mkdir(io.getPFMsDir(), "-p", false); //create dir
	bc.rm(io.getPFMsDir()); // remove PWMs if there are any
	if (activeTFs != ""){
		bc.callPythonScriptCheckActiveMotifs(io.getSourceDir(), io.getActiveTFs(), io.getPFMs(), io.getPFMsDir(), io.getEnsembleIDGeneName(), io.getActivityThreshold(), io.getOutputDir()); 
	}else{
		bc.callPythonScriptSplitPFMs(io.getSourceDir(), io.getPFMs(), io.getPFMsDir(), io.getOutputDir()); 
		//for snp selex data (see also callBashCommands for more details) 
		//bc.callPythonScriptSplitPFMsSELEX(io.getSourceDir(), io.getPFMs(), io.getPFMsDir(), io.getOutputDir()); 
		
	}
//	bool scalesPerMotifs = false;
	unordered_map<string, double> scales; 
	if (io.getScaleFile() != ""){
		io.readScaleValues(io.getScaleFile(),scales);
//		scalesPerMotifs = true;
	}

//	for(auto& i : scales){
//		cout << i.first << " " << i.second << endl;
//	}

	//convert PFMs internally to PWMs
	vector<string> PWM_files = readDirectory(io.getPFMsDir().c_str()); //stores names of all pwm files
	vector<double> freq = readFrequence(io.getFrequence()); //TODO datenspezifisch bestimmen gibts da ein ebdtool function?
	vector<Matrix<double>> PWMs = PFMsToPWMs(PWM_files,freq, io.getPFMsDir());

	//determine pvalues for PFMs 
	string motif = "";
	Matrix<double> transition_matrix(4,4);
	ifstream transition_file(io.getSourceDir() + "/necessaryInputFiles/transition_matrix.txt");
        transition_file>> transition_matrix;

	unordered_map<string, vector<double>> all_pvalues;
	#pragma omp parallel for private(motif) num_threads(io.getNumberThreads())
	for(int j = 0; j < PWM_files.size(); ++j){ // iterate over all given motifs
		motif = PWM_files[j].substr(0 , PWM_files[j].size() -4); // set motif
		//cout << "motif name: " << motif << endl;
		pvalue  pvalue_obj(PWMs[j].ncol(), EPSILON); // EPSILON entspricht accuracy !!! 
		vector<double> pvalues = pvalue_obj.calculatePvalues(PWMs[j], freq, transition_matrix);
		#pragma omp critical (storePvalues)
		{  
			all_pvalues[motif] = pvalues;
		}
	}
	

	//stores all sequences and the according header
	string line = "";
	string helper1 (1000, 'a'); //initialize as string of length 1000 with only a's
	string helper2 (105, 'a'); //sequences is of length 100
	vector<string> headers (numberSNPs, helper1); //to avoid reallocation
	vector<string> sequences (numberSNPs, helper2);
	//unordered_map<double, double> preLog; //stores log of pvalues, avoid to recalculate them

	//vector<unordered_map<double, double>> vec_overall_preLog;

	int sigEntries = 0, sigEntriesOverlappingREM = 0;

	ifstream pro_seq(io.getSNPfastaFile()); //open fasta file
	for(int n = 0; n < numberSNPs; ++n){
		getline(pro_seq, line, '\n'); //header
		headers[n] = line.substr(1);
		getline(pro_seq, line, '\n'); //header
		sequences[n] = line;
	}
	pro_seq.close(); // close fasta file

	bool writeOutput = false;
	bool outputMax =  io.getMaxOutput();
	ofstream output;
	if (io.getOutputAll() != ""){
		writeOutput = true;
		output = io.openFile(io.getOutputAll(), false);
	}
	//if (io.getMaxOutput() != ""){
	//	outputMax = io.openFile(io.getMaxOutput(), true);
	//	writeOutputMax = true;
	//}

	int numMotifs = PWM_files.size();
	cout << "numMotifs: " << numMotifs << endl;
	vector<string> motifNames; 
	vector<int> lenMotifs;
	for(int i = 0; i < numMotifs; ++i){ //store motif names only once
		motifNames.push_back(PWM_files[i].substr(0 , PWM_files[i].size() -4)); //set motif
		lenMotifs.push_back(PWMs[i].ncol());
	}	
	//determine average counts for TFs, genes and REMs
	unordered_map<string, double> realData_TFs;
	for (auto& elem : motifNames){ //initialize the TF map
		realData_TFs[elem] = 0;
	}

	vector<tuple<double,string, string, string>>  resultAllSNPs; // pvalue , motifName, outputPart1 outputPart2
	int rounds = io.getRounds();

	#pragma omp parallel for private(motif) num_threads(io.getNumberThreads())
	for (int i = 0; i <numberSNPs; ++i){ //iterates over all sequences which contain the SNP(each sequence: 50bp + SNP + 50bp)
		unordered_map<double, double> preLog; //stores log of pvalues, avoid to recalculate them
		//define variables
		vector<tuple<double, string>> overallResultMax; //stores max result
		//unordered_map<string, double> TFBindingValues;
		unordered_map<string, string> helperOverallResult;
		string chr = "",  mutSeq = "", direction = "", wildtypeSeq = "", genes = "", REMs_ = ""; //stores if the SNP annotation was for the forward or the reversed strand
		int posMut = 0, pos = 0; // position of the snp within the genome an relative within the given sequence

		//create mutated sequence
		vector<string> splittedHeader; //snp_pos var1 var2 peak_pos REM_posensemblId genName REMId coefficient pvalue consortium version reverse
		wildtypeSeq = sequences[i];
		mutSeq = createMutatedSeq(headers[i], wildtypeSeq, splittedHeader); // check if all SNPs are fitting to the reference this is done within create mutated seq (also check for reverse strand)
		pos = getMutatedPos(splittedHeader[0]);
		//cout << "pos: " << pos << endl;
		chr = getChr(splittedHeader[0]);
		//cout << "chr: " << chr << endl;
		posMut = 50; // since we considere 50 bp upstream and downstream the snp

		// store parts of the output that is same per snp in the following to avoid multiple access to splittedHeader
		string part1 = splittedHeader[0] + '\t' + splittedHeader[1] + '\t' + splittedHeader[2] + '\t' + splittedHeader[3] + '\t' + splittedHeader[4] + '\t' + splittedHeader[5];
		string part2 = "";
		if (REMs == ""){
			part2 = splittedHeader[6] + '\n';
		}else{
			part2 = splittedHeader[6] + '\t' +  splittedHeader[7] + '\t'  + splittedHeader[8] + '\t' + splittedHeader[9] + '\t' + splittedHeader[10] + '\t' +  splittedHeader[11] + '\t' +  splittedHeader[12] + '\t' + splittedHeader[13] + '\t' + splittedHeader[14] +  '\t'  +  splittedHeader[15] +  '\n';
		}
		
		//write parts of the output
		string currentOutput = "", currentMaxOutput = "", currentResult = "";
		if (!mutSeq.empty()){//if mutSeq is not a fit, the stringis empty
			writeHeadersOutputFiles(writeOutput, currentOutput, headers[i], wildtypeSeq, mutSeq, splittedHeader, pos, chr);
		}
		//#pragma omp parallel for private(motif) num_threads(io.getNumberThreads())
		for(int j = 0; j < numMotifs; ++j){ // iterate over all given motifs
			int sigHit = 0; //for MatrixSigHits
				
			//define vaiables
			//necessary for determinig differential binding affinity
			int posSigHit = 0, maxPosSigHit = 0, maxLenMotif = 0;
			string direction = "", maxOrientation = "";
			double diffBindingAff = 0.0; //stores result of the differential binding affinity
			double diffBinding = 0.0, log_ = 0.0, maxDiffBinding = 1.0, maxVal1 = 0.0, maxVal2 = 0.0, maxLog = 0.0, result = 0.0;//maxDiffBinding for determinig best hit per seq
			//variables store the information for the best hit er seq
			vector<double>::iterator it; //iterater to determine min pos in first_seq and second_seq
			vector<double> firstSeq; //stores binding affinity normal seq
			vector<double> secondSeq; // stores binding affinity mutated seq
		
			motif = motifNames[j]; //set motif
			//cout << "motif: " << motif << endl;
			int lenMotif = lenMotifs[j];
			//cout << lenMotif << endl;
			vector<double> pvalues = all_pvalues[motif]; //pvalues for current motif

			//calculate probabiblity for sequences of the current TF	
			if (probProSeq(PWMs[j], wildtypeSeq,  pvalues, io.getPvalue(), posMut, firstSeq) == 1){ //first_seq contains wildtype sequence
				sigHit  = 1;
			}
			//calculate pvalue of the mutated sequence	
			if (probProSeq(PWMs[j], mutSeq, pvalues, io.getPvalue(), posMut, secondSeq) == 1){ //second_seq contain mutated sequence
				sigHit = 1;
			}
			//int maxDiff_index = 0;
			if (sigHit == 1){ //check if there is a sig hit in one of the kmers otherwise skip diffBind score and pvalue correction
				for (int l = 0; l < firstSeq.size(); ++l){ 
			
					if (firstSeq[l] <= io.getPvalue() or secondSeq[l] <= io.getPvalue()){
						posSigHit = pos + floor(l/2); //deterime position of the hit	

						//determine differntial binding
						//if (scalesPerMotifs == true){
							//cout << motif << "\t" << lenMotif << "\t"<< scales[motif] << endl;
							//diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l], preLog, helper_preLog, log_, scales[motif], lenMotif*2);
							//diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l], preLog, log_, scales[motif], lenMotif*2, rounds);
						diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l], log_, scales[motif], lenMotif*2);
						//}else{
							//diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l], preLog, helper_preLog, log_, SCALES[lenMotif], lenMotif*2);
						//	diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l],  log_, SCALES[lenMotif], lenMotif*2);
						//}
						//cout <<   log_ << "\t" << lenMotif << "\t" << motif << endl;
						//stores info maximal diff binding affinity
						if (diffBinding <= maxDiffBinding){// and (firstSeq[l] <= io.getPvalue() or secondSeq[l] <= io.getPvalue())){
							maxDiffBinding = diffBinding;
							maxPosSigHit = posSigHit;
							maxOrientation = forwardOrReverse(l);
							maxVal1 = firstSeq[l];
							maxVal2 = secondSeq[l];
							maxLog = log_;
							maxLenMotif = lenMotif;
							//maxDiff_index = l;
						}
						//add current SNP info to output (currentOutput)
						if (writeOutput){ // and (firstSeq[l] <= io.getPvalue() or secondSeq[l] <= io.getPvalue())){
							currentOutput.append(motif + '\t' + chr + ":" + to_string(posSigHit) + forwardOrReverse(l) + '\t' + to_string(firstSeq[l]) +  "\t" +  to_string(secondSeq[l]) +  "\t" + to_string(log_) + '\t' +  to_string(diffBinding)  + '\n');
						}
					}
				}

				if (outputMax == true){
					cout << maxLog << "\t" << lenMotif << "\t" << motif << endl;
				}
				// determined maxLog pro motifs (pro SNP)
				overallResultMax.push_back(make_tuple(maxDiffBinding, motif)); //add  pvalue to overallResultMAx
				string value = "";
				if (maxOrientation  == "(f)"){
					value =  chr + ':' + to_string(maxPosSigHit- maxLenMotif+1) + "-" + to_string(maxPosSigHit+1) + '\t' +  maxOrientation + '\t' + to_string( pos - (maxPosSigHit- maxLenMotif+1) +1)  + '\t' + to_string(maxVal1) + '\t' + to_string(maxVal2) + '\t' + to_string(maxLog) + '\t';
				}else{
					value =  chr + ':' + to_string(maxPosSigHit- maxLenMotif+1) + "-" + to_string(maxPosSigHit+1) + '\t' +  maxOrientation + '\t' +  to_string(maxPosSigHit+1 - pos) + '\t'+  to_string(maxVal1) + '\t' + to_string(maxVal2) + '\t' + to_string(maxLog) + '\t';	
				}

				helperOverallResult[motif] = value;
			}else{
				//cout << "hier" << endl;
				if (io.getPvalueDiff() == 1){ //for ASB/non-ASB testing (we need all snps in the output file) 

					
					#pragma omp critical 
					resultFile << splittedHeader[0] + '\t' + splittedHeader[1] + '\t' + splittedHeader[2] + '\t' + splittedHeader[3] + '\t'  + splittedHeader[4] + '\t' + splittedHeader[5] + '\t' + motif + "\t-\t-\t-\t0.0\t0.0\t0.0\t1.0\t1.0\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n"; 

				}
			}
		}
	
		//sort overall_result per snp and all TFs
		string m = "";
		double p = 0.0;
		double c_p = 0.0; //stores fdr corrected pvalue
		sort(overallResultMax.begin(), overallResultMax.end(), sortbyth);
		//cout << overallResultMax.size() << endl;
	
		for (int k = 0; k < overallResultMax.size(); ++k){ // k is the number of motifs -1 (starts from zero)
			p = get<0>(overallResultMax[k]); //pvalue as double
			m = get<1>(overallResultMax[k]); //motif name
				
			#pragma omp critical 
			//resultFile << part1 << '\t' << motif << '\t' << value   << std::scientific << maxDiffBinding << "\tneedToBeComputeds\t" << part2;
			resultAllSNPs.push_back(make_tuple( p,m,  part1 + '\t' + m + '\t' +  helperOverallResult[m],  part2 ));
		}
		
		if (writeOutput){ // and (firstSeq[l] <= io.getPvalue() or secondSeq[l] <= io.getPvalue())){
			#pragma omp critical 
			output << currentOutput;
		}

		// add preLog to overall_preLog
		//#pragma omp critical 
		//vec_overall_preLog.push_back(preLog);

	}


	//sort resultAllSNPs based on p-value 
	sort(resultAllSNPs.begin(), resultAllSNPs.end(),sortby);
	//check sorting ;
	//for (auto& elem : resultAllSNPs){ // get vector per entry
	//	cout << to_string(get<0>(elem)) + '\t' + get<1>(elem)<< endl;
	//}

	//compute fdr	
	double snpsXmotifs = numMotifs * numberSNPs; 
	cout << "number of tests: "  << snpsXmotifs << endl;
	double rank = resultAllSNPs.size(); 

	cout << "rank bevor correction: " << rank << endl;
	double previous_pvalue = 1.0;

	//bool sigEntriesOnce = true;// necessary to count the sig entrie only once even if there are more than on TF affected
	//bool overlappingREMsOnce = true; //same here

	for (auto& elem : resultAllSNPs){ // get vector per entry
	
		//compute fdr
		double p = get<0>(elem);
		string m = get<1>(elem);
	//	if (sigEntriesOnce == true){
	//		sigEntries+=1;
	//		sigEntriesOnce = false;
	//	}
		double helper = min(p*(snpsXmotifs/rank), previous_pvalue);
	//	cout << get<1>(elem) << "\tp-value:\t" << p << "\tcorrected pvalue\t" << helper << "\trank\t" << rank << "\tsnpsXmotifs/rank\t"<< snpsXmotifs/rank << endl;
		previous_pvalue = helper;
		rank--;
		//if (helper <= io.getPvalueDiff()){ // and p <= io.getPvalueDiff()){ //cutoff based on fdr corrected pvalue
		// for jayas analysis
		if (p <= io.getPvalueDiff()){ // and p <= io.getPvalueDiff()){ // cutoff based on not fdr corrected pvalue for Jayas data (also for background sampling)
			#pragma omp critical
			realData_TFs[m]+=1;// count number of TF hits seen in original data

			resultFile << get<2>(elem) << std::scientific << p << '\t' << std::scientific << helper << '\t' << get<3>(elem); 
		}
	}
	cout << "rank after correction: " << rank << endl; 

	resultFile.close();
	///outputMax.close();
	output.close();

	//open info file
	info = io.openFile(io.getInfoFile(), true); //open info file

	//info << "!\tsigEntryOverlappingREMs: " << sigEntriesOverlappingREM << endl; //write number of entries which are significant and overlap with a REM
	//info << "!\tsigOnes: " << sigEntries << endl; // write number of entries which are significant

	//write number motifs to output file 
	info << "!\tnumMotifs: " << numMotifs << endl; // number of used motisf

	
	
	//write TFs and co to files
	string outputDir = io.getOutputDir();
	ofstream outputTFs;
	outputTFs.open(outputDir + "/TF_count.txt");
	double helper_g = 0;
	double helper_l = 0;
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
	
		//add values to preLog
		//cout << "add computed log values to overall_preLog" << endl;
		//unordered_map<double, double> overall_preLog; 
		//for(auto& elem : vec_overall_preLog){
		//	for(auto& e : elem){
		//		overall_preLog[e.first] = e.second;
		//	}
		//}
		//vec_overall_preLog.clear(); //empty helper_preLog

		cout << "start random sampling" << endl;


		//update used motifs	
		vector<string> randomSampling_motifNames; 
		vector<int> randomSampling_lenMotifs;
		vector<Matrix<double>> randomSampling_PWMs;
		string c_motif = "";
		int randomSampling_numMotifs = 0;
		for(int i = 0; i < numMotifs; ++i){ //store motif names only once
			
			c_motif = motifNames[i];
			if (realData_TFs[c_motif] > io.getMinTFCount()){
				randomSampling_numMotifs++;
				randomSampling_motifNames.push_back(c_motif); //set motif
				randomSampling_lenMotifs.push_back(lenMotifs[i]);
				randomSampling_PWMs.push_back(PWMs[i]);
			}
		}	
		cout << "number TFs: " << numMotifs << endl; 
		cout << "number TFs considered background sampling: " << randomSampling_numMotifs << endl; 

		//initialize variables
		vector<double> MAF;
		unordered_map<double, int> MAF_counter;
		double pvalue = io.getPvalue(), pvalueDiff = io.getPvalueDiff();
		string randomDir = outputDir + "sampling";
		bc.mkdir(randomDir, "-p" , false);
		determineMAFsForSNPs(io.getSNPBedFile(), MAF); //read MAF distribution from input SNP file
		sort(MAF.begin(), MAF.end(), std::less<double>()); //sort MAF with default operation <

		rsIDsampler s(0.01, io.getdbSNPs(), MAF, bc); //initialize snp sampler
		MAF_counter = s.splitMAFinBins(); // split original MAF distribution in bins
		//for (auto& elem : MAF_counter){ // get vector per entry
		//	cout << elem.first << " " << elem.second << endl;
		//}
		cout << "before sampling" << endl;
		vector<string> SNP_filenames = s.determineRandomSNPs(MAF_counter, randomDir, rounds, io.getNumberThreads(), io.getSourceDir(), io.getSeed()); //ddetermine random SNPs for number of rounds based on dbSNP file
		cout << "after sampling" << endl;

		string SNP_file = "", SNPs_overlappingREMs = "", bedFile = "", fastaFile = "", currentRound = "";
		ofstream randomResult;
		ifstream fasta;
		#pragma omp parallel for private (currentRound, SNP_file, SNPs_overlappingREMs, bedFile, fastaFile, randomResult) num_threads(io.getNumberThreads())
		for (int r = 0; r < rounds; r++){

			currentRound = to_string(r); // store i as string
			SNP_file = SNP_filenames[r];
			SNPs_overlappingREMs = randomDir + "/randomSNPsOverlapiingREMs_" + currentRound + ".bed";
			bedFile = randomDir + "/snpsRegions_" + currentRound + ".bed"; 
			fastaFile = randomDir + "/snpsRegions_" + currentRound + ".fa"; 
			randomResult.open(randomDir + "/randomResult_" + currentRound + ".txt"); //open result file
			if (REMs != ""){
		//		cout << "REMs" << endl;
				//write header output file
				randomResult << "SNP_position\tvar1\tvar2\trsID\tMAF\tpeakPosition\tTF\tTF-binding_position\tstrand\teffectedPositionInMotif\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tREM_positions\tensemblIDs\tgeneNames\tREMIds\tcoefficients\tpvalues_REM\tnormModelScore\tmeanDNase1Signal\tstdDNase1Signal\tconsortium\tversion\tSNP_strand\n";
				//intersect SNPs with REMs
				bc.intersect(REMs, SNP_file, SNPs_overlappingREMs, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
				//parse bedfile 
				io.parseRandomSNPs(SNP_file, SNPs_overlappingREMs, bedFile, r);
			}else{
				//write header output file
				randomResult << "SNP_position\tvar1\tvar2\tsrID\tMAF\tpeakPosition\tTF-binding_position\tstrand\teffectedPositionInMotif\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tSNP_strand\n";
				//parse bedfile 
				io.parseRandomSNPs(SNP_file, "", bedFile , r);
			}
			randomResult.close();
			//call getFasta
		//	cout << "getFATSA" << endl;
			bc.getFasta(bedFile, fastaFile, "-name");
			//read current fasta file
		}
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
				current_headers[n] = current_line.substr(1);
				getline(fasta, current_line, '\n'); //header
				current_sequences[n] = current_line;
			}
			fasta.close(); // close fasta file
	//		cout<< "after read fasta" << endl;
			
			vector<tuple<double,string, string, string>>  currentResultAllSNPs; 
			#pragma omp parallel for  num_threads(io.getNumberThreads())
			for (int i = 0; i < numberSNPs; ++i){ //iterates over all sequences which contain the SNP(each sequence: 50bp + SNP + 50bp)
				//unordered_map<double, double> preLog; //stores log of pvalues, avoid to recalculate them
				
				string currentResult = "";
				//define variables
				//vector<tuple<double, string, double>> current_overallResultMax; //stores max result
				vector<tuple<double, string>> current_overallResultMax; //stores max result
				unordered_map<string, string> current_helperOverallResult;
				unordered_map<string, double> current_TFBindingValues;
				string current_chr = "",  current_mutSeq = "", current_direction = "", current_genes = "", current_REMs = "", current_wildtypeSeq = ""; //stores if the SNP annotation was for the forward or the reversed strand
				int current_posMut = 0, current_pos = 0; // position of the snp within the genome an relative within the given sequence
	
				//create mutated sequence
				vector<string> current_splittedHeader; //snp_pos var1 var2 peak_pos REM_posensemblId genName REMId coefficient pvalue consortium version reverse
				current_wildtypeSeq = current_sequences[i];
				current_mutSeq = createMutatedSeq(current_headers[i], current_wildtypeSeq, current_splittedHeader); // check if all SNPs are fitting to the reference this is done within create mutated seq (also check for reverse strand)
				current_pos = getMutatedPos(current_splittedHeader[0]);
				current_chr = getChr(current_splittedHeader[0]);
				current_posMut = 50; // since we considere 50 bp upstream and downstream the snp
				for(int j = 0; j < randomSampling_numMotifs; ++j){ // iterate over all given motifs
	
				
					//define vaiables
					//necessary for determinig differential binding affinity
					int current_posSigHit = 0, current_maxPosSigHit = 0, current_maxLenMotif = 0, current_sigHit = 0; //for MatrixSigHits
					string current_direction = "", current_maxOrientation = "";
					double current_diffBindingAff = 0.0, current_diffBinding = 0.0, current_log_ = 0.0, current_maxDiffBinding = 1.0, current_maxVal1 = 0.0, current_maxVal2 = 0.0, current_maxLog = 0.0;//maxDiffBinding for determinig best hit per seq
					//variables store the information for the best hit er seq
					vector<double>::iterator current_it; //iterater to determine min pos in first_seq and second_seq
					vector<double> current_firstSeq; //stores binding affinity normal seq
					vector<double> current_secondSeq; // stores binding affinity mutated seq
		
					string current_motif = randomSampling_motifNames[j]; //set motif
					int lenMotif = randomSampling_lenMotifs[j];
					vector<double> pvalues = all_pvalues[current_motif]; //pvalues for current motif

					//avoid to calculate TF scores for TFs which are not considered in the original data -> also don't need them in the background
					//if (realData_TFs[current_motif] > 0){ // TODO: what threshold do we want to use here?
				
	
						//calculate probabiblity for sequences of the current TF
						if (probProSeq(randomSampling_PWMs[j], current_wildtypeSeq,  pvalues, pvalue, current_posMut, current_firstSeq) == 1){ //first_seq contains wildtype sequence
							current_sigHit  = 1;
						}
						//calculate pvalue of the mutated sequence	
						if (probProSeq(randomSampling_PWMs[j], current_mutSeq, pvalues, pvalue, current_posMut, current_secondSeq) == 1){ //second_seq contain mutated sequence
							current_sigHit = 1;
						}
						//determine differential binding affinity for all sig kmers in wildtype and mutated seq
						vector<pair<double, int>> current_allDiffBinding; 
						int current_maxDiff_index = 0;
						if (current_sigHit == 1){
							for (int l = 0; l < current_firstSeq.size(); ++l){ 
			
								if (current_firstSeq[l] <= pvalue or current_secondSeq[l] <= pvalue){
						
									current_posSigHit = current_pos + floor(l/2); //deterime position of the hit	


									// TODO: adapte for scale -> done but check if it works 
									//determine differntial binding
									//if (scalesPerMotifs == true){
										//cout << "motif: " << current_motif << " corresponding scale: " << scales[current_motif] << endl;
									current_diffBinding = differentialBindingAffinity(current_firstSeq[l], current_secondSeq[l], current_log_, scales[current_motif], lenMotif*2);
										//current_diffBinding = differentialBindingAffinityBackground(current_firstSeq[l], current_secondSeq[l],overall_preLog, preLog,  current_log_, scales[current_motif], lenMotif*2);
								//	}else{
								//		current_diffBinding = differentialBindingAffinity(current_firstSeq[l], current_secondSeq[l],  current_log_, SCALES[lenMotif], lenMotif*2);
										//current_diffBinding = differentialBindingAffinityBackground(current_firstSeq[l], current_secondSeq[l], overall_preLog, preLog,  current_log_, SCALES[lenMotif], lenMotif*2);
								//	}
									//current_allDiffBinding.push_back(make_pair(current_diffBinding, l));
								
									//stores info maximal diff binding affinity
									if (current_diffBinding <= current_maxDiffBinding){// and (current_firstSeq[l] <= pvalue or current_secondSeq[l] <= pvalue)){
										current_maxDiffBinding = current_diffBinding;
										current_maxPosSigHit = current_posSigHit;
										current_maxOrientation = forwardOrReverse(l);
										current_maxVal1 = current_firstSeq[l];
										current_maxVal2 = current_secondSeq[l];
										current_maxLog = current_log_;
										current_maxLenMotif = lenMotif;
										//current_maxDiff_index = l;
									}
								}	
							}
							//ouput maximal binding affinity for the current seq
							//current_overallResultMax.push_back(make_tuple(current_maxDiffBinding, current_motif, current_corrected_pvalue));
							current_overallResultMax.push_back(make_tuple(current_maxDiffBinding, current_motif));
							string value = "";
							if (current_maxOrientation  == "(f)"){
								value =  current_chr + ':' + to_string(current_maxPosSigHit- current_maxLenMotif+1) + "-" + to_string(current_maxPosSigHit+1) + '\t' +  current_maxOrientation + '\t' + to_string( current_pos - (current_maxPosSigHit- current_maxLenMotif+1) +1) + '\t' +   to_string(current_maxVal1) + '\t' + to_string(current_maxVal2) + '\t' + to_string(current_maxLog) + '\t' +  to_string(current_maxDiffBinding);
							}else{
								value =  current_chr + ':' + to_string(current_maxPosSigHit- current_maxLenMotif+1) + "-" + to_string(current_maxPosSigHit+1) + '\t' +  current_maxOrientation + '\t' + to_string(current_maxPosSigHit+ 1 - current_pos) + '\t' +   to_string(current_maxVal1) + '\t' + to_string(current_maxVal2) + '\t' + to_string(current_maxLog) + '\t' +  to_string(current_maxDiffBinding);
							}
							current_helperOverallResult[current_motif] = value;
							current_TFBindingValues[current_motif] = current_maxLog;
						}//else{
						//	current_overallResultMax.push_back(make_tuple(2.0, current_motif)); //we need all motifs to  fdr corrected the pvalue
						//}
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
							//currentResult.append(current_splittedHeader[0] + '\t' + current_splittedHeader[1] + '\t' + current_splittedHeader[2] + '\t' + current_splittedHeader[3] + '\t' + current_splittedHeader[4] + '\t' + current_splittedHeader[5] + '\t' + m + '\t' + to_string(helper) + '\t'  + current_helperOverallResult[m] + '\t' + current_splittedHeader[6] + '\n');
							//allows only one thread to write in the output 
							#pragma omp critical
							currentResultAllSNPs.push_back(make_tuple(p, m,current_splittedHeader[0] + '\t' + current_splittedHeader[1] + '\t' + current_splittedHeader[2] + '\t' + current_splittedHeader[3] + '\t' + current_splittedHeader[4] + '\t' + current_splittedHeader[5] + '\t' + m ,  current_helperOverallResult[m] + '\t' + current_splittedHeader[6] + '\n'));
						}else{
							//allows only one thread to write in the output 
							#pragma omp critical
							currentResultAllSNPs.push_back(make_tuple(p, m, current_splittedHeader[0] + '\t' + current_splittedHeader[1] + '\t' + current_splittedHeader[2] + '\t' + current_splittedHeader[3] + '\t'  + current_splittedHeader[4] + '\t' + current_splittedHeader[5] + '\t' + m + '\t' +  current_helperOverallResult[m], current_splittedHeader[6] + '\t' +  current_splittedHeader[7] + '\t'  + current_splittedHeader[8] + '\t' + current_splittedHeader[9] + '\t' + current_splittedHeader[10] + '\t' +  current_splittedHeader[11] + '\t' + current_splittedHeader[12] + '\t' + current_splittedHeader[13] + '\t' + current_splittedHeader[14] +  '\t'  + current_splittedHeader[15] + '\n'));
						}
						//}
					}
					//allows only one thread to write in the output 
					//#pragma omp critical (writeRandomResult)
					//randomResult << currentResult;
			//	}
			
			// add preLog to overall_preLog
			//#pragma omp critical 
			////vec_overall_preLog.push_back(preLog);
			}

			//sort resultAllSNPs based on p-value 
			sort(currentResultAllSNPs.begin(), currentResultAllSNPs.end(),sortby);
			//check sorting 
			//for (auto& elem : resultAllSNPs){ // get vector per entry
			//	cout << to_string(get<0>(elem)) + '\t' + get<1>(elem)<< endl;
			//}
			//compute fdr
			//cout << "number of tests: "  << snpsXmotifs << endl;
			double rank = currentResultAllSNPs.size(); 			
			double previous_pvalue = 1.0;

			for (auto& elem : currentResultAllSNPs){ // get vector per entry
	
				//compute fdr
				double p = get<0>(elem);
				string m = get<1>(elem);
				double helper = min(p*(snpsXmotifs/rank), previous_pvalue);
				previous_pvalue = helper;
				rank--;
				//if (helper <= io.getPvalueDiff()){// and p <= io.getPvalueDiff()){
				// for Jayas analysis
				if (p <= io.getPvalueDiff()){// and p <= io.getPvalueDiff()){
					#pragma omp critical
					TF_counts[m] += 1; //count per round number of sig. TF hits
					randomResult << get<2>(elem) << '\t' << to_string(helper) << '\t' << get<3>(elem); 
				}
			}
			randomResult.close();	

			//write output TF_counts for the current round
			outputTFs << r;
			for (auto& elem : motifNames){
				outputTFs << '\t' << TF_counts[elem];	
			}
			outputTFs << '\n';

			//add helper_preLog to original preLog
			//int grr = 0;
			//for(auto& elem : vec_overall_preLog){
			//	for(auto& e : elem){
			//		overall_preLog[e.first] = e.second;
			//	}
			//}
			//vec_overall_preLog.clear(); //empty helper_preLog
			//cout << "round r: " << r << " einträge helper_preLog: " << grr << endl;  
			//cout << "size preLog: " << preLog.size() << '\n'<<  endl;
			//helper_preLog.clear(); //empty helper_preLog
		}
		outputTFs.close(); //close TF_counts file
	}
	return 0;
}

/*
*input: path to a directory
*output: vector that contains the names of all files of the directory
*WATCH OUT: it contains two files that are not neccessary for further calcultions: . and .. (aktuell rausgenommen)
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

/*
/ splits string at given token (string befor token is returned, remaining line is shorter by the token)
*/
string getToken(string& line, char delim){
	int pos = 0;
	//cout << line << " " << line.find(delim) << endl;
	if(((pos = line.find(delim)) != std::string::npos) || ((pos = line.find('\n')) != std::string::npos)){
    		string token = line.substr(0, pos);
    		line.erase(0, pos + 1);
		return token;
	}else{
		throw invalid_argument ("invalid file format:" + line);
	}
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
/
/
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
		}
		MAF.push_back(maf);
	}
	return;
}

vector<double> readFrequence(string frequence){

	vector<double> freq;
	if(frequence == "")
		throw invalid_argument("path to frequence.txt is not set");
	
	ifstream file(frequence);
	if (not file.is_open())
		throw invalid_argument("Cannot open frequence.txt");
	string word = "";
	while (file >> word){
		freq.push_back(stod(word));
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

int probProSeq(Matrix<double>& PWM, string line, vector<double>& pvalues, double pvalue, int pos_snp, vector<double>& pvalue_seq){

	//cout << line << endl;
	int min_ = pos_snp -50;
	int max_ = pos_snp + 50;

	int length_motif = PWM.ncol();
	int start_seq = max(pos_snp - length_motif, min_-1);
	int end_seq = min(pos_snp, max_-1)-1;
	//cout <<"lengthMotif: " << length_motif << " start_seq: " << start_seq << " end seq: " << end_seq << endl;

        double prob = 0.0; // stores prob for each k_mer
        double prob_comp = 0.0; // stores prob for each k_mer of the complement
	int result = 0; //1 if a sig occurs, 0 if not

	string::iterator start = line.begin() - length_motif + 1;
	string::iterator end = line.begin();

	double rounder = 1/EPSILON;
	int size_pvalues = pvalues.size() -1;

	//store pos for vector that contains the pvalues
	int pos = 0; 
	int pos_comp = 0;
	int check = 0; // checks if there is a N in the sequence
	int counter = min_ - length_motif;
	
	double pvalue_prob = 0.0;
	double pvalue_prob_comp = 0.0;
        for(int i = 0; i < line.size() ; ++i){ // da index  bei 0 anfaengt oder?			
		if (line[i] == 'N' || line[i] == 'n'){
			check = length_motif; 
		}
		if (line[i] != 'N' && line[i] != 'n' && check == 0 && start >= line.begin() && counter >= start_seq && counter <= end_seq){
			prob = probKmer(PWM, start, end); //calculate probability of k_mer

			prob_comp = probKmerComplement(PWM,start, end); // determine prob for complement

			pos = round(size_pvalues -(rounder*prob)); //berechnet an welcher position im vector der pvalue für prob steht
			pos_comp = round(size_pvalues -(rounder*prob_comp)); //berechnet an welcher position im vector der pvalue für prob steht

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
double probKmer(Matrix<double>& PWM,string::iterator start, string::iterator end){
	//cout << "forward" << endl;
	
        double prob = 0.0; //muss 0.0 sein bei plus
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
double probKmerComplement(Matrix<double>& PWM,string::iterator start, string::iterator end){

	//cout << "reverse" << endl;
        double prob = 0.0; //muss 0.0 sein bei plus
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

int getSeqStart(string header){

	int del = header.find(":");
	int del2 = header.find("-");
	int start = stoi(header.substr(del+ 1, del2-del));
//	cout << "start: " << start << endl;
	return start;
}

int getSeqEnd(string header){

	int del = header.find(":");
	int del2 = header.find("-");
	int end = stoi(header.substr(del2+1));
//	cout << "end: " << end << endl;
	return end;
}

string createMutatedSeq(string header, string&  seq, vector<string>& splittedHeader){

	map<char , char> complement = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}}; //find complement of a base
	splittedHeader = parseHeader(header, ";"); //chr:start, var1, var2 
	char var1 = toupper(splittedHeader[1][0]); //cast to char
	char var2 = toupper(splittedHeader[2][0]); //cast to char
	int pos = getMutatedPos(splittedHeader[0]);
	string chr = getChr(splittedHeader[0]);	
	string direction;
	string result = seq;
	int start = pos - 50;
	int end = pos + 50;

	int mut = pos - start; //position of the mutation in relation to the sequence length (should be at position 50)
//	cout << mut << endl;
	//check for forward strand
	if (toupper(seq[mut]) == var1){ //forward strand var1 in reference seq
	//	cout << "var1 is reference genome forward" << endl;
		direction = "forward";
		result[mut] = var2;
	} 
	if (toupper(seq[mut]) == var2){ // forward strand var2 in reference seq
	//	cout << "var2 is reference genome forward" << endl;
		direction = "forward";
		//result[mut] = var1;
		seq[mut] = var1; // since result is a copy from the original seq ->  result[mut] = var2
	//check for reverse strand
	}/*else if (toupper(seq[mut]) == complement[var1]){ //reverse strand var1 in reference seq
		cout << "var1 is reference genome reverse" << endl;
		direction = "reverse";
		result[mut] = complement[var2];
	}else if (toupper(seq[mut]) == complement[var2]){ //reverse strand var2 in refernece seq
		cout << "var2 is reference genome reverse" << endl;
		direction = "reverse";
		result[mut] = complement[var1];
	}else{
		cout << "letter at position 50 (0-based): " << seq[mut] << " "  << chr << ":" << pos << " seq: " << seq << " var1: " << var1 << " var2: " << var2 << " header " << header << endl;
	//	throw invalid_argument("mutation passt nicht!");
		cout << "da stimmt was nicht !!!!" << endl;
	}*/
	splittedHeader.push_back(direction);
	return result;
}

double cdf_laplace(double my, double beta, double value){

	double result = 0.0;
	if (value < my){
		result = 0.5 * (exp(((value - my)/beta)));
	}else{ //value >= location
		result = 1 - (0.5 * (exp(-((value - my)/beta))));
	}
	//cout << "pvalue diffBind: " << result << endl;
	return result;
}

double cdf_laplace_abs_max(double scale, double numberKmers, double value){

	//cout << "scale: " << scale << endl;
	//cout << "numberKmers: " << numberKmers << endl;
	//cout << "value: " << value << endl;

	//compute cdf of laplace abs maximum 
	double result = 0;
	if (scale > 0.0){ 
		double cdf = pow((1 - exp(-(abs(value))/ scale )), numberKmers);
	//	cout << "CDF: " << cdf << endl;
		result = 1 - cdf;
	//	cout << "1 - cdf: " << result << endl;
	}else{ // if scale is not defind for this motif length 
		result = 1.0;
	}
	return result;
}



double differentialBindingAffinityBackground(double first_elem, double second_elem, unordered_map<double, double>& pre_log, unordered_map<double, double>& helper_preLog, double& log_, double& scale, int numberKmers){

	double div = 0.0;
	log_ = 0.0;
	double pvalue = 0.0;

	double log_first = 0.0, log_second = 0.0;

	auto found = pre_log.find(first_elem);

	//cout <<"first elem: " << first_elem << endl;
	//cout <<"second elem: " << second_elem << endl;
	if (found != pre_log.end()){ //check if the log of div was calculated bevore
		log_first = found->second;
	//	#pragma omp atomic
	//	counter++;
	}else{
		log_first = log(first_elem);
		helper_preLog[first_elem] = log_first;		
	}
	found = pre_log.find(second_elem);
	if (found != pre_log.end()){ //check if the log of div was calculated bevore
		log_second = found->second;
	//	#pragma omp atomic
	//	counter++;
	}else{
		log_second = log(second_elem);
		helper_preLog[second_elem] = log_second;		
	}
	//}
	log_ = log_first-log_second;
	//cout <<"log diffBinding: " << log_ << endl;

	//determine pvalue old versions
	//symmetric laplace(0,1) 
	//pvalue = 2* cdf_laplace(0.0, 1.0, -abs(log_));
	//pvalue = 2* cdf_laplace(0.0, 0.368, -abs(log_));
	//determine pvalue
	pvalue = cdf_laplace_abs_max(scale, numberKmers, log_);
	//cout << "pvalue: " << pvalue << endl;
	return pvalue;	
}


//double differentialBindingAffinity(double first_elem, double second_elem, unordered_map<double, double>& pre_log, unordered_map<double, double>& helper_preLog, double& log_, double& scale, int numberKmers){
//double differentialBindingAffinity(double first_elem, double second_elem, unordered_map<double, double>& pre_log,  double& log_, double& scale, int numberKmers, int& round){
double differentialBindingAffinity(double first_elem, double second_elem,  double& log_, double& scale, int numberKmers){

	double div = 0.0;
	log_ = 0.0;
	double pvalue = 0.0;

	double log_first = 0.0, log_second = 0.0;
	
	//if (round > 0){

	//	auto found = pre_log.find(first_elem);

		//cout <<"first elem: " << first_elem << endl;
		//cout <<"second elem: " << second_elem << endl;
	//	if (found != pre_log.end()){ //check if the log of div was calculated bevore
	//		log_first = found->second;
	//	}else{
	//		log_first = log(first_elem);
	//		pre_log[first_elem] = log_first;		
	//	}
	//	found = pre_log.find(second_elem);
	//	if (found != pre_log.end()){ //check if the log of div was calculated bevore
	//		log_second = found->second;
	//	}else{
	//		log_second = log(second_elem);
	//		pre_log[second_elem] = log_second;		
	//	}
	//}else{
	log_first = log(first_elem);
	log_second = log(second_elem);
	//}

	log_ = log_first-log_second;
	//cout <<"log diffBinding: " << log_ << endl;

	//determine pvalue old versions
	//symmetric laplace(0,1) 
	//pvalue = 2* cdf_laplace(0.0, 1.0, -abs(log_));
	//pvalue = 2* cdf_laplace(0.0, 0.368, -abs(log_));
	//determine pvalue
	pvalue = cdf_laplace_abs_max(scale, numberKmers, log_);
	//cout << "pvalue: " << pvalue << endl;
	return pvalue;	
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

	string var1 = splittedHeader[1]; //cast to char
	string var2 = splittedHeader[2]; //cast to char
	int start = pos - 50;
	int end = pos + 50;
	int mut = pos - start; //position of the mutation in relation to the sequence length (should be at position 51)

	if (writeOutput){
		currentOutput.append("\n#\tinfo: " + header + "\n#\twildtyp: " + seq + "\n" + "#\tmutated: " +  mutSeq + "\n\n.\tpos\tBindAff_Var1\tBindAff_Var2\tlog(div)\tpvalue_diffBindAff\n");
	}
	//header output_max
	//if (writeMaxOutput){
	//	currentMaxOutput.append("snp\t" + chr+ ":" + to_string(pos) + "-" + to_string(pos + 1)+ '\t' +  var1 +  '\t' +  var2 +  "\n.\tmotif\tpos\tBindAff_Var1\tBindAff_Var2\tlog(div)\tpvalue_diffBindAff\n");
	//}
}
		/*vector<double> keys;
		//int numSNPs = 0;
		#pragma omp parallel num_threads(2)
		{
			#pragma omp single nowait
			{
				cout << "start store dbSNPs" << endl;
				dbSNPs = s.storeDbSNPs(); // read dbSNPs 
				cout << "end store dbSNPs" << endl;
			}
			#pragma omp single nowait
			{


				//Step 1: for given input SNPs determine MAF and create plot with R which shows the distribution
				cout << "start write histogram and read MAF_counter" << endl;

				ofstream histogram;
				string file = outputDir + "histogram.txt";
				histogram.open(file);
				histogram << "SNP\tMAF\n"; //write header of the file
				for (auto i : MAF){
					histogram << "input\t" << i << "\n"; 
				}
				histogram.close();
				//plot histogram
				bc.callHistogram(file, outputDir + "/histogram_inputSNPs.pdf", sourceDir);
		
				//Step: 2 determine background
				MAF_counter = s.splitMAFinBins(); // split original MAF distribution in bins
				for(auto& i : MAF_counter){ //determine keys of the unordered map
					keys.push_back(i.first);
					//numSNPs += i.second;
				}
				cout << "MAF" << endl;
				for(auto& i : MAF_counter){
					cout << i.first << "\t" << i.second << endl;
				}
				cout << "end write histogram and read MAF_counter" << endl;
			}
		}*/
