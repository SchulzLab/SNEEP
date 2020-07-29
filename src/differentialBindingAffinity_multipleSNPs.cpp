#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <unordered_map>
#include <stdexcept>
#include <algorithm>
#include <omp.h>
#include <bitset>
#include <ctime> //for time
#include <chrono> //for time

//own classes
#include "pvalue_copy.hpp"
#include "Matrix_new.hpp"
#include "callBashCommand.hpp"
#include "HandleInOutput.hpp"


//neccessary for function strcmp
#include <stdio.h>
#include <string.h>

//important for reading a directory
#include <dirent.h>
#include <sys/types.h>
#include <getopt.h> //einlesen der argumente

//for log
#include <math.h>

using namespace std;

//------------------------------
//Global variables default values 
int COMPLEMENT[] = {0,4,0,3,1,0,0,2,0,0,0,0,0,0,5};// considers also N
int POSITION[] = {0,1,0,2,4,0,0,3,0,0,0,0,0,0,5}; // considers also N 
double EPSILON = 0.001;
int NUM_THREADS = 1; //number of thread the programm is running on
int MAXIMUM_LENGTH_SEQUENCE = 1000; // number of the maximum lenght of the input sequnces
string FREQUENCE = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/Phase2/frequence.txt";

//functions
vector<string> readDirectory(const char *path); //stores all files of a directory
vector<double> readFrequence(string frequence); //read the base pair frequence of the genome or 0.25 per base
vector<Matrix<double>> PFMsToPWMs( vector<string>& PWM_files, vector<double>& freq, string path_to_pwms); //determines PWMs based on the input count matrices
int probProSeq(Matrix<double>& PWM, string& line, vector<double>& pvalues, double pvalue, int& pos_snp, vector<double>& prob_sequences);//determines the binding affinity pro sequence
double probKmer(Matrix<double>& PWM,string::iterator start, string::iterator end); // determines binding affinity for each kmer
double probKmerComplement(Matrix<double>& PWM,string::iterator start, string::iterator end); // same for the complement
int getMutatedPos(string mut_pos); //calculates which position in the input sequence is the mutated position 
string getChr(string& mut_pos); //return the chr of  chr:start-end 
int getSeqStart(string& header); //same but return the start
int getSeqEnd(string& header); //return the end
string createMutatedSeq(string header, string seq, vector<string>& splittedHeader); //determines the mutated sequence based on the given SNP
double cdf_laplace(double location, double scale, double value); //determine laplace distribution
double differentialBindingAffinity(double& first_elem, double& second_elem, unordered_map<double, double>& pre_log, double& log_); //calculates differential binding affinity
string forwardOrReverse(int distance_); //determines on which strand the motif bindes
vector<string> parseHeader(string header, string delim); //parse the header of the fasta file
void writeHeadersOutputFiles(InOutput io, string header, string seq, string mutSeq, bool app); // write the headers to the output files

int main(int argc, char *argv[]){

	//handle input 
	InOutput io; //create object
	io.parseInputPara(argc, argv); //set all paths according to input 

	//create instance bashCommand
	BashCommand bc(io.getGenome()); //constructor bashcommand class
	bc.mkdir(io.getOutputDir(), "-p", true); //make outputDir
	ofstream info = io.openFile(io.getInfoFile(), false); //open info file
	info << io << endl; // write settings 
	info.close();

	io.checkIfSNPsAreUnique(); // removes SNPs from inputSNP list which are not unique  and stores them in info file
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
		resultFile << "SNP_position\tvar1\tvar2\tpeak_position\tTF\tTF-binding_position\tstrand\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tREM_positions\tensemblIDs\tgeneNames\tREMIds\tcoefficients\tpvalues_REM\tnormModelScore\tmeanDNase1Signal\tstdDNase1Signal\tconsortium\tversion\tSNP_strand\n";
		string output = io.getOverlappingREMs();
		//bc.intersect(REMs, SNPFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
		bc.intersect(REMs, SNPFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
	}else{
		
		//write header output file
		resultFile << "SNP_position\tvar1\tvar2\tpeak_position\tTF\tTF-binding_position\tstrand\tpvalue_BindAff_var1\tpvalue_BindAff_var2\tlog_pvalueBindAffVar1_pvalueBindAffVar2\tpvalue_DiffBindAff\tSNP_strand\n";
	}
	//determine overlapping peaks
	string footprintFile = io.getFootprints();
	if (footprintFile != ""){ //if footprint or region file is given
		string output = io.getOverlappingFootprints();
		//determine intersection
		//bc.intersect(SNPFile,footprintFile, output, "-wao"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
		bc.intersect(SNPFile,footprintFile, output, "-wa -wb"); //result stored in outputDir +  SNPsOverlappingFootrpints.bed
		//parse bedfile 
		io.parseSNPsBedfile(output, entriesSNPFile); //skip insertions that are longer than 1 and remember additionl number of entries per line (the first 5 entries are not counted, since they are necessary), result sored in OutputDir + snpRegions.bed, use getSNPBedFile() 

	}else{
		//parse bedfile 
		io.parseSNPsBedfile(SNPFile, entriesSNPFile);
	}
	//call getFasta
	bc.getFasta(io.getSNPBedFile(), io.getSNPfastaFile(), "-name");
	cout << "after getFastaa" << endl;

	//check which TFs are active and write this TFs seperated in files  (done with a python script) also determine frequence matrix from the count matrices
	string activeTFs = io.getActiveTFs();
	bc.mkdir(io.getPFMsDir(), "-p", false); //create dir
	bc.rm(io.getPFMsDir()); // remove PWMs if there are any
	if (activeTFs != ""){
		bc.callPythonScriptCheckActiveMotifs(io.getSourceDir(), io.getActiveTFs(), io.getPFMs(), io.getPFMsDir(), io.getEnsembleIDGeneName(), io.getActivityThreshold()); 
	}else{
		bc.callPythonScriptSplitPFMs(io.getSourceDir(), io.getPFMs(), io.getPFMsDir()); 
	}
	//convert PFMs internally to PWMs
	vector<string> PWM_files = readDirectory(io.getPFMsDir().c_str()); //stores names of all pwm files
	vector<double> freq = readFrequence(io.getFrequence()); //TODO datenspezifisch bestimmen gibts da ein ebdtool function?
	vector<Matrix<double>> PWMs = PFMsToPWMs(PWM_files,freq, io.getPFMsDir());


	//determine pvalues for PFMs 
	string motif = "";
	Matrix<double> transition_matrix(4,4);
	ifstream transition_file("necessaryInputFiles/transition_matrix.txt"); //TODO: datenspezifisch bestimmen
        transition_file>> transition_matrix;
//	cout << "transition matrix:\n" << transition_matrix << endl;
	unordered_map<string, vector<double>> all_pvalues;
	for(int j = 0; j < PWM_files.size(); ++j){ // iterate over all given motifs
		motif = PWM_files[j].substr(0 , PWM_files[j].size() -4); // set motif
		pvalue  pvalue_obj(PWMs[j].ncol(), EPSILON); // EPSILON entspricht accuracy !!! 
		vector<double> pvalues = pvalue_obj.calculatePvalues(PWMs[j], freq, transition_matrix);
		all_pvalues[motif] = pvalues;
	}

	//stores all sequences and the according header
	string line = "";
	vector<string> sequences, headers;
	unordered_map<double, double> preLog; //stores log of pvalues, avoid to recalculate them

	int sigEntries = 0; //counts significant entires 
	int sigEntriesOverlappingREM = 0;

	ifstream pro_seq(io.getSNPfastaFile()); //open fasta file
	while (!(pro_seq.eof())){
		getline(pro_seq, line, '\n');
		if (line[0] == '>')
			line = line.substr(1, line.size());
		else
			line = line.substr(0, line.size());
		headers.push_back(line);
		getline(pro_seq, line, '>'); //reads until next header (starts with >)	
		line.pop_back(); //delets new line symbol (\n)
		sequences.push_back(line);
	}
	pro_seq.close(); // close fasta file

	for (int i = 0; i <sequences.size(); ++i){ //iterates over all sequences which contain the SNP(each sequence: 50bp + SNP + 50bp)
//		cout << headers[i] << " " <<  sequences[i] << endl;

		//define variables
		vector<tuple<double, string>> overallResultMax; //stores max result
		unordered_map<string, string> helperOverallResult;
		string chr = "",  mutSeq = "", direction = ""; //stores if the SNP annotation was for the forward or the reversed strand
		int posMut = 0, pos = 0; // position of the snp within the genome an relative within the given sequence

		//create mutated sequence
		vector<string> splittedHeader; //snp_pos var1 var2 peak_pos REM_posensemblId genName REMId coefficient pvalue consortium version reverse
		mutSeq = createMutatedSeq(headers[i], sequences[i], splittedHeader); // check if all SNPs are fitting to the reference this is done within create mutated seq (also check for reverse strand)
		pos = getMutatedPos(splittedHeader[0]);
		chr = getChr(splittedHeader[0]);
		posMut = 50; // since we considere 50 bp upstream and downstream the snp
		
		if (!mutSeq.empty()){//if mutSeq is not a fit, the stringis empty
			//write parts of the output
			if (i == 0){
				writeHeadersOutputFiles(io, headers[i], sequences[i], mutSeq, false);
			}else{
				writeHeadersOutputFiles(io, headers[i], sequences[i], mutSeq, true);
			}

		//open MaxDiffBinding file and AllOutputFile
		ofstream output, outputMax;
		if (io.getOutputAll() != ""){
			output = io.openFile(io.getOutputAll(), true);
		}
		if (io.getMaxOutput() != ""){
			outputMax = io.openFile(io.getMaxOutput(), true);
		}
		
		for(int j = 0; j < PWM_files.size(); ++j){ // iterate over all given motifs
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
		
				motif = PWM_files[j].substr(0 , PWM_files[j].size() -4); //set motif
			//	cout << "motif: " << motif << endl;
	
				int lenMotif = PWMs[j].ncol();
//				cout << PWMs[j] << endl;
				vector<double> pvalues = all_pvalues[motif]; //pvalues for current motif

				//calculate probabiblity for sequences of the current TF
				if (probProSeq(PWMs[j], sequences[i],  pvalues, io.getPvalue(), posMut, firstSeq) == 1){ //first_seq contains wildtype sequence
						sigHit  = 1;
				}
				//calculate pvalue of the mutated sequence	
				if (probProSeq(PWMs[j], mutSeq, pvalues, io.getPvalue(), posMut, secondSeq) == 1){ //second_seq contain mutated sequence
					sigHit = 1;
				}
				//determine differential binding affinity for all sig kmers in wildtype and mutated seq
				for (int l = 0; l < firstSeq.size(); ++l){ 
			
					if (firstSeq[l] <= io.getPvalue() or secondSeq[l] <= io.getPvalue()){
						
						//posSigHit = seq_start + floor(l/2); //deterime position of the hit	
						posSigHit = pos + floor(l/2); //deterime position of the hit	

						//determine differntial binding
						diffBinding = differentialBindingAffinity(firstSeq[l], secondSeq[l], preLog, log_);
							
						//stores info maximal diff binding affinity
						if (diffBinding <= maxDiffBinding){
							maxDiffBinding = diffBinding;
							maxPosSigHit = posSigHit;
							maxOrientation = forwardOrReverse(l);
							maxVal1 = firstSeq[l];
							maxVal2 = secondSeq[l];
							maxLog = log_;
							maxLenMotif = lenMotif;
						}
						//write output to the file
						if (io.getOutputAll() != ""){
							output << motif << '\t' << chr << ":" << posSigHit << forwardOrReverse(l) << '\t' << firstSeq[l]<< "\t" << secondSeq[l]<< "\t" << log_ << '\t' <<  diffBinding<< '\n';
						}
					}
				}
				//ouput maximal binding affinity for the current seq
				overallResultMax.push_back(make_tuple(maxDiffBinding, motif));
				string value =  chr + ':' + to_string(maxPosSigHit- maxLenMotif+1) + "-" + to_string(maxPosSigHit+1) + '\t' +  maxOrientation + '\t' + to_string(maxVal1) + '\t' + to_string(maxVal2) + '\t' + to_string(maxLog) + '\t' +  to_string(maxDiffBinding);
				helperOverallResult[motif] = value;

			/*	cout << "pvalues of the first sequence:" << endl;
				for(auto& s: firstSeq)
					cout << s << "\t";
				cout << endl;
	
				cout << "pvalues of the second sequence:" << endl;
				for(auto& s: secondSeq)
					cout << s << "\t";
				cout << endl;*/
			}
	
			//sort overall_result per seq and all TFs
			string m = "";
			double p = 0.0;
			sort(overallResultMax.begin(), overallResultMax.end());
			//prints overall result max per sequence for all TFs
			//for(auto& s : overall_result_max)
			//	cout << get<0>(s) << " " << get<1>(s) << endl;
			for (int i = 0; i < overallResultMax.size(); ++i){
				p = get<0>(overallResultMax[i]); //pvalue
				m = get<1>(overallResultMax[i]); //motif name
				if (p <= io.getPvalueDiff()){
					sigEntries+=1;
					if (io.getMaxOutput() != ""){
						outputMax << i +1 << '\t' << m << '\t' << helperOverallResult[m] + '\n';	
					}
					//for (int k =0; k< splittedHeader.size();++k){
					if (REMs == ""){
						resultFile << splittedHeader[0] << '\t' << splittedHeader[1] << '\t' << splittedHeader[2] << '\t' << splittedHeader[3] << '\t' << m << '\t' <<  helperOverallResult[m] << '\t' << splittedHeader[4] << '\n';
					}else{
						if (splittedHeader[5] != "."){
							sigEntriesOverlappingREM += 1;
						}
						resultFile << splittedHeader[0] << '\t' << splittedHeader[1] << '\t' << splittedHeader[2] << '\t' << splittedHeader[3] << '\t' << m << '\t' <<  helperOverallResult[m] << '\t' << splittedHeader[4] << '\t' <<  splittedHeader[5] << '\t'  << splittedHeader[6] << '\t' << splittedHeader[7] << '\t' << splittedHeader[8] << '\t' <<  splittedHeader[9] << '\t' <<  splittedHeader[10] << '\t' << splittedHeader[11] << '\t' << splittedHeader[12]<< '\t'  << splittedHeader[13] << '\t' << splittedHeader[14] << "\t" << splittedHeader[15] << '\n'; 
					}
				}
			}
		if (io.getOutputAll() != ""){
			output.close();
		}
		if (io.getMaxOutput() != ""){
			outputMax << "\n";
			outputMax.close();
		}
		}
	}
	resultFile.close();

	//open info file
	info = io.openFile(io.getInfoFile(), true); //open info file

	info << "!\tsigEntryOverlappingREMs: " << sigEntriesOverlappingREM << endl; //write number of entries which are significant and overlap with a REM
	info << "!\tsigOnes: " << sigEntries << endl; // write number of entries which are significant

	//write number motifs to output file 
	info << "!\tnumMotifs: " << PWM_files.size() << endl; // number of used motisf

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
		ifstream PWM(path_to_pwms + PWM_files[k]);//open actual PWM file
		PWM >> PWM_matrix; //read PWM file as matrix
		// determine PWM, round matrix and shift it in such a way that only values > 0 are included in the matrix -> important for pvalue calculation
		for (int i = 1; i<= PWM_matrix.ncol(); i++){
			for (int j = 1; j <= 4; j++){
				value = log10(PWM_matrix(j,i)/freq[j-1]) - log10(EPSILON/freq_max) + EPSILON;
				PWM_matrix(j,i) =  (round(value * rounder ) / rounder);
			}
		}
		PWMs[k] = PWM_matrix; // store matrix in map
	}
	return PWMs;
}

int probProSeq(Matrix<double>& PWM, string& line, vector<double>& pvalues, double pvalue, int& pos_snp, vector<double>& pvalue_seq){

	int min_ = pos_snp -50;
	int max_ = pos_snp + 50;

	int length_motif = PWM.ncol();
	int start_seq = max(pos_snp - length_motif, min_-1);
	int end_seq = min(pos_snp, max_-1)-1;
//	cout <<"lengthMotif: " << length_motif << " start_seq: " << start_seq << " end seq: " << end_seq << endl;

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
		if (line[i] == 'N'){
			check = length_motif; 
		}
		if (line[i] != 'N' && check == 0 && start >= line.begin() && counter >= start_seq && counter <= end_seq){
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
//			cout << "prob: " << prob << " pos : " << pos << " pvalue: " << pvalue_prob << endl;
//			cout << "probComp: " << prob_comp << " posComp : " << pos_comp << " pvalueComp: " << pvalue_prob_comp << endl;
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
	
        double prob = 0.0; //muss 0.0 sein bei plus
	int counter = 1;
	for(auto i = start; i != end+1; i++){ //iterators over k_mer
//		cout << *i;
		prob+= PWM(POSITION[((*i) & 0x0f)], counter); //determine which letter we do consider and look up the accodirng prob in the PWM
		counter ++;
	}
//	cout << endl;
//	cout << prob << endl;
        return prob;
}
//same as for probKMer, only different iterates in reverse order and determne complement for actual letter
double probKmerComplement(Matrix<double>& PWM,string::iterator start, string::iterator end){

        double prob = 0.0; //muss 0.0 sein bei plus
	int counter = 1;
	for(auto i = end; i != start-1; i--){// iterator in reverse order
		//cout << COMPLEMENT[((*i) & 0x0f)];
//		cout << *i;
		prob+= PWM(COMPLEMENT[((*i) & 0x0f)], counter); //determine which letter we do consider and look up the accodirng prob in the PWM
		counter ++;
	}
//	cout << endl;
//	cout << prob << endl;
        return prob;
}

int getMutatedPos(string mut_pos){

	int del = mut_pos.find(":");
	int pos = stoi(mut_pos.substr(del+1));
//	cout << "pos: " << pos << endl;
	return pos;
}

string getChr(string& mut_pos){

	int del = mut_pos.find(":");
	string chr = mut_pos.substr(0, del);
//	cout << "chr: " << chr << endl;
	return chr;
}

int getSeqStart(string& header){

	int del = header.find(":");
	int del2 = header.find("-");
	int start = stoi(header.substr(del+ 1, del2-del));
//	cout << "start: " << start << endl;
	return start;
}

int getSeqEnd(string& header){

	int del = header.find(":");
	int del2 = header.find("-");
	int end = stoi(header.substr(del2+1));
//	cout << "end: " << end << endl;
	return end;
}

string createMutatedSeq(string header, string seq, vector<string>& splittedHeader){

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
		direction = "forward";
		result[mut] = var2;
	}else if (toupper(seq[mut]) == var2){ // forward strand var2 in reference seq
		direction = "forward";
		result[mut] = var1;
	//check for reverse strand
	}else if (toupper(seq[mut]) == complement[var1]){ //reverse strand var1 in reference seq
		direction = "reverse";
		result[mut] = complement[var2];
	}else if (toupper(seq[mut]) == complement[var2]){ //reverse strand var2 in refernece seq
		direction = "reverse";
		result[mut] = complement[var1];
	}else{
		cout << "letter at position 50 (0-based): " << seq[mut] << " "  << chr << ":" << pos << " seq: " << seq << " var1: " << var1 << " var2: " << var2 << " header " << header << endl;
	//	throw invalid_argument("mutation passt nicht!");
		cout << "da stimmt was nicht !!!!" << endl;
	}
	splittedHeader.push_back(direction);
	return result;
}

double cdf_laplace(double location, double scale, double value){

	double result = 0.0;
	if (value < location){
		result = 0.5 * (exp(((value - location)/scale)));
	}else{ //value >= location
		result = 1 - (0.5 * (exp(-((value - location)/scale))));
	}
	return result;
}

double differentialBindingAffinity(double& first_elem, double& second_elem, unordered_map<double, double>& pre_log, double& log_){

	double div = 0.0;
	log_ = 0.0;
	double pvalue = 0.0;

	div = first_elem/second_elem;

	if (pre_log.find(div) != pre_log.end()){ //check if the log of div was calculated bevore
		log_ = pre_log[div];
	}else{
		log_ = log(div);
		pre_log[div] = log_;		
	}

	//determine pvalue
	//symmetric laplace(0,1) 
	pvalue = 1 - cdf_laplace(0.0, 1.0, abs(log_));
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

void writeHeadersOutputFiles(InOutput io, string header, string seq, string mutSeq, bool app){
	
	ofstream output, output_max;
	if (app){
		if (io.getOutputAll() != "")
			output = io.openFile(io.getOutputAll(), true);
		if (io.getMaxOutput() != "")
			output_max = io.openFile(io.getMaxOutput(), true);
	}else{
		if (io.getOutputAll() != ""){
			output = io.openFile(io.getOutputAll(), false);
			output << io << endl;
		}
		if (io.getMaxOutput() != "")
			output_max = io.openFile(io.getMaxOutput(), false);
	}

	vector<string> splittedHeader = parseHeader(header, ";"); //chr:start, var1, var2 
	char var1 = toupper(splittedHeader[1][0]); //cast to char
	char var2 = toupper(splittedHeader[2][0]); //cast to char
	int pos = getMutatedPos(splittedHeader[0]);
	string chr = getChr(splittedHeader[0]);
	int start = pos - 50;
	int end = pos + 50;
	int mut = pos - start; //position of the mutation in relation to the sequence length (should be at position 51)

	if (io.getOverlappingFootprints() != ""){
		if (io.getOutputAll() != "")
			output <<"\n#\tfootprint info: "<<  header << endl; //TODO oder besser nur peak info
	}
	//write to seq file
	if (io.getOutputAll() != "")
		output <<"#\twildtyp: " <<  seq << endl;
	if (io.getMutatedSequences() != ""){
		ofstream output2;
		if (app){ 
			output2 = io.openFile(io.getMutatedSequences(), true);
		}else{
			output2 = io.openFile(io.getMutatedSequences(), false);
		}
		output2 <<">" << chr <<  ":" << pos  << " " << var1 << " " << var2 << "\nwil: " <<   seq <<"\nmut: " << mutSeq <<  endl;
		output2.close();
	}
	if (io.getOutputAll() != ""){
		output <<"#\tmutated: " <<  mutSeq << endl;
		output << "\n.\t";
		output << "pos\tBindAff_Var1\tBindAff_Var2\tlog(div)\tpvalue_diffBindAff\n";
		output << endl;
		output.close();
	}
	//header output_max
	if (io.getMaxOutput() != ""){
		output_max << "snp\t" << pos << '\t' << var1 << '\t'<< var2 << '\n';
		output_max << ".\tmotif\tpos\tBindAff_Var1\tBindAff_Var2\tlog(div)\tpvalue_diffBindAff\n";
		output_max.close();
	}
	
}
