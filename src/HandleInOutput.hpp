#ifndef HANDLEINOUTPUT_HPP
#define HANDLEINOUTPUT_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <iostream> 
#include <ostream>
#include <fstream>
#include <unistd.h>
#include <ctime> //for time
#include <chrono> //for time
#include <getopt.h> //einlesen der argumente
#include <unordered_map>
#include <unordered_set>
#include <random>

// include own functions
#include "callBashCommand.hpp"

using namespace std;

class InOutput{

	public: 
	//Constructor
	InOutput();
	//Deconstructor
	~InOutput();
	
	//functions
	void parseInputPara(int argc, char *argv[]);
	friend ostream& operator<< (ostream& os, InOutput& io);
	ofstream openFile(string path, bool app);
	void parseSNPsBedfile(string inputFile, int number);
	string getToken(string& line, char delim);
	void callHelp();
	int CountEntriesFirstLine(string inputFile, char delim);
	void checkIfSNPsAreUnique();
	void parseRandomSNPs(string inputFile, string REMsOverlappFile, string outputFile, int seed);
	void readScaleValues(string scaleFile, unordered_map<string, double>& scales);
	void getInfoSNPs(vector<string>& leadSNPs, unordered_map<string, vector<string>>& proxySNPs);
	void checkUniqAgain();
	int getNumberSNPs(string inputFile);
	void fileFormatVCF();

	//getter
	double getPvalue();
	double getPvalueDiff();
	string getFrequence();
	string getFootprints();
	//string getMutatedSequences();
	bool  getMaxOutput();
	string getActiveTFs();
	string getREMs();
	string getOutputAll();
	string getOutputDir();
	string getPFMs();
	string getSNPs();
	string getOverlappingFootprints();
	string getOverlappingREMs();
	string getSNPBedFile();
	string getSNPfastaFile();
	string getInfoFile();
	string getPFMsDir();
	string getEnsembleIDGeneName();
	double getActivityThreshold();
	//string getSourceDir();
	string getBedFileInDels();
	string getMappingGeneNames();
	string getResultFile();
	string getNotConsideredSNPs();
	string getGenome();
	string getdbSNPs();
	int getRounds();
	string getCodingRegions();
	int getConsideredSNPs();
	string getBackgroundSequences();
	int getNumberThreads();
	int getSeed();
	int getMinTFCount();
	string getScaleFile();
	bool getGeneBackground();
	bool getTfBackground();
	string getTransitionMatrix();
	string getRandomSNPs();

	private: //glaube das sollte nicht private sein
	int num_threads = 1; //-n
	double pvalue = 0.5; //-p	
	double pvalue_diff = 0.01; //-c
	string frequence = ""; //-b
	string footprint = ""; //-f path to footprint file
	string outputDir = "SNEEP_output/"; //-o
	string allOutput =  ""; //-d 
	bool maxOutput =  false; // -m
	//new
	string activeTFs = ""; // -t path to geneExpression file
	string REMs = ""; // -r bed-like REM file
	string mappingGeneNames = "";
	string fasta = ""; //stores fasta seq of snps
	string bed = ""; //bed file snps 
	string bed_notUniq = ""; // bed file not uniq after REM or footprint intersection
	string bed_notUniq_sorted = ""; //bed file not uniq after REM or footprint intersction but sorted
	string info = "";
	string overlappingFootprints = ""; //path to file where the overlapping footprints are stored
	string overlappingREMs = ""; // same for REMs
	string PFMsDir = ""; // output dir wherer each PFM is stored seperatly
	string bedFileInDels = "";
	string resultFile =  "";
	string PFMs; //must be given, transfac file
	string snpsNotUnique; // must be given
	string snpsNotUniqueSorted; // since inplace sort is not possible we need an addtional file here
	string snps = "";
	string notConsideredSNPs = "";
	string ensembleGeneName = "";
	double thresholdTFActivity = 0.0;
	string transition_matrix = "";
	//string sourceDir = "";
	//string genome = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/hg38.fa";
	//string genome = "/home/nbaumgarten/hg38.fa";
	string genome = "";
	int samplingRounds = 0;
	//string codingRegions = "";
	int consideredSNPs  = 0;
	string backgroundSeq = "";
	string dbSNPs = "";
	int seed = 1;
	int minTFCount = 0;
	bool geneBackground = false;
	bool tfBackground = false;
	string scaleFile = "";
	string randomSNPs = "";
};

//construtor
InOutput::InOutput()
{
//	cout << "constructor" << endl;
}
//deconstructor
InOutput::~InOutput()
{
}

void InOutput::parseInputPara(int argc, char *argv[]){
	
	int opt = 0;
	while ((opt = getopt(argc, argv, "o:n:p:c:b:af:mt:r:e:d:g:j:k:l:q:uvhx:i:")) != -1) {
       		switch (opt) {
		case 'o':
			outputDir = optarg;
			cout << "-o outputDir: " << outputDir << endl;
			break;
		case 'n':
			num_threads = stoi(optarg);
			cout << "-t number threads " << num_threads << endl;
			break;
		case 'p':
			pvalue = stod(optarg);
			cout << "-p use pvalue: " << pvalue <<  endl;
			break;
		case 'c':
			pvalue_diff = stod(optarg);
			cout << "-c use pvalue_diff: " << pvalue_diff <<  endl;
			break;
		case 'b':
			frequence = optarg;
			cout << "-b frequency: " << frequence << endl;
			break;
		case 'a':
			allOutput = outputDir + "AllDiffBindAffinity.txt";
			cout << "-a AllDiffBindAff: " << allOutput << endl;
			break;
		case 'f':
			footprint = optarg;
			cout << "-f footprint/region file: " << footprint << endl;
			break;
		case 'm':
			maxOutput = true;
			cout << "-m maxDiffBindAffinity: " << endl;
			break;
		case 't':
			activeTFs = optarg;
			cout << "-t activeTFs: " << activeTFs << endl;
			break;
		case 'r':
			REMs = optarg;
			cout << "-r REMs: " << REMs << endl;
			break;
		case 'e':
			ensembleGeneName = optarg;
			cout << "-e ensemble_geneName: " << ensembleGeneName << endl;
			break;
		case 'd':
			thresholdTFActivity = stod(optarg);
			cout << "-d threshold TF activity: " << thresholdTFActivity << endl;
			break;
		case 'g':
			mappingGeneNames = optarg;
			cout << "-g ensemblID to GeneName mapping: " << mappingGeneNames << endl;
			break;
		case 'l':
			seed = stoi(optarg);
			cout << "-l seed: " << seed << endl;
			break;
		case 'j':
			samplingRounds = stoi(optarg);
			cout << "-j number of randmoly sampled backgrounds: " << samplingRounds << endl;
			break;
		case 'k':
			dbSNPs = optarg;
			cout << "-k path to dbSNPs: " << dbSNPs << endl;
			break;
		case 'q':
			minTFCount = stoi(optarg);
			cout << "-q min TF count: " << minTFCount << endl;	
			break;
		case 'u':
			geneBackground = true;
			cout << "-u perform gene background analysis: " << geneBackground << endl;	
			break;
		case 'v':
			tfBackground = true;
			cout << "-v perform TF enrichment  analysis: " << tfBackground << endl;	
			break;
		case 'i':
			randomSNPs = optarg;
			cout << "-i RandomSNPs are given: " << randomSNPs << endl;	
			break;
		case 'x':
			transition_matrix = optarg;
			cout << "-w transition matrix: " << transition_matrix << endl;
			break;
		case 'h':
			callHelp();
			throw invalid_argument("help function end"); 
			//break;

		default:
			throw invalid_argument("invalid input parameter");
        	}
    	}

	fasta = outputDir + "snpRegions.fa"; //stores fasta seq of snps
	bed = outputDir + "snpRegions.bed"; //bed file snps 
	bed_notUniq = outputDir + "snpsRegions_notUniq.bed"; // bed file but noy uniq 
	bed_notUniq_sorted = outputDir + "snpsRegions_notUniq_sorted.bed"; // bed file but noy uniq but sorted
	info = outputDir + "info.txt";
	notConsideredSNPs = outputDir + "notConsideredSNPs.txt";
	overlappingREMs = outputDir + "overlappingREMs.bed"; // same for REMs
	PFMsDir = outputDir + "PFMs/"; // output dir wherer each PFM is stored seperatly
	//allOutput =  outputDir + "AllDiffBindAffinity.txt"; //-d 
	//maxOutput =  outputDir +  "MaxDiffBindingAffinity.txt"; // -m
	bedFileInDels = outputDir + "InDels.bed";
	resultFile = outputDir + "result.txt";
	if (footprint != ""){
		overlappingFootprints = outputDir + "SNPsOverlappingFootprints.bed";
	}
	if (REMs != ""){
		overlappingREMs = outputDir + "SNPsOverlappingREMs.bed";
	}
	if (samplingRounds > 0){
		backgroundSeq = outputDir + "backgroundSequences.bed";
	}
	if (optind + 4 > argc)  // there should be 4 more non-option arguments
		throw invalid_argument("missing motif file in TRANSFAC format, bed-like SNP file, genome file and/or scaleFile \n for help use  -h"); // TODO: besser in motif file umwandeln
	if (samplingRounds > 0 and dbSNPs.length() == 0){
		throw invalid_argument("for a background analysis the path to the sorted dbSNP file is requiered\n for help use -h"); // both parameters need to be set
	}
	if (activeTFs.length() > 0 and (thresholdTFActivity == 0.0  or ensembleGeneName.length() == 0 )){
		throw invalid_argument("either -d or -e is not set but requiered\n for help -h"); 
	} 
	if (samplingRounds == 0 and (tfBackground ||  geneBackground)){
		throw invalid_argument("number of backgroudn rounds -j (and dbSNP file -k) musst be specificed\n for help use -h"); 
	}

	PFMs = argv[optind++];
	cout <<"PFM dir: " << PFMs << endl;	
	snpsNotUnique = argv[optind++];
	snpsNotUniqueSorted = outputDir + "sortedSNPsNotUnique.txt";  
	snps = outputDir + "SNPsUnique.bed";
	cout <<"SNP file: " << snpsNotUnique << endl;
	genome = argv[optind++];
	cout << "genome file: " << genome << endl;
	scaleFile = argv[optind++];
	cout << "scale file: " << scaleFile << endl;



}

ostream& operator<< (ostream& os, InOutput& io){

	//determine time
        auto time = chrono::system_clock::now();
        time_t end_time = chrono::system_clock::to_time_t(time);

	os << "#\tdate and time: " << ctime(&end_time) << 
	"#\t-o outputDir: " << io.outputDir << 
	"\n#\t-p p-value threshold motifHits: " << io.pvalue << 
	"\n#\t-c p-value threshold diffBindAff: " << io.pvalue_diff << 
	"\n#\t-b file of background freq: " << io.frequence << 
	"\n#\t-f footprint/region file: " << io.footprint << 
	"\n#\t-m maxOutput: " << 
	"\n#\t-t activeTFs: " << io.activeTFs << 
	"\n#\t-r REMs: " << io.REMs <<
	"\n#\t-a allDiffBindAffinities: " << io.allOutput <<  
	"\n#\t-n number threads: " << io.num_threads << 
	"\n#\t-e ensemblID geneName mapping TFs: " << io.ensembleGeneName << 
	"\n#\t-d threshold TF activity: " << io.thresholdTFActivity <<
	"\n#\t-g EnsemblID to GeneName mapping REMs: " << io.mappingGeneNames << 
	"\n#\t-j rounds of background sampling: " << io.samplingRounds <<
	"\n#\t-k path to dbSNPs: " << io.dbSNPs <<
	"\n#\t-l start seed for random sampling: " << io.seed <<
	"\n#\t-q min TF count: " << io.minTFCount << 
	"\n#\t-u perform gene background analysis: " << io.geneBackground <<	
	"\n#\t-v perform TF enrichment  analysis: " << io.tfBackground <<	
	"\n#\t-x transition matrix for binding affinity p-value: " << io.transition_matrix <<	
	"\n#\tPFMs: " << io.PFMs << 
	"\n#\tSNPs file: " << io.snpsNotUnique << 
	"\n#\tpath to genome: " << io.genome <<
	"\n#\tpath to scale file: " << io.scaleFile <<
	"\n#\tinfo file: " << io.info;
	return os;
}

ofstream InOutput::openFile(string path, bool app){

        ofstream output;
        if (path != ""){
		//if (app){
                output.open(path,std::ios_base::app);
		//}else{
		//	output.open(path);
		//}
                if(!output.is_open()){
                        throw invalid_argument ("cannot open file:" + path);
                }   
        }   
        return output;
}

//gets the number of SNPs in the fastaa file 
// necessary to count them here, since we can lose snps which are out of the region of a chromosome if the genome version does not map to the dbSNP version

int InOutput::getNumberSNPs(string inputFile){


	ifstream input(inputFile); //either overlappingPeak file or snpFile
	int counter = 0;
	string line = "";
	while (getline(input, line, '\n')){
		counter++;
	}
	return counter/2;
}


// check if the snp file is in vcf format and if so parse in bed-like format
void InOutput::fileFormatVCF(){

	// is snp file in bed-like format or VCF format?
	string formatedSNPFile =  outputDir + "formatedSNPs.txt";
	string helper = snpsNotUnique; // copy of the file name
	char delim = '.';
	// check if file ending is vcf
	getToken(helper, delim); // split file name at point 
	//cout << helper << endl;
	if (helper == "vcf" or helper == "VCF"){
		//string formatedSNPFile =  outputDir + "formatedSNPs.txt";
		// call python script to parse vcf to bed-like format
		BashCommand bc_(getGenome()); //constructor bashcommand class
		bc_.callFormatingScript(snpsNotUnique, formatedSNPFile);
		snpsNotUnique = formatedSNPFile;
	}
	cout << "SNP file after VCF check: " << snpsNotUnique << endl;

}

//TODO:kommentieren
void InOutput::parseSNPsBedfile(string inputFile, int entriesSNPFile){

	double counterOverlappingREMs = 0;
	double counterOverlappingPeaks = 0;
	ofstream info_ = openFile(getInfoFile(), true);

	ifstream input(inputFile); //either overlappingPeak file or snpFile
	bool REMs = false;
	unordered_map<string, vector<string>> infoREMs; //chr:start:end_var1_var2 -> chr:start-end of REM, linked gene ensembl id, gene name, activity, coefficient
	char delim = '\t';
	int start = 0, end = 0, pos = 0;
	string chr = "", var1 = "", var2 = "", token = "", line = "", ensembl = ""; //stores current line
	if (getREMs() != ""){
		REMs = true;
		ifstream overlappingREMs(getOverlappingREMs());
		int entries = CountEntriesFirstLine(getREMs(),'\t');
		//read mapping ensembl id -> gene Name
		unordered_map<string, string> mappingGeneNames;
		ifstream mapping(getMappingGeneNames());
		while (getline(mapping, line, '\n')){
			ensembl = getToken(line, ',');
			mappingGeneNames[ensembl] = line;
		} 
		while (getline(overlappingREMs, line, '\n')){
			line = line + '\n';
			vector<string> helper; //stores info per line
			string chr = getToken(line, delim);
			string start = getToken(line, delim);
			string end = getToken(line, delim);
			string comb = chr + ":" + start + "-" + end;
//			cout << "comb: " << comb << endl;
			helper.push_back(comb);
			for (int i = 3; i <= entries-1; ++i){ //read only entries REMs
				pos = line.find(delim);
				var2 = line.substr(0, pos);
				helper.push_back(var2); //REM,ensembl id, geneName,  REMid,coefficient,pvalue, normModelScore, meanDNase1Signal, stdDNase1Signal,consortium, version
				if (i == 3){//determine gene name
					helper.push_back(mappingGeneNames[var2]);
				}
				line.erase(0, pos + 1);
			}
		/*	for(auto& elem : helper){
				cout << elem << " ";
			}
			cout << "\n";*/
			//read SNP Info
			chr =  getToken(line,delim);
			start =  getToken(line,delim);
			end =  getToken(line,delim);
			string var1 =  getToken(line,delim);
			string var2 =  getToken(line,delim);
			//string id = getToken(line,delim);
			//string MAF = getToken(line,delim);
			string key = chr + ":" + start + "-" + end + "_" + var1 + "_" + var2;
			//cout << "key REMs: " << key << endl;

			if (infoREMs.count(key)>0){ //key exists
				vector<string> existingInfo = infoREMs[key];
				for(int i = 0; i < existingInfo.size(); ++i){ //add new info to existing entries
					existingInfo[i] = existingInfo[i] + "," + helper[i]; 
				}
			/*	cout << "existing Info" << endl;
				for(auto& elem : existingInfo){
					cout << elem << " ";
				}
				cout << "\n"; */
				infoREMs[key] = existingInfo; //update key info
			}else{
				infoREMs[key] = helper;
			}
		}
		overlappingREMs.close();
	}
	ofstream output(bed_notUniq);
	ofstream output2(getBedFileInDels());//store InDels

	string id = "", MAF = "";

	while (getline(input, line, '\n')){
		line = line + '\n';	
		//getline(input, line, '\n'); //getLine
		//extract information
		chr = getToken(line, delim);
		//cout << "chr: " << chr << endl;
		start = stoi(getToken(line, delim)) - 50;
		//cout << "start: " << start << endl;
		end = stoi(getToken(line, delim)) + 50;
		//cout << "end: " << end << endl;
		var1 = getToken(line, delim);
		///cout << "var1: " << var1 << endl;
		var2 = getToken(line, delim);
		id = getToken(line,delim);
		MAF = getToken(line,delim);
		//cout << "var2: " << var2 << endl;
		string key = chr + ":" + to_string(start+50) + "-" + to_string(end-50) + "_" + var1 + "_" + var2;
		//cout << "key: " << key << endl;
		if ((var1 != "*") and (var2 != "*") and (var1.length() == 1) and (var2.length() == 1)){ 
			counterOverlappingPeaks++;
			output << chr << '\t' <<  start << '\t' << end << '\t' << chr << ":" << to_string(start+50) << "-"<< to_string(end-50) << ";" << var1 << ";" << var2 << ";" << id << ";" << MAF;// << '\t' << var1<< '\t' <<  var2 << '\n';
			//store skipped entries as header of the fasta file
			if (getOverlappingFootprints() == inputFile){
				vector<string> entriesLine;
				/*for(int i = 5; i <= entriesSNPFile-1; ++i){ //skip entries which are not important
					pos = line.find(delim);
					line.erase(0, pos + 1);
				}*/
				//read peak info
	//			cout << line << endl;
	//			cout << "raed peak info before" << endl;
				output << ";" << getToken(line, delim) << ":" << getToken(line, delim) << "-" << getToken(line, delim);
	//			cout << "raed peak info after" << endl;
			}else{
				output << ";.";
			}
			// add REM info
			if (REMs){
				if (infoREMs.count(key) == 0){
					output << ";.;.;.;.;.;.;.;.;.;.;."; 
					//fuer dennis Hi-c file
					//output << ";.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;."; 
				}else{
					counterOverlappingREMs++;
					vector<string> helper = infoREMs[key];
					for(auto& elem : helper){
				//		cout << elem << " ";
						output <<  ";" << elem;
					}
					//cout << '\n';
				}
			}
			output << '\n';
	
		}else{//write InDels in a file, without specific header
			output2 << chr << '\t' <<  start << '\t' << end << '\t' << chr << ":" << start << "," << var1 << "," << var2;// << '\t' << var1<< '\t' <<  var2 << '\n';
		}
	}
	if (getOverlappingFootprints() == inputFile){
		info_ << "!\toverlapPeak: " << counterOverlappingPeaks << '\n'; //info file
	}else{
		info_ << "!\toverlapPeak: -\n"; //info file
	}
	consideredSNPs = counterOverlappingPeaks;
//	cout << "peaks considered: " << consideredSNPs << endl;

	if (getREMs() != ""){
		info_ << "!\toverlapREMAllSNPs: " << counterOverlappingREMs << '\n';
	}else{
		info_ << "!\toverlapREMAllSNPs: -\n";
	}
	info_.close();
	input.close();
	output.close();
	output2.close();
	return;
}

string InOutput::getToken(string& line, char delim){
	int pos = 0;
	//cout << line.find(delim) << endl;
	if(((pos = line.find(delim)) != std::string::npos) || ((pos = line.find('\n')) != std::string::npos)){
    		string token = line.substr(0, pos);
    		line.erase(0, pos + 1);
		return token;
	}else{
		throw invalid_argument ("invalid file format:" + line);
	}
}

void InOutput::callHelp(){

	cout << "Call program with ./src/differentialBindingAffinity_multipleSNPs\noptinal parameters:\n" << 
	"-o outputDir (default SNEEP_output/, if you want to specific it, it must be done as first argument)\n" <<
	"-n number threads (default 1)\n" <<
	"-p pvalue for motif hits (default 0.5)\n"<<
	"-c pvalue differential binding (default 0.01)\n" <<
	"-b base frequency for PFMs -> PWMs ( /necessaryInputFiles/frequency.txt)\n" <<
	"-a if flag is set,  all computed differential bindinding affinities are stored in <outputDir>/AllDiffBindAffinity.txt\n"<<
	"-f additional footprint/open chromatin region file in bed file format\n" <<
	"-m if flag is set, the  maximal differential binding affinity per SNP is printed\n"<<
	"-t file where expression values of TFs are stored (e.g RNA-seq in a tab-seperated format e.g. ensemblID\texpression-value)\n" <<
	"-d threshold TF activity (must be given if -t is given)\n"<<
	"-e tab-seperated file containing ensemblID to gene name mapping of the TFs (must be given if -t is given)\n"<<
	"-r bed-like file with epigenetic interactions\n"<<
	"-g path to file containing ensemblID to gene name mapping, must be given if -r is given (,-seperated)(mapping for all genes within EpiRegio)\n" <<
	"-j rounds sampled background (default 0)\n"
	"-k path to sorted dbSNP file (if our provided file is used only SNPs in coding regions are considered)\n"
	//"-i path to the source GitHub dir (default .)\n"<<
	"-l start seed (default 1)\n" <<
	"-q minimal TF count which needs to be exceeded to be considered in random sampling (default 0)\n" << 
	"-u gene background analysis is performed (defaul false), -j must be set \n" <<
	"-v perform TF enrichment  analysis (default  false), -j must be set\n" <<	
	"-x transition matrix for binding affinity p-value, (default all transitions are equally likely) (necessaryInputFiles/transitionMatrix.txt)" <<
	"-h help\n" <<
	"transfac PFM file,  bed-like SNP file, path to genome file (fasta format) and scale file (see necessaryInputFiles/estimatedScalesPerMotif_1.9.txt for human data)  must be given"<<endl;
}

int InOutput::CountEntriesFirstLine(string inputFile, char delim){

	ifstream input(inputFile); //open file
	
	int pos = 0, entries = 0;
	string line = "";
	getline(input, line, '\n');
	while ((pos = line.find(delim)) != std::string::npos) {
		entries++;
		line.erase(0, pos + 1);
	}
	//add last elem which ends with \n and not with delim
	entries++;
	//cout << "entries: " << entries << endl;
	return entries;
}

void InOutput::getInfoSNPs(vector<string>& leadSNPs, unordered_map<string,vector<string>>& proxySNPs){

	ifstream input(snps); //open file
	string chr, start, end, allele1, allele2, rsID, maf, info, line, helper;
	//read file
	char delim = '\t';	
	string snp = "";
	int counter = 0;
	int c = 0;
	//vector<string> c_m;
	while (getline(input, line, '\n')){
		//cout << line << endl;
		chr =  getToken(line,delim);
		start =  getToken(line,delim);
		end =  getToken(line,delim);
		allele1 =  getToken(line,delim);
		allele2 =  getToken(line,delim);
		rsID =  getToken(line,delim);
		maf =  getToken(line,delim);
		info =  line;
		//if (rsID == "."){
		//	cout << chr << " " << start << " " << end <<  " " << rsID << endl;
		//}
		
		snp = chr + ":" + start + "-" + end + '\t' + allele1 + '\t' + allele2 + '\t' + rsID + '\t' + maf;
		
		if (info =="-"){ //identified lead snp
			leadSNPs.push_back(snp);

			//if (maf != "-"){
			//	c_m.push_back(rsID);
			//}
		}

		// store per lead snp the proxy snps
		
		else{ // the snp is a proxySNP
		// is there more than one lead SNPs associated to the current SNP?
			counter = count(info.begin(), info.end(), ',');
			c += counter;
			//cout << counter << " " << info << endl;
			//cout << counter << endl;
			for (int i=0; i<counter; ++i){
				//cout << i << endl;
				// split and add rsID
				helper = getToken(info, ',');
				if (proxySNPs.find(helper) != proxySNPs.end()){ // lead snp is already in map
					vector<string> a = proxySNPs[helper]; // add proxy SNPs to all lead SNPs
					a.push_back(rsID);
					proxySNPs[helper] = a;

				}else{
					vector<string> a {rsID};
					proxySNPs[helper] = a;
				}
			}
			//last rsID or of no , is found
			if (proxySNPs.find(info) != proxySNPs.end()){ // lead snp is already in map
				vector<string> a = proxySNPs[info]; // add proxy SNPs to all lead SNPs
				a.push_back(rsID);
				proxySNPs[info] = a;

			}else{
				vector<string> a {rsID};
				proxySNPs[info] = a;
			}
		}
	}
}

// check if snps are still uniq after intersection with REMs or footprints
void InOutput::checkUniqAgain(){

	BashCommand bc_(getGenome()); //constructor bashcommand class
	bc_.sort("-k1,1 -k2,2n", bed_notUniq, bed_notUniq_sorted);
	bc_.uniq("", bed_notUniq_sorted ,getSNPBedFile());
}

// to avoid to hold the file in memory, we rather sort the file and can so directly identify duplicated entries 
void InOutput::checkIfSNPsAreUnique(){

	// sort the file using a bash command
	BashCommand bc_(getGenome()); //constructor bashcommand class
	bc_.sort("-k1,1 -k2,2n", snpsNotUnique, snpsNotUniqueSorted);
	cout << "sorting done" << endl;

	ifstream input(snpsNotUniqueSorted); //open file
	ofstream output = openFile(snps, false);
	ofstream info_ = openFile(getInfoFile(), true);
	ofstream notConsidered = openFile(getNotConsideredSNPs(), false);

	double counterAllSNPs = 0; //counts all SNPs of the input file
	double counterUnique = 0; //counts relevant SNPs
	string line = "", originalLine = "", token = "";
	char delim = '\t';
	//vector<string> SNPs;
	int pos = 0;
	//read file
	string previous_pos = "000";
	unordered_set<string> seen_alleles;
	while (getline(input, line, '\n')){
		counterAllSNPs++;
		originalLine = line;
		// get info from file (chr start end allele1 allele2
		pos = line.find(delim);
		string chr_ = line.substr(0, pos);
		line.erase(0, pos + 1);
		pos = line.find(delim);
		string start = line.substr(0, pos);
		line.erase(0, pos + 1);
		pos = line.find(delim);
		string end = line.substr(0, pos);
		line.erase(0, pos + 1);
	
		string current_pos =  chr_ + ":" + start + "-" + end;

		pos = line.find(delim);
		string a1 = line.substr(0, pos);
		line.erase(0, pos + 1);
		pos = line.find(delim);
		string a2 = line.substr(0, pos);
		line.erase(0, pos + 1);
		if((a1 == "A"  or a1 == "C"  or a1 == "G"  or a1 == "T"  or a1 == "a"  or a1 == "c"  or a1 == "g"  or a1 == "t") and (a2 == "A"  or a2 == "C"  or a2 == "G"  or a2 == "T"  or a2 == "a"  or a2 == "c"  or a2 == "g"  or a2 == "t")){ //check if valid SNP
			if (current_pos == previous_pos){ // only if the positions are the same check if alleles are the same
				string alleles = a1 + "-" + a2;

				if (count(seen_alleles.begin(), seen_alleles.end(), alleles) == 0){
					//SNPs.push_back(helper);
					seen_alleles.insert(alleles);
					output << originalLine << '\n';
					counterUnique++;
				}
			}else{
				previous_pos = current_pos;
				seen_alleles.clear();
				output << originalLine << '\n';
				counterUnique++;
			}
		}else{
			notConsidered << originalLine + '\n';

		}
	}
	info_ << "!\tnumSNPSInputFile: " << counterAllSNPs << "\n" << "!\tnumRemovedDuplicates: " << counterAllSNPs - counterUnique << '\n'; 
	notConsidered.close();
	input.close();
	output.close();
	info_.close();
}


//TODO:kommentieren
void InOutput::parseRandomSNPs(string inputFile, string REMsOverlappFile, string outputFile, int seed){

//	double counterOverlappingREMs = 0;
//	double counterOverlappingPeaks = 0;
//	ofstream info_ = openFile(getInfoFile(), true);

	ifstream input(inputFile); //either overlappingPeak file or snpFile
	bool REMs = false;
	unordered_map<string, vector<string>> infoREMs; //chr:start:end_var1_var2 -> chr:start-end of REM, linked gene ensembl id, gene name, activity, coefficient
	char delim = '\t';
	int start = 0, end = 0, pos = 0;
	string chr = "", var1 = "", var2 = "", token = "", line = "", ensembl = ""; //stores current line
	//cout << "before doing anything" << endl;
	if (getREMs() != ""){
		int entries = 12; //TODO: stimmt das?
		REMs = true;
		ifstream overlappingREMs(REMsOverlappFile);
		//int entries = CountEntriesFirstLine(getREMs(),'\t');
		//read mapping ensembl id -> gene Name
		unordered_map<string, string> mappingGeneNames;
		ifstream mapping(getMappingGeneNames());
		while (getline(mapping, line, '\n')){
			ensembl = getToken(line, ',');
			mappingGeneNames[ensembl] = line;
		} 
		while (getline(overlappingREMs, line, '\n')){
			line = line + '\n';
			vector<string> helper; //stores info per line
			string chr = getToken(line, delim);
			string start = getToken(line, delim);
			string end = getToken(line, delim);
			string comb = chr + ":" + start + "-" + end;
//			cout << "comb: " << comb << endl;
			helper.push_back(comb);
			for (int i = 3; i <= entries-1; ++i){ //read only entries REMs
				pos = line.find(delim);
				var2 = line.substr(0, pos);
				helper.push_back(var2); //REM,ensembl id, geneName,  REMid,coefficient,pvalue, normModelScore, meanDNase1Signal, stdDNase1Signal,consortium, version
				if (i == 3){//determine gene name
					helper.push_back(mappingGeneNames[var2]);
				}
				line.erase(0, pos + 1);
			}
		/*	for(auto& elem : helper){
				cout << elem << " ";
			}
			cout << "\n";*/
			//read SNP Info
			chr =  getToken(line,delim);
			start =  getToken(line,delim);
			end =  getToken(line,delim);
			string var1 =  getToken(line,delim);
			string var2 =  getToken(line,delim);
			//string id = getToken(line,delim);
			//string MAF = getToken(line,delim);
			string key = chr + ":" + start + "-" + end + "_" + var1 + "_" + var2;
			//cout << "key REMs: " << key << endl;

			if (infoREMs.count(key)>0){ //key exists
				vector<string> existingInfo = infoREMs[key];
				for(int i = 0; i < existingInfo.size(); ++i){ //add new info to existing entries
					existingInfo[i] = existingInfo[i] + "," + helper[i]; 
				}
			/*	cout << "existing Info" << endl;
				for(auto& elem : existingInfo){
					cout << elem << " ";
				}
				cout << "\n"; */
				infoREMs[key] = existingInfo; //update key info
			}else{
				infoREMs[key] = helper;
			}
		}
		overlappingREMs.close();
	}
	//cout << "done withe REM file" << endl;
	ofstream output(outputFile);
	//ofstream output2(getBedFileInDels());//store InDels
	string id = "", MAF = "";
	int counter_commas = 0;
	mt19937 generator(seed); // seed muss fuer jeden thread ein andere sein 	
	int randomNum = 0; //sampled unifrom distributed number
	int counterSNPs = 0;

	while (getline(input, line, '\n')){
		counterSNPs++;
		line = line + '\n';	
		//getline(input, line, '\n'); //getLine
		//extract information
		chr = getToken(line, delim);
		//cout << "chr: " << chr << endl;
		start = stoi(getToken(line, delim)) - 50;
		//cout << "start: " << start << endl;
		end = stoi(getToken(line, delim)) + 50;
		//cout << "end: " << end << endl;
		var1 = getToken(line, delim);
	//	cout << "var1: " << var1 << endl;
		var2 = getToken(line, delim);
		id = getToken(line,delim);
		MAF = getToken(line,delim);
		//cout << "var2: " << var2 << endl;
		string key = chr + ":" + to_string(start+50) + "-" + to_string(end-50) + "_" + var1 + "_" + var2;
		//cout << "key: " << key << endl;

		//if there are multiple options for the mutant base, pick randomly one
		if (var2.length() != 1){
			counter_commas = count(var2.begin(), var2.end(), ',');
			vector<string> helper; 
			for (int i = 0; i < counter_commas; i++){
				helper.push_back(getToken(var2, ','));
			}
			helper.push_back(var2); //add last element (without a comma at the end)
			uniform_int_distribution<int> distribution(0, counter_commas); //specifiy distribution of the random number
			randomNum = distribution(generator); // generat random number

			var2 = helper[randomNum];
		//	if (var2 == "N" || var2 == "n"){
		//		var2 = helper[0]; 
		//	}
		}

		//if ((var1 != "*") and (var2 != "*") and (var1.length() == 1) and (var2.length() == 1)){ 
		output << chr << '\t' <<  start << '\t' << end << '\t' << chr << ":" << to_string(start+50) << "-"<< to_string(end-50) << ";" << var1 << ";" << var2 << ";" << id << ";" << MAF << ";.";// << '\t' << var1<< '\t' <<  var2 << '\n';
		// add REM info
		if (REMs){
			if (infoREMs.count(key) == 0){
				output << ";.;.;.;.;.;.;.;.;.;.;."; 
				//fuer dennis Hi-c file
				//output << ";.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;.;."; 
			}else{
				vector<string> helper = infoREMs[key];
				for(auto& elem : helper){
			//		cout << elem << " ";
					output <<  ";" << elem;
				}
				//cout << '\n';
			}
		}
		output << '\n';
	}
	input.close();
	output.close();
	return;// counterSNPs;
}




void InOutput::readScaleValues(string scaleFile, unordered_map<string, double>& scales){

	string line = "";
	string motif = "";
	double scale = 0.0;
	string skip = "";
	ifstream input(scaleFile); //open scaleFile
	getline(input, line, '\n'); // skip header
	while (getline(input, line, '\n')){
		motif = getToken(line, '\t');
		skip = getToken(line, '\t'); // scale newton
		skip = getToken(line, '\t'); // MSE
		scale = stod(getToken(line, '\t'));
		scales[motif] = scale;
	}
	input.close();
	return;
}


//getter
double InOutput::getPvalue(){
	return this->pvalue;
}

double InOutput::getPvalueDiff(){
	return this->pvalue_diff;
}

string InOutput::getFrequence(){
	return this->frequence;
}

string InOutput::getFootprints(){
	return this->footprint;
}
string InOutput::getScaleFile(){
	return this->scaleFile;
}
bool InOutput::getMaxOutput(){
	return this->maxOutput;
}

string InOutput::getActiveTFs(){
	return this->activeTFs;
}

string InOutput::getREMs(){
	return this->REMs;
}

string InOutput::getOutputAll(){
	return this->allOutput;
}

string InOutput::getOutputDir(){
	return this->outputDir;
}

string InOutput::getPFMs(){
	return this->PFMs;
}

string InOutput::getSNPs(){
	return this->snps;
}

string InOutput::getOverlappingFootprints(){
	return this->overlappingFootprints;
}

string InOutput::getOverlappingREMs(){
	return this->overlappingREMs;
}

string InOutput::getSNPBedFile(){
	return this->bed;
}

string InOutput::getSNPfastaFile(){
	return this->fasta;
}
string InOutput::getInfoFile(){
	return this->info;
}
string InOutput::getPFMsDir(){
	return this->PFMsDir;
}
string InOutput::getEnsembleIDGeneName(){
	return this->ensembleGeneName;
}
double InOutput::getActivityThreshold(){
	return this->thresholdTFActivity;
}
string InOutput::getBedFileInDels(){
	return this->bedFileInDels;
}
string InOutput::getMappingGeneNames(){
	return this->mappingGeneNames;
}
string InOutput::getResultFile(){
	return this->resultFile;
}
string InOutput::getNotConsideredSNPs(){
	return this->notConsideredSNPs;
}
string InOutput::getGenome(){
	return this->genome;
}
int InOutput::getRounds(){
	return this->samplingRounds;
}
int InOutput::getConsideredSNPs(){
	return this->consideredSNPs;
}
int InOutput::getSeed(){
	return this->seed;
}
string InOutput::getBackgroundSequences(){
	return this->backgroundSeq;
}
string InOutput::getdbSNPs(){
	return this->dbSNPs;
}
int InOutput::getNumberThreads(){
	return this->num_threads;
}
int InOutput::getMinTFCount(){
	return this->minTFCount;
}
bool InOutput::getGeneBackground(){
	return this->geneBackground;
}
bool InOutput::getTfBackground(){
	return this->tfBackground;
}
string InOutput::getTransitionMatrix(){
	return this->transition_matrix;
}
string InOutput::getRandomSNPs(){
	return this->randomSNPs;
}
#endif/*HANDLEINOUTPUT_HPP*/
