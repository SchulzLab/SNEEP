#ifndef HANDLEINOUTPUT_HPP
#define HANDLEINOUTPUT_HPP

#include <string>
#include <vector>
#include <algorithm>
#include <cctype> //tolower
#include <iostream> 
#include <ostream>
#include <fstream>
#include <unistd.h>
#include <ctime> //for time
#include <chrono> //for time
#include <getopt.h> //parse the command line arguments
#include <unordered_map>
#include <unordered_set>
#include <random>

// include own functions
#include "callBashCommand.hpp"
#include "stringUtils.hpp"

// number of bases on each side of the SNV in the extracted sequences (sequence = SNV_FLANK bp + SNV + SNV_FLANK bp, SNV at index SNV_FLANK)
const int SNV_FLANK = 50;

using namespace std;

class InOutput{

	public: 
	//Constructor
	InOutput();
	//Destructor
	~InOutput();
	
	//functions
	void parseInputPara(int argc, char *argv[]);
	friend ostream& operator<< (ostream& os, InOutput& io);
	ofstream openFile(string path, bool app);
	void parseSNPsBedfile(string inputFile, int number);
	void callHelp();
	int CountEntriesFirstLine(string inputFile, char delim);
	void checkIfSNPsAreUnique();
	void parseRandomSNPs(string inputFile, string REMsOverlappFile, string outputFile, int seed);
	unordered_map<string, vector<string>> readOverlappingREMs(string overlapFile);
	string remColumns(unordered_map<string, vector<string>>& infoREMs, const string& key, int numREMFields);
	void readScaleValues(string scaleFile, unordered_map<string, double>& scales);
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
	int getNumberThreads();
	int getSeed();
	int getMinTFCount();
	string getScaleFile();
	bool getGeneBackground();
	bool getTfBackground();
	string getTransitionMatrix();
	string getRandomSNPs();
	bool getGCMatching();

	private:
	int num_threads = 1; //-n
	double pvalue = 0.5; //-p	
	double pvalue_diff = 0.01; //-c
	string frequence = ""; //-b
	string footprint = ""; //-f path to footprint file
	string outputDir = "SNEEP_output/"; //-o
	string allOutput =  ""; //-a 
	bool writeAllOutput = false; //-a, the path is set after parsing all arguments (depends on -o)
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
	string dbSNPs = "";
	int seed = 1;
	int minTFCount = 0;
	bool geneBackground = false;
	bool tfBackground = false;
	string scaleFile = "";
	string randomSNPs = "";
	bool gcMatching = false; //-s match GC content in the background sampling
	string gcMatchingWarning = ""; //reason why -s is ignored
};

//construtor
InOutput::InOutput()
{
//	cout << "constructor" << endl;
}
//destructor
InOutput::~InOutput()
{
}

/*
* value of a flag as int / double; the whole value must be a number, otherwise a clear error is thrown
*/
int toInt(const string& value, char flag){
	try{
		size_t pos = 0;
		int result = stoi(value, &pos);
		if (pos == value.size()){
			return result;
		}
	}catch (const exception& e){
	}
	throw invalid_argument("invalid value for -" + string(1, flag) + ": " + value + " (integer expected)");
}

double toDouble(const string& value, char flag){
	try{
		size_t pos = 0;
		double result = stod(value, &pos);
		if (pos == value.size()){
			return result;
		}
	}catch (const exception& e){
	}
	throw invalid_argument("invalid value for -" + string(1, flag) + ": " + value + " (number expected)");
}

void InOutput::parseInputPara(int argc, char *argv[]){
	
	int opt = 0;
	while ((opt = getopt(argc, argv, "o:n:p:c:b:af:mt:r:e:d:g:j:k:l:q:uvhx:i:s:")) != -1) {
       		switch (opt) {
		case 'o':
			outputDir = optarg;
			cout << "-o outputDir: " << outputDir << endl;
			break;
		case 'n':
			num_threads = toInt(optarg, 'n');
			cout << "-n number threads: " << num_threads << endl;
			break;
		case 'p':
			pvalue = toDouble(optarg, 'p');
			cout << "-p use pvalue: " << pvalue <<  endl;
			break;
		case 'c':
			pvalue_diff = toDouble(optarg, 'c');
			cout << "-c use pvalue_diff: " << pvalue_diff <<  endl;
			break;
		case 'b':
			frequence = optarg;
			cout << "-b frequency: " << frequence << endl;
			break;
		case 'a':
			writeAllOutput = true;
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
			thresholdTFActivity = toDouble(optarg, 'd');
			cout << "-d threshold TF activity: " << thresholdTFActivity << endl;
			break;
		case 'g':
			mappingGeneNames = optarg;
			cout << "-g ensemblID to GeneName mapping: " << mappingGeneNames << endl;
			break;
		case 'l':
			seed = toInt(optarg, 'l');
			cout << "-l seed: " << seed << endl;
			break;
		case 'j':
			samplingRounds = toInt(optarg, 'j');
			cout << "-j number of randmoly sampled backgrounds: " << samplingRounds << endl;
			break;
		case 'k':
			dbSNPs = optarg;
			cout << "-k path to dbSNPs: " << dbSNPs << endl;
			break;
		case 'q':
			minTFCount = toInt(optarg, 'q');
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
		case 's':{
			string value = optarg;
			transform(value.begin(), value.end(), value.begin(), ::tolower);
			if (value == "true" or value == "1"){
				gcMatching = true;
			}else if (value == "false" or value == "0"){
				gcMatching = false;
			}else{
				throw invalid_argument("-s must be true or false (or 1 or 0): " + string(optarg));
			}
			cout << "-s match GC content in background sampling: " << gcMatching << endl;
			break;
		}
		case 'x':
			transition_matrix = optarg;
			cout << "-x transition matrix: " << transition_matrix << endl;
			break;
		case 'h':
			callHelp();
			throw invalid_argument("help function end"); 
			//break;

		default:
			throw invalid_argument("invalid input parameter, for help use -h");
        	}
    	}

	// all output paths are set after parsing, so -o can be given at any position and with or without a final /
	if (outputDir.empty()){
		outputDir = "./";
	}
	if (outputDir.back() != '/'){
		outputDir += '/';
	}
	if (writeAllOutput){
		allOutput = outputDir + "AllDiffBindAffinity.txt";
		cout << "-a AllDiffBindAff: " << allOutput << endl;
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
	if (optind + 4 > argc)  // there should be 4 more non-option arguments
		throw invalid_argument("missing motif file in TRANSFAC format, bed-like SNP file, genome file and/or scaleFile \n for help use  -h"); // TODO: besser in motif file umwandeln
	if (REMs != ""){ // the REM columns of result.txt are fixed (see header), so the interaction file needs exactly 12 columns
		if (!ifstream(REMs)){
			throw invalid_argument("cannot open interaction file (-r): " + REMs);
		}
		int columns = CountEntriesFirstLine(REMs, '\t');
		if (columns != 12){
			throw invalid_argument("the interaction file (-r) " + REMs + " has " + to_string(columns) + " columns, but 12 tab-separated columns are expected (chr, start, end, ensemblID, regionID and 7 further columns, '.' or '-' if not available)");
		}
	}
	if (samplingRounds > 0 and dbSNPs.length() == 0 and randomSNPs.length() == 0){
		throw invalid_argument("for a background analysis the path to the sorted dbSNP file is requiered\n for help use -h"); // both parameters need to be set
	}
	if (activeTFs.length() > 0 and (thresholdTFActivity == 0.0  or ensembleGeneName.length() == 0 )){
		throw invalid_argument("either -d or -e is not set but requiered\n for help -h"); 
	} 
	if (samplingRounds == 0 and (tfBackground ||  geneBackground)){
		throw invalid_argument("number of background rounds -j (and dbSNP file -k) must be specified\n for help use -h"); 
	}
	if (gcMatching and randomSNPs != ""){ // no sampling, the random SNPs are given
		gcMatchingWarning = "-s is ignored, since the random SNPs are given with -i";
	}else if (gcMatching and samplingRounds == 0){
		gcMatchingWarning = "-s is ignored, since no background sampling is performed (-j)";
	}
	if (gcMatchingWarning != ""){
		gcMatching = false;
		cout << "WARNING: " << gcMatchingWarning << endl;
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
	"\n#\t-s match GC content in background sampling: " << io.gcMatching << (io.gcMatchingWarning != "" ? " (WARNING: " + io.gcMatchingWarning + ")" : "") << 
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
	bool REMs = (getREMs() != "");
	unordered_map<string, vector<string>> infoREMs; //chr:start-end_var1_var2 -> REM info (see readOverlappingREMs)
	int numREMFields = 0;
	char delim = '\t';
	int start = 0, end = 0;
	string chr = "", var1 = "", var2 = "", line = ""; //stores current line
	if (REMs){
		infoREMs = readOverlappingREMs(getOverlappingREMs());
		numREMFields = CountEntriesFirstLine(getREMs(),'\t') - 1; // REM position, gene name and all columns of the interaction file after start and end
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
		start = stoi(getToken(line, delim)) - SNV_FLANK;
		//cout << "start: " << start << endl;
		end = stoi(getToken(line, delim)) + SNV_FLANK;
		//cout << "end: " << end << endl;
		var1 = getToken(line, delim);
		///cout << "var1: " << var1 << endl;
		var2 = getToken(line, delim);
		id = getToken(line,delim);
		MAF = getToken(line,delim);
		//cout << "var2: " << var2 << endl;
		string key = chr + ":" + to_string(start+SNV_FLANK) + "-" + to_string(end-SNV_FLANK) + "_" + var1 + "_" + var2;
		//cout << "key: " << key << endl;
		if ((var1 != "*") and (var2 != "*") and (var1.length() == 1) and (var2.length() == 1)){ 
			counterOverlappingPeaks++;
			output << chr << '\t' <<  start << '\t' << end << '\t' << chr << ":" << to_string(start+SNV_FLANK) << "-"<< to_string(end-SNV_FLANK) << ";" << var1 << ";" << var2 << ";" << id << ";" << MAF;// << '\t' << var1<< '\t' <<  var2 << '\n';
			//store skipped entries as header of the fasta file
			if (getOverlappingFootprints() == inputFile){
				vector<string> entriesLine;
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
				if (infoREMs.count(key) > 0){
					counterOverlappingREMs++;
				}
				output << remColumns(infoREMs, key, numREMFields);
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

void InOutput::callHelp(){

	cout << "Call program with ./src/differentialBindingAffinity_multipleSNPs\noptional parameters:\n" << 
	"-o outputDir (default SNEEP_output/); must be empty or contain a former SNEEP output (info.txt), which is deleted\n" <<
	"-n number threads (default 1)\n" <<
	"-p pvalue for motif hits (default 0.5)\n"<<
	"-c pvalue differential binding (default 0.01)\n" <<
	"-b base frequency for PFMs -> PWMs ( /necessaryInputFiles/frequency.txt)\n" <<
	"-a if flag is set,  all computed differential binding affinities are stored in <outputDir>/AllDiffBindAffinity.txt\n"<<
	"-f additional footprint/open chromatin region file in bed file format\n" <<
	"-m if flag is set, the  maximal differential binding affinity per SNP is printed\n"<<
	"-t file where expression values of TFs are stored (e.g RNA-seq in a tab-separated format e.g. ensemblID\texpression-value)\n" <<
	"-d threshold TF activity (must be given if -t is given)\n"<<
	"-e tab-separated file containing ensemblID to gene name mapping of the TFs (must be given if -t is given)\n"<<
	"-r bed-like file with epigenetic interactions\n"<<
	"-g path to file containing ensemblID to gene name mapping, must be given if -r is given (,-separated)(mapping for all genes within EpiRegio)\n" <<
	"-j rounds sampled background (default 0)\n" <<
	"-k path to sorted dbSNP file (required for -j, unless -i is given)\n" <<
	"-i directory with already sampled random SNPs (randomSNPs_<round>.txt, e.g. sampling/ of a former run), used instead of sampling them (-j must be set)\n" <<
	"-l start seed (default 1)\n" <<
	"-q minimal TF count which needs to be exceeded to be considered in random sampling (default 0)\n" << 
	"-s true or false, match the GC content (+- 30bp around the SNV) in addition to the MAF in the random sampling, requires a dbSNP file with GC content (-k) (default false)\n" <<
	"-u gene background analysis is performed (default false), -j must be set \n" <<
	"-v perform TF enrichment  analysis (default  false), -j must be set\n" <<	
	"-x transition matrix for binding affinity p-value, (default all transitions are equally likely) (necessaryInputFiles/transitionMatrix.txt)\n" <<
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

				if (seen_alleles.count(alleles) == 0){
					//SNPs.push_back(helper);
					seen_alleles.insert(alleles);
					output << originalLine << '\n';
					counterUnique++;
				}
			}else{
				previous_pos = current_pos;
				seen_alleles.clear();
				seen_alleles.insert(a1 + "-" + a2); // alleles of the first line at this position (otherwise one duplicate is kept)
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
	bool REMs = (getREMs() != "");
	unordered_map<string, vector<string>> infoREMs; //chr:start-end_var1_var2 -> REM info (see readOverlappingREMs)
	int numREMFields = 0;
	char delim = '\t';
	int start = 0, end = 0;
	string chr = "", var1 = "", var2 = "", line = ""; //stores current line
	if (REMs){
		infoREMs = readOverlappingREMs(REMsOverlappFile);
		numREMFields = CountEntriesFirstLine(getREMs(),'\t') - 1; // same columns as for the input SNPs
	}
	//cout << "done withe REM file" << endl;
	ofstream output(outputFile);
	//ofstream output2(getBedFileInDels());//store InDels
	string id = "", MAF = "";
	int counter_commas = 0;
	mt19937 generator(seed); // different seed per round (each round can run in its own thread) 	
	int randomNum = 0; //sampled unifrom distributed number

	while (getline(input, line, '\n')){
		line = line + '\n';	
		//getline(input, line, '\n'); //getLine
		//extract information
		chr = getToken(line, delim);
		//cout << "chr: " << chr << endl;
		start = stoi(getToken(line, delim)) - SNV_FLANK;
		//cout << "start: " << start << endl;
		end = stoi(getToken(line, delim)) + SNV_FLANK;
		//cout << "end: " << end << endl;
		var1 = getToken(line, delim);
	//	cout << "var1: " << var1 << endl;
		var2 = getToken(line, delim);
		id = getToken(line,delim);
		MAF = getToken(line,delim);
		//cout << "var2: " << var2 << endl;
		string key = chr + ":" + to_string(start+SNV_FLANK) + "-" + to_string(end-SNV_FLANK) + "_" + var1 + "_" + var2;
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
		output << chr << '\t' <<  start << '\t' << end << '\t' << chr << ":" << to_string(start+SNV_FLANK) << "-"<< to_string(end-SNV_FLANK) << ";" << var1 << ";" << var2 << ";" << id << ";" << MAF << ";.";// << '\t' << var1<< '\t' <<  var2 << '\n';
		// add REM info
		if (REMs){
			output << remColumns(infoREMs, key, numREMFields);
		}
		output << '\n';
	}
	input.close();
	output.close();
	return;
}




/*
* reads the REMs overlapping SNPs (bedtools intersect -wa -wb of the interaction file and a SNP file)
* returns per SNP (key chr:start-end_var1_var2) the REM info: REM position (chr:start-end), ensemblID, gene name and all
* further columns of the interaction file; if a SNP overlaps several REMs, the entries are comma separated
*/
unordered_map<string, vector<string>> InOutput::readOverlappingREMs(string overlapFile){

	unordered_map<string, vector<string>> infoREMs;
	char delim = '\t';
	string line = "", ensembl = "";
	int entries = CountEntriesFirstLine(getREMs(),'\t'); // columns of the interaction file
	//read mapping ensembl id -> gene Name
	unordered_map<string, string> mappingGeneNames;
	ifstream mapping(getMappingGeneNames());
	while (getline(mapping, line, '\n')){
		ensembl = getToken(line, ',');
		mappingGeneNames[ensembl] = line;
	} 
	ifstream overlappingREMs(overlapFile);
	while (getline(overlappingREMs, line, '\n')){
		line = line + '\n';
		vector<string> helper; //stores info per line
		string chr = getToken(line, delim);
		string start = getToken(line, delim);
		string end = getToken(line, delim);
		helper.push_back(chr + ":" + start + "-" + end);
		for (int i = 3; i <= entries-1; ++i){ //read only entries REMs
			string field = getToken(line, delim);
			helper.push_back(field); //ensembl id, REMid, coefficient, pvalue, normModelScore, meanDNase1Signal, stdDNase1Signal, consortium, version
			if (i == 3){//determine gene name
				helper.push_back(mappingGeneNames[field]);
			}
		}
		//read SNP Info
		chr =  getToken(line,delim);
		start =  getToken(line,delim);
		end =  getToken(line,delim);
		string var1 =  getToken(line,delim);
		string var2 =  getToken(line,delim);
		string key = chr + ":" + start + "-" + end + "_" + var1 + "_" + var2;

		if (infoREMs.count(key)>0){ //key exists: add new info to existing entries
			vector<string>& existingInfo = infoREMs[key];
			for(size_t i = 0; i < existingInfo.size(); ++i){
				existingInfo[i] = existingInfo[i] + "," + helper[i]; 
			}
		}else{
			infoREMs[key] = helper;
		}
	}
	overlappingREMs.close();
	return infoREMs;
}

/*
* REM columns of a SNP for the bed file (";" + field each), ";." per field if the SNP overlaps no REM
*/
string InOutput::remColumns(unordered_map<string, vector<string>>& infoREMs, const string& key, int numREMFields){
	string result = "";
	auto found = infoREMs.find(key);
	if (found == infoREMs.end()){
		for (int i = 0; i < numREMFields; ++i){
			result += ";.";
		}
	}else{
		for(auto& elem : found->second){
			result += ";" + elem;
		}
	}
	return result;
}

void InOutput::readScaleValues(string scaleFile, unordered_map<string, double>& scales){

	string line = "";
	string motif = "";
	double scale = 0.0;
	string skip = "";
	string originalLine = "";
	int lineNumber = 1;
	ifstream input(scaleFile); //open scaleFile
	if (!input){
		throw invalid_argument("cannot open scale file: " + scaleFile);
	}
	getline(input, line, '\n'); // skip header
	while (getline(input, line, '\n')){
		lineNumber++;
		originalLine = line;
		line += '\n'; // getToken needs a delimiter after the last token
		try{
			motif = getToken(line, '\t');
			skip = getToken(line, '\t'); // scale newton
			skip = getToken(line, '\t'); // MSE
			scale = stod(getToken(line, '\t'));
		}catch (const exception& e){
			throw invalid_argument("invalid scale file " + scaleFile + ", line " + to_string(lineNumber) + " (expected tab-separated: motif, estimatedScale, MSE, optimizedScale, ...): " + originalLine);
		}
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
int InOutput::getSeed(){
	return this->seed;
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
bool InOutput::getGCMatching(){
	return this->gcMatching;
}
#endif/*HANDLEINOUTPUT_HPP*/
