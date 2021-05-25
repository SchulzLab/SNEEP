#ifndef SAMPLERANDOMRSIDS2_HPP
#define SAMPLERANDOMRSIDS2_HPP

#include <string>
#include <stdexcept>
#include <iostream> 
#include <algorithm>
#include <fstream>
#include <sstream>
#include <random>
#include <array>
#include <unordered_map>
#include <vector>

//for parallelization 
#include <omp.h>


//own classes
#include "callBashCommand.hpp"

int MAXIMAL_ROUNDS = 1000;

using namespace std;

class rsIDsampler{

	public:
	//Constructor
	rsIDsampler(string pathTodbSNPFile_, vector<double>& MAF_, BashCommand& bc_);
	rsIDsampler(double binwidth_, string pathTodbSNPFile_, vector<double>& MAF_, BashCommand& bc_);
	//Deconstructor
	~rsIDsampler();

	//functions

	unordered_map<double, int> splitMAFinBins();
	string getToken(string& line, char delim);
	vector<string> determineRandomSNPs(unordered_map<double,int>& MAF_counter, string outputDir, int rounds, int numThreads, string sourceDir, int seed);
	void sampleSNPs(int counter, string outputFile, vector<string>& SNPs, int seed, int size, string histogramFile);
//	void determineRandomSamples(unordered_map<double,int>& MAF_counter, string outputFile, unordered_map<double, vector<string>>& dbSNPs, int seed, string histogram, vector<double>& keys);
//	unordered_map<double, vector<string>> storeDbSNPs();

	//getter
	double getBinwidth();
	string getdbSNPFile();
	vector<double> getMAF();

	private: // all variables with a getter should be private
	double binwidth;
	BashCommand bc; //bash command object	
	string dbSNPFile = "";
	vector<double> MAF;
};
//constructor default
rsIDsampler::rsIDsampler(string pathTodbSNPFile_, vector<double>& MAF_, BashCommand& bc_)
:binwidth(0.01), dbSNPFile(pathTodbSNPFile_),MAF(MAF_), bc(bc_)
{
//	BashCommand bc;
	cout << "constructor with default setting (num_reg = 1 and len = 100)" << endl;
}

// constructor normally used
rsIDsampler::rsIDsampler(double binwidth_, string pathTodbSNPFile_, vector<double>& MAF_, BashCommand& bc_)
:binwidth(binwidth_), dbSNPFile(pathTodbSNPFile_),MAF(MAF_), bc(bc_)
{
	cout << "constructor with user defind settings" << endl;
}

//deconstructor
rsIDsampler::~rsIDsampler()
{
}

/*
// reads SNPs bin per bin from the dbSNP file and determines random SNPs for the given rounds per bin
*/
vector<string> rsIDsampler::determineRandomSNPs(unordered_map<double,int>& MAF_counter, string outputDir, int rounds, int numThreads, string sourceDir, int seed){

	double bin = -1.0, currentMAF = 0.0; 
	int counter = 0,  size = 0; //alternative 53330
	vector<string> currentSNPs;
	string line = "";
	
	//create vector that holds rounds as string (not as int as in the for loop)
	vector<string> SNP_files (rounds, "");
	vector<string> histogram_files (rounds, "");
	string num = "";
	ofstream h; 
	for(int r = 0; r < rounds; ++r){
		num = to_string(r);
		SNP_files[r] = outputDir + "/randomSNPs_" + num + ".txt";
		histogram_files[r] = outputDir + "/histogram_" + num + ".txt";
		//write header of the histogram file
		h.open(histogram_files[r]);
		h << "SNP\tMAF\n"; 
		h.close();
	}
	
	ifstream inputFile(dbSNPFile); //open dbSNPFile
	while (getline(inputFile, line, '\n')){
		//cout << "line: " << line << endl;
		currentMAF = stod(getToken(line, '\t'));
		if (currentMAF <= bin){
			currentSNPs.push_back(line);
		}else{
			//determine MAF counter
			if (MAF_counter.count(bin) >0){ //checks if the bin exist in MAF_counter
				counter = MAF_counter[bin];
				size = currentSNPs.size();
				//cout << "bin " << bin << endl;
				for(int r = 0; r < rounds; ++r){
					sampleSNPs(counter, SNP_files[r] , currentSNPs, seed + r, size, histogram_files[r]);
				}
			}	
			//increase bin, clear currentSNPs and increase seed
			currentSNPs.clear();
			seed+=MAXIMAL_ROUNDS; //increase seed per bin
			if (bin == -1){ //set to 0 if bin was -1
				bin = 0.0;
			}else{
				bin += binwidth;
			}
			//add new SNP to currentSNPs
			currentSNPs.push_back(line);
		}
	}
	//sample random SNPs for the last bin
	if (MAF_counter.count(bin) >0){ //checks if the bin exist in MAF_counter
		//cout << "bin " << bin << endl;
		counter = MAF_counter[bin];
		size = currentSNPs.size();
		for(int r = 0; r < rounds; ++r){
			sampleSNPs(counter, SNP_files[r], currentSNPs, seed + r, size, histogram_files[r]);
		}
	}
	//plot histograms for random SNPs
	#pragma omp parallel for num_threads(numThreads)
	for(int r = 0; r < rounds; ++r){
		bc.callHistogram(histogram_files[r], outputDir +  "/histrogram_" + to_string(r) + ".pdf", sourceDir);// plot for control
	}
		
	return SNP_files;
}
/*
// determine random SNPs from current bin
*/
void rsIDsampler::sampleSNPs(int counter, string outputFile, vector<string>& SNPs, int seed, int size, string histogramFile){

	mt19937 generator(seed); // seed muss fuer jeden thread ein andere sein 	
	uniform_int_distribution<int> distribution(0, size -1); //specifiy distribution of the random number
	int randomNum = 0; //sampled unifrom distributed number
	string currentSNPs, helper, SNP;
	currentSNPs.reserve(counter * 80); //80 is the expected length of a snp string 
	helper.reserve(counter*10); //10 is the expected length of a MAF

	vector<int> randomNumbers;// stores random number we already considered

	for (int j = 0; j < counter; j++){

		//sample random number
		randomNum = distribution(generator); // generat random number
//		cout << "random number: " << randomNum << endl;
		while (find(randomNumbers.begin(), randomNumbers.end(), randomNum) != randomNumbers.end()){ //randomNum already seen
			randomNum = distribution(generator); // generat random number
		}
		randomNumbers.push_back(randomNum); //add randomNum to already used ones
		SNP = SNPs[randomNum] + '\n';
		currentSNPs.append(SNP);
		helper.append("randomSNPs\t" + SNP.substr(SNP.rfind('\t')+1));

	}
	//write sampled SNPs to file and to histogram file
	ofstream output;
	output.open(outputFile, std::ofstream::app);	
	output << currentSNPs;
	output.close();
	ofstream h; 
	h.open(histogramFile, std::ofstream::app);
	h <<  helper; 
	h.close();
}

/*
/ splits the MAFs from the ordiginal dataset in bins and count how many entries we observed per bin
/ return unordered map key: bin and value the count of how many MAFs lie within this bin
*/
unordered_map<double,int> rsIDsampler::splitMAFinBins(){

	double bin = -1;
	int counter = 0;
	unordered_map<double,int> MAF_counter;
	for (auto& i : MAF){
		if (i <= bin){
			counter++;
		}else{
			if (counter != 0){//we don't want to store bins with 0 entries
				MAF_counter[bin] = counter; //store number of MAFs for the current bin
			}
			counter = 0; //reset counter
			if (bin == -1){ //set to 0 if bin was -1
				bin = 0.0;
			}else{
				bin += binwidth;
			}
			while (i > bin){ //check if the next bin is also emty
			//	cout << "while" << endl;
				bin += binwidth;
			}
			counter++; //if i <= bin increase counter
		}
	}
	//add last bin 
	if (counter != 0){
		MAF_counter[bin] = counter;
	}
	return MAF_counter;
}

/*
/ Input: line and delim as char
/ Ouput: return the substring until the delim and it cuts this token from line
*/
string rsIDsampler::getToken(string& line, char delim){
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


//getter
double rsIDsampler::getBinwidth(){
	return this->binwidth;
}
string rsIDsampler::getdbSNPFile(){
	return this->dbSNPFile;
}
vector<double> rsIDsampler::getMAF(){
	return this->MAF;
}


/*
/ reads the (prepared and sorted) dbSNP file into an unordered_map 
/ keys are the maximal bin (maximal MAF value) and values are the srtings (chr11:56297697-56297698_rs1000000008_T_A_-1)
/ takes around 4 min for bin 0.01 and 4 min for bin 0.001 
*/
/*unordered_map<double, vector<string>> rsIDsampler::storeDbSNPs(){

	unordered_map<double, vector<string>> dbSNPs; //store all entries from dbSNPFile
	double bin = -1.0, currentMAF = 0.0;
	vector<string> currentSNPs;
	string line = "";
	int emptyBins = 0;

	ifstream inputFile(dbSNPFile); //open dbSNPFile
	while (getline(inputFile, line, '\n')){
		//cout << "line: " << line << endl;
		currentMAF = stod(getToken(line, '\t'));
		getToken(line, '\t'); //skip rsID
		if (currentMAF <= bin){
			currentSNPs.push_back(line);
		}else{
			if (currentSNPs.size() == 0){ //check if there are empty bins
				emptyBins++;
				cout << "Es gibt ein leeres Bin: " << bin << endl;
			}
			//store currentSNPs with according bin value
			dbSNPs[bin] = currentSNPs;	
			//increase bin and clear currentSNPs
			currentSNPs.clear();
			if (bin == -1){ //set to 0 if bin was -1
				bin = 0.0;
			}else{
				bin += binwidth;
			}
			while (currentMAF > bin){ //check if the next bin is also emty
				dbSNPs[bin] = currentSNPs;	
				bin += binwidth;
				emptyBins++;
				cout << "Es gibt ein leeres Bin: " << bin << endl;
			}
			//add new SNP to currentSNPs
			currentSNPs.push_back(line);
		}
	}
	//add last bin to dbSNPs
	dbSNPs[bin] = currentSNPs;	
	cout << "number empty bins: " << emptyBins << endl;	
	inputFile.close(); //close dbSNP file

	return dbSNPs;
}*/

/*
/
/ samples random SNPs with a similar MAF distribution than the orginal data set 
/ outputFile is a bed file with the sampled SNPs
/
*/
/*void rsIDsampler::determineRandomSamples(unordered_map<double,int>& MAF_counter, string outputFile, unordered_map<double, vector<string>>& dbSNPs, int seed, string histogram, vector<double>& keys){

	cout << "seed: " << seed << endl;
	ofstream output;//open output file
	output.open(outputFile);	
	ofstream h; 
	h.open(histogram);
	h << "SNP\tMAF\n"; 

	double bin = 0.0;
	int counter = 0, size = 0;
	vector<string> SNPs;
	string currentSNP = "";
	vector<int> randomNumbers;// stores random number we already considered
	mt19937 generator(seed); // seed muss fuer jeden thread ein andere sein 	
	int randomNum = 0; //sampled unifrom distributed number
//	vector<string> sampledSNPs; 
//	vector<string> helper;

	#pragma omp parallel for  private(randomNumbers, bin, counter, SNPs, size, randomNum, currentSNP) num_threads(5)
	for (int i = 3; i < MAF_counter.size(); i++){
		bin = keys[i];
		counter = MAF_counter[bin];
		//extract SNPs per bin
		SNPs = dbSNPs[bin];
		//sample random SNPs from vector
		size = SNPs.size();
		if (size == 0){
			cout << "size is 0" << endl;
			//TODO: decide what to do
		}
		randomNumbers.clear();
		uniform_int_distribution<int> distribution(0, size -1); //specifiy distribution of the random number
		for (int j = 0; j < counter; j++){

			//sample random number
			randomNum = distribution(generator); // generat random number
//			cout << "random number: " << randomNum << endl;
			while (find(randomNumbers.begin(), randomNumbers.end(), randomNum) != randomNumbers.end()){ //randomNum already seen
				randomNum = distribution(generator); // generat random number
//				cout << "try again " << randomNum << endl;
			}
//			cout << "picked number: " << randomNum << endl;
			randomNumbers.push_back(randomNum); //add randomNum to already used ones
//			cout << "randomSNP: " << SNPs[randomNum] << endl;;
			currentSNP = SNPs[randomNum];
			//cout << "currentSNP: " << currentSNP << endl;
			//helper.push_back(currentSNP);
			string helper = currentSNP.substr(currentSNP.rfind('\t')+1);
			#pragma omp critical
			{
				output << currentSNP << "\n";
				h << "randomSNP\t" << helper << '\n';
			}

		}
//		#pragma omp critical
//		{
//			sampledSNPs.insert(sampledSNPs.end(), helper.begin(), helper.end());
//		}
	}
	
//	for (auto& elem : sampledSNPs){
//		output << elem << "\n";
//		h << "randomSNP\t" << elem.substr(elem.rfind('\t')+1) << '\n';
//	}
	h.close();
	output.close();
}
*/

/*
// initialized the set of sequences/
// requieres an output file and a path to a genome file
void rsIDsampler::initializeSeqSet(string output_file){	
	//create gemoes file for bedtools random using the fai-genome file assumes that this file usually exists 
	string genome = bc.getGenomeFile();
	genomes_file = genome.substr(0, genome.find(".")) + ".genomes";
	bc.sed("/_/d", genome + ".fai", genomes_file); //remove entries we are not intressed in 
	number_regions_extended = num_reg + (num_reg* 0.2); // increase the number of regions by 20 % (TODO: is this enough?)
	bc.bedtoolsRandom(len, number_regions_extended, seed, genomes_file, output_file);
//	size = number_regions_extended * LETTERS; // for bed file: ~ 35 chars per line
//	size2 = (LETTERS_FASTA*num_reg) + (num_reg*len ) ; // for FASTA file 50 chars for the header and len many for the fasta seq

//	char buffer[size] = {0};
//	b.startBuffering(buffer, size);

//	BedRandom *r= new BedRandom(genomes_file, number_regions_extended, seed, true, len);
//	delete r;
//	b.stopBuffering();
	//cout << string(buffer) << endl;
//	b.writeBufferToFile(output_file, buffer);	
	return;
}

void rsIDsampler::shuffleSequences(string randomSequences, string outputShuffle){
//string rsIDsampler::shuffleSequences(string randomSequences, string outputShuffle){

	//TODO:  bei jedem aufruf den seed um 5 erhuoehen sonst kommen immer die gleichen sequenzen raus
	string genome = bc.getGenomeFile();
	genomes_file = genome.substr(0, genome.find(".")) + ".genomes";
	//bedtools consol call
//	bedtools direct call
//	auto t1 = std::chrono::high_resolution_clock::now();
	bc.bedtoolsShuffle(randomSequences, genomes_file, coding_regions, seed,"-noOverlapping", outputShuffle);
//	char buffer[size] = {0};
//	string empty = "";
//	b.startBuffering(buffer, size);
//	BedShuffle *shuffle = new BedShuffle(randomSequences, genomes_file, coding_regions, empty, true, true, false, false, 0.0, 123, true, false, 1000, true, false);
//	delete shuffle;
//	b.stopBuffering();	

//	return string(buffer);
	return;
}

//function allow to store the curl output 
size_t CurlWrite_CallbackFunc_StdString(void *contents, size_t size, size_t nmemb, std::string *s)
{
    size_t newLength = size*nmemb;
    try
    {
        s->append((char*)contents, newLength);
    }
    catch(std::bad_alloc &e)
    {
        //handle memory problem
        return 0;
    }
    return newLength;
}

string rsIDsampler::sampleSNPs(string sequenceFile, int seedDis, string name, string outputDir, vector<double>& MAF){
//string rsIDsampler::sampleSNPs(int seedDis, string name, string outputDir, string buffer){

	string url = "", region = "", seq = ""; //stores lines from file
	int numExtractSNPs = 0, numSNPs = 0; // number of SNPs overlapping with the sampled region
	mt19937 generator(seedDis); // seed muss fuer jeden thread ein andere sein 	
	int randomNum = 0; //sampled unifrom distributed number
	nlohmann::json SNP; // randomly selected SNP
	vector<int> sampledNumbers;
	string chr = "", SNP_id = "", wildtype = "", mutant = "";

	ofstream bedFile;//open output file
	string bedFilePath = outputDir + name + ".bed";
	bedFile.open(bedFilePath);
	
	ifstream sequences(sequenceFile); //open shuffled sequence  file
	while (getline(sequences, seq, '\n')){ // for each line extract SNP
//	int n = count(buffer.begin(), buffer.end(), '\n'); // count number of lines in buffer
//	cout << "n's: " <<  n << endl;
//	for(int i = 0; i < n; i++){
		if (numSNPs == num_reg){ //already identified enough SNPs
			cout << "identified number of requiered SNPs" << endl;
			break;
		}	
//		seq = getToken(buffer, '\n');
//		cout << "--------------" << endl;
//		cout << seq << endl;
		nlohmann::json jsonObj;
		region = getRegion(seq);
		url = "https://rest.ensembl.org/overlap/region/human/" + region + "?feature=variation";
		//cout << "url: " << url << endl;

		jsonObj = REST_API_overlap(url); //REST_API call for the given region
//		cout << jsonObj << endl;
		numExtractSNPs = jsonObj.size();

//NOTE: das dauert zu lange (und überlappt vielleicht mit coding regions), vielleicht besser ein paar zusätzliche regionen zu samplen und aufzuhören sobald man genüegnd snps gefunden hat
/*		while (numExtractSNPs == 0){
			//extend region	
			region = extendRegion(region, 25); // need to test if this is enough (extends region in total by 50 base pairs)
			url = "https://rest.ensembl.org/overlap/region/human/" + region + "?feature=variation";
			cout << "extended url: " << url << endl;
			jsonObj = REST_API_overlap(url);
			numExtractSNPs = jsonObj.size();
		}
*/
/*
//		cout << "numExtracted SNPs: " << numExtractSNPs << endl;
		if (numExtractSNPs > 0){
			//if snp longer than one base pair throw another one -> there it can also be that all of them are more than i base pair long 
			uniform_int_distribution<int> distribution(0, numExtractSNPs -1);

			randomNum = distribution(generator); // generat random number
//			cout << randomNum << endl;
			nlohmann::json SNP = jsonObj[randomNum];
//			cout << "SNP: " << SNP << endl;
			nlohmann::json allels =  jsonObj[randomNum]["alleles"];
			wildtype = allels[0];
			mutant = allels[1];
//			cout << allels << endl;
//			cout << "hier: " << wildtype.size() << endl;
			sampledNumbers = {};
			bool SNP_found = true;
			//check if all found snps are not valid
			
			while  ((wildtype.size() > 1) || (mutant.size() > 1) || (allels.size() > 2) ||  (wildtype == "-") || (mutant == "-")){ //checks if a valid SNP was selected
				sampledNumbers.push_back(randomNum);
				if (sampledNumbers.size() < (numExtractSNPs/2)){
					randomNum = distribution(generator);
//					cout << "sample again a random number " << randomNum << endl;
					if (find(sampledNumbers.begin(), sampledNumbers.end(),randomNum) == sampledNumbers.end()){ //new SNP has not been selected before
						nlohmann::json helper = jsonObj[randomNum];
//						cout << "SNP: " << helper << endl;
						allels =  helper["alleles"];
						wildtype = allels[0];
						mutant = allels[1];
						SNP = helper;
					}else{
						cout << "number already sampled " <<  randomNum << endl;
					}
				}else{
					cout << "too many tries, we stop here" << endl; 
					SNP_found = false;
					break;
				}				
			}
			if (SNP_found == true){
				numSNPs++;
	//			cout << SNP << endl;
				int start = SNP["start"];
				start -=1; // dbSNPs are 1-based
				int end = SNP["end"];
				chr= SNP["seq_region_name"];
				SNP_id = SNP["id"];
				cout << "SNP_id: " << SNP_id << endl;
				determineMAF(SNP_id, MAF);
				
				bedFile << "chr" << chr << '\t' <<  start - 50 << '\t' << end + 50<< '\t' << "chr" << chr << ':' <<  start << '-' << end << ";" << wildtype << ";" << mutant << ";" << SNP_id <<  '\n';
			}
		}
	}

	bedFile.close(); //really important otherwise getFasta is not working!
	string fastaFile = outputDir + name + ".fa";

	//call fasta function 
//	char buffer2[size2] = {0};
//	b.startBuffering(buffer2, size2);
	//Bed2Fa *getFasta= new Bed2Fa(bc.getGenomeFile(), bedFilePath, fastaFile, true, false, false, false, false, false, false, true);
	//delete getFasta;
//	b.stopBuffering();

	bc.getFasta(bedFilePath, fastaFile, "-name");

	return fastaFile;	
}


/*
*
* extend regions at both ends by the given len parameter
* input: region in chr:start-end format and number
* 
*/
/*
string rsIDsampler::extendRegion(string region, int len){

	string chr = getToken(region, ':');
	int start = stoi(getToken(region, '-')) - len;
	int end = stoi(region) + len;
	
	if (start < 0){
		start = 0;
	}
	cout << "chr: " << chr << " start: " << start << " end: " << end << endl;

	return ( chr + ":" + to_string(start) + "-" + to_string(end));
}

nlohmann::json rsIDsampler::REST_API_overlap(string url){

	//initialzise curl for REST api request 
	CURL* curl;
	CURLcode res;
	curl_global_init(CURL_GLOBAL_DEFAULT);
	struct curl_slist *hs=NULL;
	hs = curl_slist_append(hs, "Content-Type: application/json");
	string s;

	curl = curl_easy_init();

	if (curl){

		curl_easy_setopt(curl, CURLOPT_HTTPHEADER, hs); //set header such that curl return json -H option
		//curl_easy_setopt( curl, CURLOPT_URL, "https://rest.ensembl.org/overlap/region/human/7:140424943-140424959?feature=variation" ); //rest api request 
		curl_easy_setopt( curl, CURLOPT_URL, url.c_str()); //rest api request, must be a c string, since curl is a c library 

		curl_easy_setopt(curl, CURLOPT_SSL_VERIFYPEER, 0L); //only for https
        	curl_easy_setopt(curl, CURLOPT_SSL_VERIFYHOST, 0L); //only for https

		// Perform the request, res will get the return code (0 for everthing is ok)
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, CurlWrite_CallbackFunc_StdString);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, &s); //s stores result of the request as string 
		//curl_easy_setopt (curl, CURLOPT_VERBOSE, 1L); //remove this to disable verbose output
    		res = curl_easy_perform(curl);
		//cout << res << endl;

		// Check for errors  
		if(res != CURLE_OK)
			fprintf(stderr, "curl_easy_perform() failed: %s\n", curl_easy_strerror(res));
				
		curl_easy_cleanup(curl); //cleanup
  	}
	if (s.find("Service Unavailable") != std::string::npos){
		cout <<"service unavailable case: " <<  s << endl;
		s = "[]";
	}
	// parse string to json object 
	nlohmann::json jsonObj;
	stringstream(s) >> jsonObj; //needs to parse string with request result to json for easy access, json object can be handle as a STL container
	//cout << jsonObj << endl;
	return jsonObj;
}

string rsIDsampler::getRegion(string& line){
	string chr =  getToken(line, '\t');
	chr.erase(0,3);
	string start = getToken(line, '\t');
	string end = getToken(line, '\t');

	return (chr + ":" + start + "-" + end);

}


string rsIDsampler::getToken(string& line, char delim){
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

/*
*function to determine MAFs for the input SNPs
*
* currently it identifies per region the SNP id
* TODO: auch andere optionen ermöglich -> IDs können schon gegeben sein oder MAF (wenn es keine snps sind die in dbSNP auftauchen)
*/
/*
vector<double> rsIDsampler::determineMAFsForSNPs(string bedFile, string& outputDir){

	string url = "", line = "", rsID = "", w_1 = "", m_1 = "", w = "" , m = "", chr = "", start = "", maf = "";
	int numSNPs = 0;		
	bool skip = false, match = false;
	vector<double> MAF;

	//open file where we store the ids 
	ofstream fileIds;
	fileIds.open(outputDir + "snpsWithIds.txt");
	//extract rsID via region query (ensembl)
	//read bedfile (the one without duplicates and indels)
	ifstream bedfile(bedFile); //open sequence file
	while (getline(bedfile, line, '\n')){ // for each line extract overlapping SNPs
		//cout << "hihi" << endl;
		skip = false;
		getToken(line, '\t'); //skip first 3 entries
		getToken(line, '\t'); 
		getToken(line, '\t');
		chr = getToken(line, ':');
		fileIds << chr << '\t' << getToken(line, '-') << '\t'; //skip start
		start = getToken(line, ';'); //1-based in dbSNP
	
		url = "https://rest.ensembl.org/overlap/region/human/" + chr + ":" + start + ":" + start  + "?feature=variation";
		w = getToken(line, ';'); //wildtype base 
		m = getToken(line, ';'); //mutated base
		rsID = getToken(line, ';');
		maf = getToken(line, ';');
		if (maf != "."){ //MAF is alreadt determined -> notung to do
			fileIds << start << '\t' << w << '\t' << m << '\t' << rsID << '\t' << maf << '\n';

		}else if (rsID != "."){ //rsID is known -> just determine MAF
			fileIds << start << '\t' << w << '\t' << m << '\t' << rsID << '\t';
			determineMAF(rsID, MAF, fileIds);

		}else{ //determine rsID and MAF (if possible)
			fileIds << start << '\t' << w << '\t' << m << '\t';
		//	cout << "wildtype: " << w << " mutant: " << m << endl;
		//	cout << url << endl;	
			nlohmann::json jsonObj = REST_API_overlap(url);
		//	cout << "json object: " << jsonObj << endl;
			numSNPs = jsonObj.size();
			if (numSNPs == 0){ //TODO: stop if json obj is empty -> sieht das so aus?
				fileIds << "-\t-\n";
				continue;
			}else{
				nlohmann::json SNP;
				nlohmann::json helper;
				//check if the allels are the same as in the input file
				//TODO: check if it is an indel
				if (numSNPs > 1){//did we identify more than one entry for the position?
					for (int j = 0; j < numSNPs; j++){
						skip = false;
						helper =  jsonObj[j]["alleles"];
						//cout << "helper: " << helper << endl;
						for (int i = 0; i < helper.size(); i++){ //check if the identified SNP is an indel
							if (string(helper[i]).size() > 1 or helper[i] == "-"){ //if so mark snp as indel
								//cout << "longer than 1: " << helper[i] << endl;
								skip = true;
							}
						}
						if (skip == false){ //check if we identified no indel
							SNP = jsonObj[j]; 
							break;
						}
					}
				}else{
					SNP = jsonObj[0]; 
				}
				if (skip == true){ //means no SNP was selected
					cout << jsonObj << endl;
					throw invalid_argument ("no SNP selected!!");
				}	
				nlohmann::json allels = SNP["alleles"];
				//needs to consider that there are more than 2 allels. which is the correct one?
				match = false; //did one of the allele combis fit?
				w_1 = allels[0];
				//cout << "wiltype 1: " << w_1 << endl;
				for (int i = 1; i< allels.size(); i++){
					m_1 = allels[i];
					if ((w == w_1 && m == m_1) || (w == m_1 && m == w_1)){ //this is the combination from the input file
						//determine MAF
						rsID = SNP["id"];
						fileIds << rsID << "\t";
						determineMAF(rsID, MAF, fileIds);
						match = true;
						break;
					}
				}
				if (match == false){
					fileIds << "-\t-\n" << endl;
					cout << SNP << endl;
					throw invalid_argument ("allels are not fitting: " + w  + " " + m);
				}			
			}
		}
	}
	fileIds.close();
	return MAF;
}

void rsIDsampler::determineMAF(string rsID, vector<double> &MAF, ofstream& fileIds){

	string url = "https://rest.ensembl.org/variation/human/" + rsID + "?content-type=application/json";
//	cout << url << endl;
	nlohmann::json jsonObj = REST_API_overlap(url);
//	cout << jsonObj << endl;
	double maf = jsonObj["MAF"];
	MAF.push_back(maf);
	fileIds << maf << "\n";
	return;	
}

void rsIDsampler::determineMAF(string rsID, vector<double> &MAF){

	string url = "https://rest.ensembl.org/variation/human/" + rsID + "?content-type=application/json";
//	cout << url << endl;
	nlohmann::json jsonObj = REST_API_overlap(url);
	cout << jsonObj << endl;
	cout << jsonObj["MAF"] << endl;
	try{
		MAF.push_back(jsonObj["MAF"]);
	}catch (...){
		cout << "An exception occurred. Exception Nr. " << endl;
	}

	//TODO: else fall festlegen
	return;
}

//getter
int rsIDsampler::getNumRegions(){
	return this->num_reg;
}
int rsIDsampler::getLength(){
	return this->len;
}
int rsIDsampler::getSeedValue(){
	return this->seed;
}
string rsIDsampler::getCodingRegions(){
	return this->coding_regions;
}*/

#endif/*SAMPLERANDOMRSIDS2_HPP*/
