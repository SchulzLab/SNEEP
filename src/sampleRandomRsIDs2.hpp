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
//#include <omp.h>


//own classes
#include "callBashCommand.hpp"

int MAXIMAL_ROUNDS = 1500;

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
	//vector<string> determineRandomSNPs(unordered_map<double,int>& MAF_counter, string outputDir, int rounds, int numThreads, string sourceDir, int seed);
	vector<string> determineRandomSNPs(unordered_map<double,int>& MAF_counter, string outputDir, int rounds, int numThreads, int seed);
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
//	cout << "constructor with default setting (num_reg = 1 and len = 100)" << endl;
}

// constructor normally used
rsIDsampler::rsIDsampler(double binwidth_, string pathTodbSNPFile_, vector<double>& MAF_, BashCommand& bc_)
:binwidth(binwidth_), dbSNPFile(pathTodbSNPFile_),MAF(MAF_), bc(bc_)
{
//	cout << "constructor with user defind settings" << endl;
}

//deconstructor
rsIDsampler::~rsIDsampler()
{
}

/*
// reads SNPs bin per bin from the dbSNP file and determines random SNPs for the given rounds per bin
*/
//vector<string> rsIDsampler::determineRandomSNPs(unordered_map<double,int>& MAF_counter, string outputDir, int rounds, int numThreads, string sourceDir, int seed){
vector<string> rsIDsampler::determineRandomSNPs(unordered_map<double,int>& MAF_counter, string outputDir, int rounds, int numThreads, int seed){

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
		//histogram_files[r] = outputDir + "/histogram_" + num + ".txt";
		//write header of the histogram file
		//h.open(histogram_files[r]);
		//h << "SNP\tMAF\n"; 
		//h.close();
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
					//cout << r << endl;
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
	//	cout << "bin " << bin << endl;
		counter = MAF_counter[bin];
		size = currentSNPs.size();
		for(int r = 0; r < rounds; ++r){
			sampleSNPs(counter, SNP_files[r], currentSNPs, seed + r, size, histogram_files[r]);
		}
	}
	//plot histograms for random SNPs
	/*#pragma omp parallel for num_threads(numThreads)
	for(int r = 0; r < rounds; ++r){
		bc.callHistogram(histogram_files[r], outputDir +  "/histrogram_" + to_string(r) + ".pdf", sourceDir);// plot for control
	}*/
		
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
//	h.open(histogramFile, std::ofstream::app);
//	h <<  helper; 
//	h.close();
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

#endif/*SAMPLERANDOMRSIDS2_HPP*/
