#ifndef SAMPLERANDOMRSIDS2_HPP
#define SAMPLERANDOMRSIDS2_HPP

#include <string>
#include <stdexcept>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <array>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <vector>

//for parallelization
//#include <omp.h>


//own classes
#include "callBashCommand.hpp"
#include "stringUtils.hpp"

const int MAXIMAL_ROUNDS = 1500; // seed offset per bin, larger than the number of rounds

// number of flanking bases on each side of the SNV used for the GC content (window = 2 * GC_FLANK + 1 bp)
// must be the same as in getSNPInfo.cpp, which precomputes the GC content of the dbSNP SNVs
const int GC_FLANK = 30;

using namespace std;

/*
* GC content in a window of +- flank bp around position center of seq
* GC = (#C + #G) / (#A + #C + #G + #T), lower case bases are counted, N (and other symbols) are excluded
* returns -1 if the window contains no A, C, G or T, numN and windowLength are set for the warnings
* the result is rounded to 6 significant digits, as written by getSNPInfo (awk), such that equal values end up in the same bin
*/
double gcContent(const string& seq, int center, int flank, int& numN, int& windowLength){

	int start = max(0, center - flank), end = min((int)seq.size() - 1, center + flank);
	int gc = 0, acgt = 0;
	numN = 0;
	windowLength = end - start + 1;
	for (int i = start; i <= end; ++i){
		switch (seq[i]){
		case 'C': case 'c': case 'G': case 'g':
			gc++;
			acgt++;
			break;
		case 'A': case 'a': case 'T': case 't':
			acgt++;
			break;
		case 'N': case 'n':
			numN++;
			break;
		}
	}
	if (acgt == 0){
		return -1;
	}
	ostringstream rounded;
	rounded << setprecision(6) << (double)gc / acgt;
	return stod(rounded.str());
}

class rsIDsampler{

	public:
	//Constructor
	rsIDsampler(double binwidth_, string pathTodbSNPFile_, vector<double>& MAF_); // MAF matching only
	rsIDsampler(double binwidth_, double gcBinwidth_, string pathTodbSNPFile_, vector<double>& MAF_, vector<double>& GC_); // MAF x GC matching, GC_[i] = -1: SNV i is matched by MAF only
	//Destructor
	~rsIDsampler();

	//functions
	vector<string> determineRandomSNPs(string outputDir, int rounds, int seed);

	//getter
	double getBinwidth();
	string getdbSNPFile();
	vector<double> getMAF();
	vector<string> getWarnings();

	private: // all variables with a getter should be private

	// number of SNVs to sample per MAF bin: MAF only and per GC bin
	struct Request{
		int mafOnly = 0;
		map<int, int> gc; // GC bin -> count
		int total(){
			int sum = mafOnly;
			for (auto& elem : gc){
				sum += elem.second;
			}
			return sum;
		}
	};

	int binIndex(double value, vector<double>& edges, double width);
	int gcBinIndex(double gc);
	string binName(vector<double>& edges, int bin);
	void splitInBins();
	void checkGCColumn();
	void sampleBin(int seedBin, Request& request, vector<string>& lines, vector<int>& lineGCBins, int rounds, int seed, vector<string>& output, vector<long>& counts);
	void sampleSNPs(int counter, vector<string>& lines, vector<int>* pool, mt19937& generator, string& output);

	double binwidth;
	double gcBinwidth = 0.05;
	bool gcMatching = false;
	string dbSNPFile = "";
	vector<double> MAF;
	vector<double> GC;
	vector<double> mafEdges; // upper bin edges: (edge[k-1], edge[k]], edge[0] = -1 (MAF not given), edge[1] = 0
	vector<double> gcEdges; // same for GC, edge[0] = -1 (no A, C, G, T in the window)
	int numGCBins = 0; // number of GC bins incl. bin 0 (GC = -1)
	map<int, Request> requests; // MAF bin -> SNVs to sample
	vector<string> warnings;
};

// constructor MAF matching only
rsIDsampler::rsIDsampler(double binwidth_, string pathTodbSNPFile_, vector<double>& MAF_)
:binwidth(binwidth_), dbSNPFile(pathTodbSNPFile_),MAF(MAF_)
{
	splitInBins();
}

// constructor MAF x GC matching
rsIDsampler::rsIDsampler(double binwidth_, double gcBinwidth_, string pathTodbSNPFile_, vector<double>& MAF_, vector<double>& GC_)
:binwidth(binwidth_), gcBinwidth(gcBinwidth_), gcMatching(true), dbSNPFile(pathTodbSNPFile_),MAF(MAF_), GC(GC_)
{
	if (GC.size() != MAF.size()){
		throw invalid_argument("number of MAF and GC values of the input SNVs differ");
	}
	checkGCColumn();
	splitInBins();
}

//destructor
rsIDsampler::~rsIDsampler()
{
}

/*
* returns the bin index of value: smallest k with value <= edges[k]
* edges are extended as in the original implementation (-1, 0, then repeatedly += width),
* such that the bins (incl. floating point drift) are identical to the former MAF bins
*/
int rsIDsampler::binIndex(double value, vector<double>& edges, double width){
	if (edges.empty()){
		edges.push_back(-1.0);
		edges.push_back(0.0);
	}
	while (value > edges.back()){
		edges.push_back(edges.back() + width);
	}
	return lower_bound(edges.begin(), edges.end(), value) - edges.begin();
}

/*
* GC bin: 0 for GC = -1, 1 for GC = 0, then (0, gcBinwidth], ... up to (1 - gcBinwidth, 1]
*/
int rsIDsampler::gcBinIndex(double gc){
	return min(binIndex(gc, gcEdges, gcBinwidth), numGCBins - 1); // GC = 1 might exceed the last edge due to floating point drift
}

/*
* name of a bin for the warnings
*/
string rsIDsampler::binName(vector<double>& edges, int bin){
	if (bin == 0){
		return "-1 (not given)";
	}
	if (bin == 1){
		return "0";
	}
	return "(" + to_string(edges[bin - 1]) + ", " + to_string(edges[bin]) + "]";
}

/*
* counts per MAF bin (and GC bin) how many SNVs need to be sampled
*/
void rsIDsampler::splitInBins(){
	numGCBins = 2 + (int)(1.0 / gcBinwidth + 0.5); // GC = -1, GC = 0, (0, gcBinwidth], ..., (1 - gcBinwidth, 1]
	for (size_t i = 0; i < MAF.size(); ++i){
		Request& request = requests[binIndex(MAF[i], mafEdges, binwidth)];
		if (gcMatching and GC[i] != -1){
			request.gc[gcBinIndex(GC[i])]++;
		}else{
			request.mafOnly++;
		}
	}
	return;
}

/*
* the dbSNP file must contain the GC content as column 9 (see getSNPInfo.cpp)
*/
void rsIDsampler::checkGCColumn(){
	ifstream inputFile(dbSNPFile);
	if (!inputFile){
		throw invalid_argument("cannot open dbSNP file: " + dbSNPFile);
	}
	string line = "";
	getline(inputFile, line, '\n');
	if (count(line.begin(), line.end(), '\t') < 8){
		throw invalid_argument("GC content matching (-s true) requires a dbSNP file with the GC content in column 9 (created with getSNPInfo): " + dbSNPFile);
	}
	return;
}

/*
* reads the dbSNP file (sorted by MAF) bin per bin and samples random SNPs for all rounds per bin
* returns the files of the random SNPs, one per round
*/
vector<string> rsIDsampler::determineRandomSNPs(string outputDir, int rounds, int seed){

	//create vector that holds rounds as string (not as int as in the for loop)
	vector<string> SNP_files (rounds, "");
	for(int r = 0; r < rounds; ++r){
		SNP_files[r] = outputDir + "/randomSNPs_" + to_string(r) + ".txt";
	}
	vector<long> counts (rounds, 0); // number of sampled SNPs per round
	vector<string> output (rounds, "");
	ofstream outputFile;

	// first pass: sample all MAF bins that exist in the dbSNP file
	set<int> existingBins;
	vector<string> currentSNPs;
	vector<int> currentGCBins;
	int currentBin = -1, bin = 0;
	string line = "";
	ifstream inputFile(dbSNPFile); //open dbSNPFile
	if (!inputFile){
		throw invalid_argument("cannot open dbSNP file: " + dbSNPFile);
	}
	while (getline(inputFile, line, '\n')){
		bin = binIndex(stod(getToken(line, '\t')), mafEdges, binwidth);
		if (bin != currentBin){
			if (bin < currentBin){
				throw invalid_argument("dbSNP file is not sorted by MAF: " + dbSNPFile);
			}
			if (!currentSNPs.empty()){
				sampleBin(currentBin, requests[currentBin], currentSNPs, currentGCBins, rounds, seed, output, counts);
				for(int r = 0; r < rounds; ++r){ // write per MAF bin to keep the memory small
					outputFile.open(SNP_files[r], std::ofstream::app);
					outputFile << output[r];
					outputFile.close();
					output[r].clear();
				}
			}
			currentSNPs.clear();
			currentGCBins.clear();
			currentBin = bin;
			existingBins.insert(bin);
		}
		if (requests.count(bin) > 0){ // only store SNPs of bins we need
			currentSNPs.push_back(line);
			if (gcMatching){
				currentGCBins.push_back(gcBinIndex(stod(line.substr(line.rfind('\t') + 1))));
			}
		}
	}
	inputFile.close();
	//sample random SNPs for the last bin
	if (!currentSNPs.empty()){
		sampleBin(currentBin, requests[currentBin], currentSNPs, currentGCBins, rounds, seed, output, counts);
	}

	// MAF bins of the input SNPs that do not exist in the dbSNP file -> use the nearest existing MAF bin
	map<int, vector<int>> fallbackBins; // existing MAF bin -> missing MAF bins
	for (auto& request : requests){
		int missing = request.first;
		if (existingBins.count(missing) > 0){
			continue;
		}
		if (existingBins.empty()){
			throw invalid_argument("dbSNP file is empty: " + dbSNPFile);
		}
		auto upper = existingBins.lower_bound(missing);
		int target = 0;
		if (upper == existingBins.end()){
			target = *existingBins.rbegin();
		}else if (upper == existingBins.begin()){
			target = *upper;
		}else{
			int above = *upper, below = *prev(upper);
			if (above - missing == missing - below){ // tie: random choice
				mt19937 generator(seed + missing * MAXIMAL_ROUNDS);
				target = uniform_int_distribution<int>(0, 1)(generator) == 0 ? below : above;
			}else{
				target = (above - missing < missing - below) ? above : below;
			}
		}
		fallbackBins[target].push_back(missing);
		warnings.push_back("MAF bin " + binName(mafEdges, missing) + ": no dbSNP SNVs, " + to_string(request.second.total()) + " SNVs sampled from the nearest MAF bin " + binName(mafEdges, target));
	}
	// second pass (only if necessary): sample the missing MAF bins from their nearest existing bin
	if (!fallbackBins.empty()){
		currentSNPs.clear();
		currentGCBins.clear();
		currentBin = -1;
		inputFile.open(dbSNPFile);
		while (true){
			bool read = (bool)getline(inputFile, line, '\n');
			if (read){
				bin = binIndex(stod(getToken(line, '\t')), mafEdges, binwidth);
			}
			if (!read or bin != currentBin){
				if (!currentSNPs.empty()){
					for (auto& missing : fallbackBins[currentBin]){ // seeds of the missing bin (never used before)
						sampleBin(missing, requests[missing], currentSNPs, currentGCBins, rounds, seed, output, counts);
					}
				}
				currentSNPs.clear();
				currentGCBins.clear();
				currentBin = bin;
			}
			if (!read){
				break;
			}
			if (fallbackBins.count(bin) > 0){
				currentSNPs.push_back(line);
				if (gcMatching){
					currentGCBins.push_back(gcBinIndex(stod(line.substr(line.rfind('\t') + 1))));
				}
			}
		}
		inputFile.close();
	}
	for(int r = 0; r < rounds; ++r){
		outputFile.open(SNP_files[r], std::ofstream::app);
		outputFile << output[r];
		outputFile.close();
	}

	// we never want to lose a SNV: each round must have as many SNVs as the input
	for(int r = 0; r < rounds; ++r){
		if (counts[r] != (long)MAF.size()){
			throw runtime_error("random SNPs round " + to_string(r) + ": " + to_string(counts[r]) + " SNVs sampled, but " + to_string(MAF.size()) + " input SNVs");
		}
	}
	return SNP_files;
}

/*
* samples the SNVs of one MAF bin for all rounds
* seedBin: MAF bin that determines the seeds (differs from the bin of the lines for the MAF fallback)
* MAF only: seed + seedBin * MAXIMAL_ROUNDS + r (as before)
* MAF x GC: seed + (seedBin * (numGCBins + 1) + slot) * MAXIMAL_ROUNDS + r, slot = GC bin or numGCBins for MAF only
*/
void rsIDsampler::sampleBin(int seedBin, Request& request, vector<string>& lines, vector<int>& lineGCBins, int rounds, int seed, vector<string>& output, vector<long>& counts){

	string mafBin = binName(mafEdges, seedBin);
	// MAF only (MAF matching or input SNVs without GC content)
	if (request.mafOnly > 0){
		int seedBase = seed + (gcMatching ? (seedBin * (numGCBins + 1) + numGCBins) : seedBin) * MAXIMAL_ROUNDS;
		if (request.mafOnly > (int)lines.size()){
			warnings.push_back("MAF bin " + mafBin + ": " + to_string(request.mafOnly) + " SNVs needed but only " + to_string(lines.size()) + " dbSNP SNVs available, SNVs are sampled more than once");
		}
		for(int r = 0; r < rounds; ++r){
			mt19937 generator(seedBase + r);
			sampleSNPs(request.mafOnly, lines, nullptr, generator, output[r]);
			counts[r] += request.mafOnly;
		}
	}
	if (request.gc.empty()){
		return;
	}
	// MAF x GC: SNVs per GC bin (dbSNP SNVs without GC content (bin 0) are never sampled)
	vector<vector<int>> pools (numGCBins);
	for (size_t i = 0; i < lineGCBins.size(); ++i){
		if (lineGCBins[i] > 0){
			pools[lineGCBins[i]].push_back(i);
		}
	}
	for (auto& elem : request.gc){
		int gcBin = elem.first, counter = elem.second;
		int seedBase = seed + (seedBin * (numGCBins + 1) + gcBin) * MAXIMAL_ROUNDS;
		string gcBinName = binName(gcEdges, gcBin);
		vector<int> candidates; // GC bins to sample from: the bin itself or the nearest non-empty bins
		if (!pools[gcBin].empty()){
			candidates.push_back(gcBin);
		}else{
			for (int d = 1; d < numGCBins and candidates.empty(); ++d){
				if (gcBin - d > 0 and !pools[gcBin - d].empty()){
					candidates.push_back(gcBin - d);
				}
				if (gcBin + d < numGCBins and !pools[gcBin + d].empty()){
					candidates.push_back(gcBin + d);
				}
			}
			if (candidates.empty()){
				warnings.push_back("MAF bin " + mafBin + ", GC bin " + gcBinName + ": no dbSNP SNV with GC content in this MAF bin, " + to_string(counter) + " SNVs sampled by MAF only");
			}else{
				warnings.push_back("MAF bin " + mafBin + ", GC bin " + gcBinName + ": no dbSNP SNVs, " + to_string(counter) + " SNVs sampled from the nearest GC bin" + (candidates.size() > 1 ? "s (random choice per round)" : ""));
			}
		}
		for (auto& c : candidates){
			if (counter > (int)pools[c].size()){
				warnings.push_back("MAF bin " + mafBin + ", GC bin " + gcBinName + ": " + to_string(counter) + " SNVs needed but only " + to_string(pools[c].size()) + " dbSNP SNVs available, SNVs are sampled more than once");
			}
		}
		for(int r = 0; r < rounds; ++r){
			mt19937 generator(seedBase + r);
			if (candidates.empty()){
				sampleSNPs(counter, lines, nullptr, generator, output[r]);
			}else if (candidates.size() == 1){
				sampleSNPs(counter, lines, &pools[candidates[0]], generator, output[r]);
			}else{
				int choice = uniform_int_distribution<int>(0, 1)(generator);
				sampleSNPs(counter, lines, &pools[candidates[choice]], generator, output[r]);
			}
			counts[r] += counter;
		}
	}
	return;
}

/*
* samples counter SNPs from lines (or from the lines given by the indices in pool) and appends them to output
* without replacement if possible (same random numbers as the former implementation),
* otherwise each SNP is taken once and the remaining ones are sampled with replacement
*/
void rsIDsampler::sampleSNPs(int counter, vector<string>& lines, vector<int>* pool, mt19937& generator, string& output){

	int size = (pool == nullptr) ? lines.size() : pool->size();
	uniform_int_distribution<int> distribution(0, size -1); //specifiy distribution of the random number
	int randomNum = 0; //sampled unifrom distributed number
	output.reserve(output.size() + counter * 80); //80 is the expected length of a snp string

	if (counter <= size){
		unordered_set<int> randomNumbers;// stores random number we already considered
		for (int j = 0; j < counter; j++){
			randomNum = distribution(generator); // generat random number
			while (randomNumbers.count(randomNum) > 0){ //randomNum already seen
				randomNum = distribution(generator); // generat random number
			}
			randomNumbers.insert(randomNum); //add randomNum to already used ones
			output.append(lines[(pool == nullptr) ? randomNum : (*pool)[randomNum]] + '\n');
		}
	}else{
		for (int j = 0; j < size; j++){ // each SNP once
			output.append(lines[(pool == nullptr) ? j : (*pool)[j]] + '\n');
		}
		for (int j = size; j < counter; j++){ // remaining ones with replacement
			randomNum = distribution(generator);
			output.append(lines[(pool == nullptr) ? randomNum : (*pool)[randomNum]] + '\n');
		}
	}
	return;
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
vector<string> rsIDsampler::getWarnings(){
	return this->warnings;
}

#endif/*SAMPLERANDOMRSIDS2_HPP*/
