#include "sampleRandomRsIDs2.hpp"
#include <stdio.h>
#include <algorithm>

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
void determineMAFsForSNPs(string bedFile, vector<double>& MAF){
	string line = "";
	double maf = 0.0;
	ifstream bedfile(bedFile); //open sequence file
	while (getline(bedfile, line, '\n')){ // for each line extract overlapping SNPs
		getToken(line, ';'); //skip chr\tstart\tend\tchr:start-end
		getToken(line, ';'); //skip wildtype 
		getToken(line, ';'); //skip mutatant
		getToken(line, ';');//skip rsID
		try{ //throws an error when MAF is smaller than double precisoins allows -> set to 0
			maf = stod(getToken(line, ';'));
		}catch (const std::out_of_range& oor){
			maf = 0.0;	
		}
		MAF.push_back(maf);
	}
	return;
}
int main(){

	BashCommand bc("/home/nbaumgarten/hg38.fa");

	string pathdbSNPs = "../dbSNP/dbSNPs_sorted.txt";
	string bedFile = "/projects/sneep/work/SNEEP/Kessler_SNP_Analysis/test_/snpRegions.bed";
	vector<double> MAF;
	cout << "determine MAF" << endl;
	determineMAFsForSNPs(bedFile, MAF); //read MAF distribution from input SNP file
	cout << "before sort" << endl;
	sort(MAF.begin(), MAF.end(), std::less<double>()); //sort MAF with default operation <
	rsIDsampler s(0.01, pathdbSNPs, MAF, bc); //initialize snp sampler
	cout << "after constructor" << endl;
	unordered_map<double,int> MAF_counter = s.splitMAFinBins(); // split original MAF distribution in bins
	cout << "after maf_counter" << endl;
	s.determineRandomSNPs(MAF_counter, "testDir", 100, 20,  1);





//	cout << "after sampling" << endl;
//	for(int i = 0; i < 20; ++i){
//		string command = "Rscript  histogram.R testDir/histogram_" + to_string(i) + ".txt testDir/histogram_" + to_string(i) + ".pdf";
//		system(command.c_str());
//	}
	

/*	vector<double> MAF; //for testing also empty
	MAF.push_back(0.012);
	MAF.push_back(0.1);
	MAF.push_back(0.02);
	MAF.push_back(0.023);
	MAF.push_back(0.03);
	MAF.push_back(0.3);
	MAF.push_back(0.32);
	MAF.push_back(0.4);
	MAF.push_back(0);
	MAF.push_back(-1);
	sort(MAF.begin(),MAF.end(), std::less<double>());
	cout << "MAF" << endl;
	for (auto& x : MAF){
		cout << x << "\t";
	}

	vector<double> MAF_2; 
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.01);
	MAF_2.push_back(0.0109);
	MAF_2.push_back(0.5);
	MAF_2.push_back(0.5);
	MAF_2.push_back(0.5);
	MAF_2.push_back(0.5);
	sort(MAF_2.begin(),MAF_2.end(), std::less<double>());

	cout << "Test sampleRandomRSIDs" << endl;
	rsIDsampler r2(pathdbSNPs, MAF);
	cout << r2.getBinwidth() << " " << r2.getdbSNPFile() << " "  << endl;
	rsIDsampler r3( 0.01,  pathdbSNPs, MAF_2);
	cout << r3.getBinwidth() << " " << r3.getdbSNPFile() << " "  << endl;

	//test store dbSNPs

	unordered_map<double, vector<string>> dbSNPs = r3.storeDbSNPs();
	
	//pritn for control
	for (auto& x : dbSNPs){
		cout << "bin: " << x.first << "\t";
		for (auto& i : x.second){
			cout << i << "\t";	
		}
		cout << "\n";
	}
	unordered_map<double, int> MAF_counter = r3.splitMAFinBins();
	cout << "MAF_counter" << endl;
	for(auto& i : MAF_counter){
		cout << i.first << " " << i.second << endl;
	}

	//call sample random snps
	cout << "AB HIER" << endl;
	vector<double> keys;
	int numSNPs = 0;
	for(auto& i : MAF_counter){
		keys.push_back(i.first);
		numSNPs += i.second;
	}
	r3.determineRandomSamples(MAF_counter, "test/randomSNPsFile.bed", dbSNPs, 22, "test/histogram.txt", keys);
*/

}


