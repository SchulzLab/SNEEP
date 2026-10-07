#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>

//#include <curl/curl.h> //curl library (it is  c-library)
//#include <nlohmann/json.hpp> //json library allows to handle json object as an LTD container :) 
#include <algorithm>
#include <stdio.h>
#include <unordered_map>
#include <memory>

#include "callBashCommand.hpp"

using namespace std;
//using json = nlohmann::json;

// number of flanking bases on each side of the SNV used for the GC content (window = 2 * GC_FLANK + 1 bp)
//  differentialBindingAffinity_multipleSNPs uses also 30 flanking bases up and downstream (flag -s)
const int GC_FLANK = 30;

//functions
string getToken(string& line, char delim);
string checkLength(string line);
void computeGCContent(string dbSNPBed, string genome, string gcDir, int numThreads, BashCommand& bc);
string nextGCContent(unordered_map<string, unique_ptr<ifstream>>& gcFiles, string gcDir, string chr, long lineNumber);


/*
* Input: - dbSNP file (GCF file downloaded from https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/) ## stand 01.09.2023 release 156, the old version was release dbSNP 154 -> achtung gibts in hg38 und hg19 (see readme)
* 	## stand 26.02.2026 downloaded newest release dbSNP 157 -> without excluding SNVs overlapping with coding regions
*        - output dir
*	 - genome file (fasta, hg38), the index <genome>.fai must exist (samtools faidx <genome>)
*	 - number of threads (optional, default 1), used to compute the GC content of the chromosomes in parallel
*
* Output: dbSNPs_sorted.txt with columns MAF chr start end ref alt rsID MAF GC
* 	GC = (#C + #G) / (#A + #C + #G + #T) in a window of +- GC_FLANK bp around the SNV (lower case bases are counted, N is excluded),
* 	-1 if the window contains no A, C, G or T
*
*/

int main(int argc, char *argv[]){

	if (argc<4)  // there should be 3 non-option arguments
		throw invalid_argument("Usage: getSNPInfo <dbSNP file> <outputDir> <genome.fa> [numThreads]");
	string input = argv[1];
	string outputDir = argv[2];
	string genome = argv[3];
	int numThreads = 1;
	if (argc > 4){
		numThreads = stoi(argv[4]);
	}
	cout << "dbSNPFile: " << input << " outputDir: " << outputDir << " genome: " << genome << " threads: " << numThreads << endl;
	if (!ifstream(genome + ".fai")){
		throw invalid_argument("missing genome index " + genome + ".fai (create it with: samtools faidx " + genome + ")");
	}
	//create bedfile from dbSNP file
	BashCommand bc;
	bc.mkdir(outputDir, "-p", false);
	
	
	//extract header lines staring with #, split chr such that it suits our chr notation and brings file in bed format
	string command = "awk '!/(^#|^NT|^NW|NC_012920.1)/ {split($1,a,\".\"); split(a[1], b, \"_\"); if (substr(a[1], length(a[1])-1,1) != 0) print \"chr\" substr(a[1], length(a[1])-1,1) substr(a[1], length(a[1]),1) \"\t\" $2-1 \"\t\" $2 \"\t\" $3 \"_\" $4 \"_\" $5 \"_\" $8 ; else print \"chr\" substr(a[1], length(a[1]),1) \"\t\" $2-1 \"\t\" $2 \"\t\" $3 \"_\" $4 \"_\" $5 \"_\" $8; }' "  +  input + " > " + outputDir + "/dbSNP.bed" ;
	cout << command << endl;
	bc.anyCommand(command);

	// inplace replacement of chr23 and chr24 to chrX and chrY
	command = "sed -i -e 's/chr23/chrX/g' " + outputDir + "/dbSNP.bed";
	cout << command << endl;
	bc.anyCommand(command);
	
	command = "sed -i -e 's/chr24/chrY/g' " + outputDir + "/dbSNP.bed";
	cout << command << endl;
	bc.anyCommand(command);

	// compute GC content around each SNV (one file per chromosome: line number in dbSNP.bed and GC content)
	string gcDir = outputDir + "/gc";
	computeGCContent(outputDir + "/dbSNP.bed", genome, gcDir, numThreads, bc);
	unordered_map<string, unique_ptr<ifstream>> gcFiles; // open GC file per chromosome
	long lineNumber = 0; // line number in dbSNP.bed (1-based as NR in awk)
	string GC = "";

	// read input file line by  line
	string line = "", rsID = "", ref = "", alt = "", url = "", chr = "", start = "", end = "", info = "", token = "";
	//string sizeRef = "", sizeAlt = "";
	//json jsonObj;
	int pos = 0;
	double MAF;

	ofstream helperFile;
	helperFile.open(outputDir + "/dbSNPs_indels.txt");

	ofstream outputFile;
	outputFile.open(outputDir + "/dbSNPs.txt");
	//outputFile.open(outputDir + "/dbSNPs_GWAS_Catalog.txt");
	ifstream SNP_file(outputDir + "/dbSNP.bed"); //open input file
	while (getline(SNP_file, line, '\n')){
		line += '\n';
		lineNumber++;
		chr = getToken(line, '\t');
		GC = nextGCContent(gcFiles, gcDir, chr, lineNumber); // read also for indels to stay in sync with the GC file
		start = getToken(line, '\t'); 
		end = getToken(line, '\t');
		rsID = getToken(line, '_'); 
		ref = getToken(line, '_');
		//cout << "ref: " << ref << endl;
		ref = checkLength(ref); //check if it is a snp or an indel
		//cout << "ref: " << ref << endl;
		alt = getToken(line, '_');
		info = getToken(line, '\n');
		//cout << "alt: " << alt << endl;
		alt = checkLength(alt); //check if it is a snp or an indel
		//cout << "alt: " << alt << endl;

		//cout << chr << " "  << start << " " << end << " "  << rsID << " " << ref << " " <<  alt << endl;
		if (ref == "NO" or alt == "NO"){
			helperFile << rsID << '\n';
			continue;
		}else{
			//cout << "ref: " << ref << " alt: " << alt << " id: " << rsID <<  endl; 
			pos = info.find("FREQ");
			if (pos != std::string::npos){
				//cout << info << endl;
    				token = info.substr(pos, info.size());
				pos = token.find(":");
    				token = token.substr(pos+1, token.size());
				pos = token.find("|"); //check if there are more than one MAF -> considere first one
				if (pos != std::string::npos){
    					token = token.substr(0, pos);
				}				
				pos = token.find(";"); //check if it is a common SNP - if so remove info
				if (pos != std::string::npos){
    					token = token.substr(0, pos);
				}				
				//getToken(token, ','); //skip allelFreq of wildtype 
				vector<double> allMAFs;
				string helper = "";
				while(token.find(",") != std::string::npos){
					helper= getToken(token, ',');
					if (helper != "."){
						allMAFs.push_back(stod(helper));	
					}

				}
				if (token != "."){
					allMAFs.push_back(stod(token));
				}
				
				int d  = distance(allMAFs.begin(), min_element(allMAFs.begin(), allMAFs.end()));
				MAF = allMAFs[d];
				//cout << "MAF: " <<  MAF << endl;
			}else{ // FREQ not in info -> MAF is null
			//	cout << "NO FREQ" << endl;
				MAF = -1;	
			
			}
			// write output 
			outputFile << MAF  <<  "\t" << chr << "\t"  << start << "\t" << end << "\t"  << ref << "\t" <<  alt << "\t" << rsID << "\t" << MAF << "\t" << GC << endl;
		}
	}
	SNP_file.close();
	outputFile.close();
	gcFiles.clear(); // close GC files
	command = "rm -r " + gcDir;
	cout << command << endl;
	bc.anyCommand(command);
	string tmpDir = outputDir + "/tmp"; // temporary files of sort (needed for large files)
	bc.mkdir(tmpDir, "-p", false);
	command = "sort -k1,1 -T " + tmpDir + " --parallel=" + to_string(numThreads) + " -g " + outputDir + "/dbSNPs.txt" + " > " + outputDir + "/dbSNPs_sorted.txt"; //sorts resulting SNPs (-g because of e-5 schreibweise)
	cout << command << endl;
	if (system(command.c_str()) != 0){ // check exit status, otherwise a failed sort leaves an incomplete dbSNPs_sorted.txt unnoticed
		throw runtime_error("sorting " + outputDir + "/dbSNPs.txt failed");
	}
	command = "rm -r " + tmpDir;
	cout << command << endl;
	bc.anyCommand(command);

	helperFile.close();
	return 0;
}

/*
* Input: dbSNP bed file, genome file (fasta, with .fai index), output dir for the GC files, number of threads
*
* computes the GC content in a window of +- GC_FLANK bp around each SNV, separately per chromosome:
* bedtools slop (window, clipped at the chromosome ends) | bedtools nuc (base counts)
* writes per chromosome <gcDir>/gc_<chr>.txt with line number in the dbSNP bed file and GC content
* GC = (#C + #G) / (#A + #C + #G + #T), N is excluded, -1 if there is no A, C, G or T
*/
void computeGCContent(string dbSNPBed, string genome, string gcDir, int numThreads, BashCommand& bc){

	bc.mkdir(gcDir, "-p", true);
	// split SNV positions per chromosome, 4th column is the line number in the dbSNP bed file
	string command = "awk -v dir=" + gcDir + " 'BEGIN{OFS=\"\\t\"} {if (!($1 in seen)) {seen[$1] = 1; print $1 > (dir \"/chromosomes.txt\")} print $1, $2, $3, NR > (dir \"/snvs_\" $1 \".bed\")}' " + dbSNPBed;
	cout << command << endl;
	bc.anyCommand(command);

	vector<string> chromosomes;
	string chr = "";
	ifstream chrFile(gcDir + "/chromosomes.txt");
	while (getline(chrFile, chr, '\n')){
		chromosomes.push_back(chr);
	}
	chrFile.close();

	// bedtools nuc columns (4 user columns): 7 num_A, 8 num_C, 9 num_G, 10 num_T
	string gcCommand = "awk 'BEGIN{OFS=\"\\t\"} NR > 1 {acgt = $7 + $8 + $9 + $10; if (acgt > 0) print $4, ($8 + $9) / acgt; else print $4, -1}'";
	cout << "bedtools slop -i " << gcDir << "/snvs_<chr>.bed -g " << genome << ".fai -b " << GC_FLANK << " | bedtools nuc -fi " << genome << " -bed stdin | " << gcCommand << " > " << gcDir << "/gc_<chr>.txt" << endl;
	#pragma omp parallel for schedule(dynamic) num_threads(numThreads)
	for (int i = 0; i < (int)chromosomes.size(); ++i){
		string snvs = gcDir + "/snvs_" + chromosomes[i] + ".bed";
		string currentCommand = "bedtools slop -i " + snvs + " -g " + genome + ".fai -b " + to_string(GC_FLANK) + " | bedtools nuc -fi " + genome + " -bed stdin | " + gcCommand + " > " + gcDir + "/gc_" + chromosomes[i] + ".txt && rm " + snvs;
		bc.anyCommand(currentCommand);
	}
	return;
}

/*
* Input: open GC files per chromosome, dir of the GC files, chromosome and line number of the current SNV in the dbSNP bed file
*
* returns the GC content of the current SNV (next line of the GC file of the chromosome)
* throws an error if the GC file does not belong to the current line (e.g. a window skipped by bedtools)
*/
string nextGCContent(unordered_map<string, unique_ptr<ifstream>>& gcFiles, string gcDir, string chr, long lineNumber){

	string file = gcDir + "/gc_" + chr + ".txt";
	if (gcFiles.count(chr) == 0){
		gcFiles[chr] = unique_ptr<ifstream>(new ifstream(file));
	}
	string line = "";
	if (!getline(*gcFiles[chr], line, '\n')){
		throw invalid_argument("missing GC content for line " + to_string(lineNumber) + " of the dbSNP bed file in " + file);
	}
	line += '\n';
	if (stol(getToken(line, '\t')) != lineNumber){
		throw invalid_argument("GC content in " + file + " does not match line " + to_string(lineNumber) + " of the dbSNP bed file");
	}
	return getToken(line, '\n');
}


/*
* Input: line for istance from a file, and a delim like '/t'
* 
* returns the new element until the delim symbol and cuts the element from the inpit line
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

/*
* Input: ref or alt string
*
* checks if the ref or the alt snp is a indel (longer than 1)
*
*/

string checkLength(string line){
	string result = ""; 
	int counter = count(line.begin(), line.end(), ','); 
	string helper = "";
	bool added = false;
	if (counter > 0){
		for (int i = 0; i< counter; i++){
			helper = getToken(line, ',');
			//if (helper.size() > 1){
			//	return (false);
			//}
			//if (helper == "N" || helper == "n"){
			//	cout << "es gibt SNPs mit N" << endl;
			//	return (false);
			//}
			if (helper == "A" || helper == "C" || helper == "G" || helper == "T"){
				added = true;
				result += helper + ",";
				//cout << "was bist du? " << helper << endl;  
				//return false;
			}
		}

		//check last token
		if (line == "A" || line == "C" || line == "G" || line == "T"){
			added = true;	
			result += line;
		//	cout << "was bist du last token? " << line << endl;  
		//	return false;
		}else{
			result.pop_back(); // removes the last comma which is not needed
		}

		if (added == false){
			result = "NO";
		}

	}else{
		if (line.size() >1){
			result = "NO";
		//	return(false);
		}
		//TODO: muss man auch für - prüfen?
		//if (line == "-"){
		//	cout << "ja muss man" << endl;
		//}
		if (line == "A" || line == "C" || line == "G" || line == "T"){
			result = line;
			//cout << "was bist du? " << line << endl;  
			//return false;
		}else{
			result = "NO";
		}
	}
	return (result);
}
