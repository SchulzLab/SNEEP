#ifndef CALLBASHCOMMAND_HPP
#define CALLBASHCOMMAND_HPP

#include <cstdlib> // for system call
#include <string>
#include <iostream> 
#include <ostream>

using namespace std;

class BashCommand{

	public:
	//constructor
	BashCommand(); //TODO im constructor path to bedtools setzen
	BashCommand(string genome); //TODO im constructor path to bedtools setzen
	~BashCommand(); //deconstructor

	void intersect(string file_SNPs, string regions, string ouput, string options);
	void getFasta(string bed_file, string output, string options);
	void mkdir(string dir, string options, bool remove);
	void rm(string dir);
	//void callPythonScriptCheckActiveMotifs(string sourceDir, string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_id, double threshold, string outputDir);
	void callPythonScriptCheckActiveMotifs(string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_id, double threshold, string outputDir);
	//void callPythonScriptSplitPFMs(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir);
	void callPythonScriptSplitPFMs(string TransfacPFMs, string PFMsDir, string outputDir);
	//void callPythonScriptSplitPFMsSELEX(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir);
	void callPythonScriptSplitPFMsSELEX(string TransfacPFMs, string PFMsDir, string outputDir);
//	void callPythonScriptSplitSEMs(string sourceDir, string TransfacPFMs, string PFMsDir);
	void bedtoolsRandom(int len, int num, int seed, string genomes, string output);
	void cut(string options, string inputFile, string outputFile);
	void bedtoolsShuffle(string randomSequences, string genomesFile, string excludedSeq, int seed, string options, string output);
	void anyCommand(string command);
	void sed(string options, string input, string output);
	void sort(string options, string input, string output);
	void uniq(string options, string input, string output);
	//void callHistogram(string input, string output, string sourceDir);
	string getGenomeFile();
	void callFormatingScript(string file, string formatedSNPFile);

	private:
//	string path_bedtools;
	string genome_;
};


//------template definition-------
//constructor
BashCommand::BashCommand()
//:genome_(GENOME)
{	
}

BashCommand::BashCommand(string genome)
:genome_(genome)
{
}

//deconstructor
BashCommand::~BashCommand()
{
}

//options usually -wa -wb
void BashCommand::intersect(string file_SNPs, string regions, string output, string options){

	string command = "bedtools intersect " + options + " -a " + file_SNPs + " -b " + regions + " > " + output;
	//cout << command << endl;
	system(command.c_str());
	return;
}

//options usually -name
void BashCommand::getFasta(string bed_file, string output, string options){

	string command = "bedtools getfasta " + options + " -fi " + genome_ + " -bed " + bed_file + " -fo " + output;
//	cout << command << endl;
	system(command.c_str());
	return;
}
//options usually -p (no error if existing and creating all parent dirs if necessary)
void BashCommand::mkdir(string dir, string options, bool remove){
	string command = "mkdir " + options + " " + dir; 
	system(command.c_str());
	if (remove == true){
		string command = "rm -r -f " + dir + "/*";
		system(command.c_str());
	}
	return;
}	

//void BashCommand::callPythonScriptCheckActiveMotifs(string sourceDir, string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_name, double threshold, string outputDir){
void BashCommand::callPythonScriptCheckActiveMotifs(string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_name, double threshold, string outputDir){
//	cout << "sourceDir: " << sourceDir<< "/src/ activeTFs: " << activeTFs << " PFMs: " << TransfacPFMs << "PFMdir: " << PFMsDir << " ensembl names: " << ensemble_name << " threshold: " << threshold << endl; 

	//string command = "python3 ./" + sourceDir + "/src/seperatePFMsAndCheckActivity.py " + activeTFs + " "+ TransfacPFMs + " " +  PFMsDir + " " + ensemble_name + " " + to_string(threshold);
	//string command = "python3 "  + sourceDir + "/src/seperatePFMsAndCheckActivity.py " + activeTFs + " "+ TransfacPFMs + " " +  PFMsDir + " " + ensemble_name + " " + to_string(threshold) + " " +  outputDir + "/motifInfo.txt";
	string command = "seperatePFMsAndCheckActivity.py " + activeTFs + " "+ TransfacPFMs + " " +  PFMsDir + " " + ensemble_name + " " + to_string(threshold) + " " +  outputDir + "/motifInfo.txt";
//	cout << "command: " << command << endl;

	system(command.c_str());
	return;
}

//void BashCommand::callPythonScriptSplitPFMs(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir){
void BashCommand::callPythonScriptSplitPFMs( string TransfacPFMs, string PFMsDir, string outputDir){
	//string command = "python3 ./" + sourceDir + "/src/seperatePFMs.py " + TransfacPFMs + " " + PFMsDir;
	//string command = "python3 " + sourceDir + "/src/seperatePFMs.py " + TransfacPFMs + " " + PFMsDir + " " + outputDir + "/motifInfo.txt" ;
	string command = "seperatePFMs.py " + TransfacPFMs + " " + PFMsDir + " " + outputDir + "/motifInfo.txt" ;
	system(command.c_str());
	return;
}

//only for snp selex data from f1000 paper https://f1000research.com/articles/11-33#ref6
//void BashCommand::callPythonScriptSplitPFMsSELEX(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir){
void BashCommand::callPythonScriptSplitPFMsSELEX(string TransfacPFMs, string PFMsDir, string outputDir){
	//string command = "python3 ./" + sourceDir + "/src/seperatePFMs.py " + TransfacPFMs + " " + PFMsDir;
	//string command = "python3 " + sourceDir + "/src/seperatePFMs_SNP_SELEX.py " + TransfacPFMs + " " + PFMsDir + " " + outputDir + "/motifInfo.txt" ;
	string command = "seperatePFMs_SNP_SELEX.py " + TransfacPFMs + " " + PFMsDir + " " + outputDir + "/motifInfo.txt" ;
	cout << command << endl;
	system(command.c_str());
	return;
}
/*
void BashCommand::callPythonScriptSplitSEMs(string sourceDir, string TransfacPFMs, string PFMsDir){
	string command = "python3 ./" + sourceDir + "/seperateSEMs.py " + TransfacPFMs + " " + PFMsDir;
	cout << command << endl;
	system(command.c_str());
	return;
}
*/
void BashCommand::rm(string dir){
	string command = "rm -f " + dir + "/*";
	system(command.c_str());
	return;
}

void BashCommand::bedtoolsRandom(int len, int num, int seed, string genomes, string output){

	string command = "bedtools random -l " + to_string(len) + " -n " + to_string(num) + " -seed " + to_string(seed) + " -g " + genomes + ">" +  output;
	cout << command << endl;
	system(command.c_str());
	return;
}

void BashCommand::cut(string options, string inputFile, string outputFile){

	string command = "cut " + options + " " + inputFile + ">"  + outputFile;
	//cout << command << endl;
	system(command.c_str());
	return;
}

void BashCommand::bedtoolsShuffle(string randomSequences, string genomesFile, string excludedSeq, int seed, string options, string output){

	string command = "";	
	if (excludedSeq != "no"){
		if (options != "no"){
			command = "bedtools shuffle -excl " + excludedSeq  + " -i " + randomSequences + " -g " +  genomesFile + " " + options + " > " + output;
		} else {
			command = "bedtools shuffle -excl " + excludedSeq  + " -i " + randomSequences + " -g " +  genomesFile +  " > " + output; 
		}
	} else { 
		if (options != "no"){
			command = "bedtools shuffle  -i " + randomSequences + " -g " +  genomesFile + " " + options + " > " + output;
		} else {
			command = "bedtools shuffle  -i " + randomSequences + " -g " +  genomesFile + " > " + output;
			cout << command << endl;
		}
	}
	//cout << command << endl;
	system(command.c_str());
	return;
}

//void BashCommand::callHistogram(string input, string output, string sourceDir){
//
//	string command = "Rscript " + sourceDir + "/src/histogram.R " + input + " " + output;
//	
//	system(command.c_str());
//
//}

void BashCommand::callFormatingScript(string file,string  formatedSNPFile){

	string command = "formatVCF.py " + file + " " + formatedSNPFile;
	//cout << command << endl;
	system(command.c_str());
	return;
}

void BashCommand::anyCommand(string command){
	system(command.c_str());
	//cout << command<< endl;
	return;
}

void BashCommand::sed(string options, string input, string output){

	//system(command.c_str());
	string command = "sed " + options + " " + input + ">" + output;
	return;
}

void BashCommand::sort(string options, string input, string output){
	// for chr positions options should be -k1,1 -k2,2n 
//	string command = "sort " + options + " -o " + input + "{,}";

	string command = "sort " + options  + " " +  input  + "> " + output;
	cout << command << endl;
	system(command.c_str());	
	return;
}
void BashCommand::uniq(string options, string input, string output){
	string command = "uniq " + options + " " + input + " > " + output;
	cout << command << endl;
	system(command.c_str());	
	return;
}

string BashCommand::getGenomeFile(){
	return this->genome_;
}

#endif/*CALLBASHCOMMAND_HPP*/
