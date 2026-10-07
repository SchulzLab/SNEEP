#ifndef CALLBASHCOMMAND_HPP
#define CALLBASHCOMMAND_HPP

#include <cstdlib> // for system call
#include <string>
#include <iostream> 
#include <ostream>
#include <stdexcept>
#include <sys/wait.h> // exit status of system()

using namespace std;

class BashCommand{

	public:
	//constructor
	BashCommand(); //TODO set the path to bedtools in the constructor
	BashCommand(string genome); //TODO set the path to bedtools in the constructor
	~BashCommand(); //destructor

	void intersect(string file_SNPs, string regions, string ouput, string options);
	void getFasta(string bed_file, string output, string options);
	void mkdir(string dir, string options, bool remove);
	void rm(string dir);
	//void callPythonScriptCheckActiveMotifs(string sourceDir, string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_id, double threshold, string outputDir);
	void callPythonScriptCheckActiveMotifs(string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_id, double threshold, string outputDir);
	//void callPythonScriptSplitPFMs(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir);
	void callPythonScriptSplitPFMs(string TransfacPFMs, string PFMsDir, string outputDir);
	//void callPythonScriptSplitPFMsSELEX(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir);
//	void callPythonScriptSplitSEMs(string sourceDir, string TransfacPFMs, string PFMsDir);
	void anyCommand(string command);
	void run(string command); // runs a shell command, stops if it fails
	void sort(string options, string input, string output);
	void uniq(string options, string input, string output);
	//void callHistogram(string input, string output, string sourceDir);
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

//destructor
BashCommand::~BashCommand()
{
}

//options usually -wa -wb
void BashCommand::intersect(string file_SNPs, string regions, string output, string options){

	string command = "bedtools intersect " + options + " -a " + file_SNPs + " -b " + regions + " > " + output;
	//cout << command << endl;
	run(command);
	return;
}

//options usually -name
void BashCommand::getFasta(string bed_file, string output, string options){

	string command = "bedtools getfasta " + options + " -fi " + genome_ + " -bed " + bed_file + " -fo " + output;
//	cout << command << endl;
	run(command);
	return;
}
//options usually -p (no error if existing and creating all parent dirs if necessary)
void BashCommand::mkdir(string dir, string options, bool remove){
	string command = "mkdir " + options + " " + dir; 
	run(command);
	if (remove == true){
		string command = "rm -r -f " + dir + "/*";
		run(command);
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

	run(command);
	return;
}

//void BashCommand::callPythonScriptSplitPFMs(string sourceDir, string TransfacPFMs, string PFMsDir, string outputDir){
void BashCommand::callPythonScriptSplitPFMs( string TransfacPFMs, string PFMsDir, string outputDir){
	//string command = "python3 ./" + sourceDir + "/src/seperatePFMs.py " + TransfacPFMs + " " + PFMsDir;
	//string command = "python3 " + sourceDir + "/src/seperatePFMs.py " + TransfacPFMs + " " + PFMsDir + " " + outputDir + "/motifInfo.txt" ;
	string command = "seperatePFMs.py " + TransfacPFMs + " " + PFMsDir + " " + outputDir + "/motifInfo.txt" ;
	run(command);
	return;
}

/*
void BashCommand::callPythonScriptSplitSEMs(string sourceDir, string TransfacPFMs, string PFMsDir){
	string command = "python3 ./" + sourceDir + "/seperateSEMs.py " + TransfacPFMs + " " + PFMsDir;
	cout << command << endl;
	run(command);
	return;
}
*/
void BashCommand::rm(string dir){
	string command = "rm -f " + dir + "/*";
	run(command);
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
	run(command);
	return;
}

void BashCommand::anyCommand(string command){
	run(command);
	//cout << command<< endl;
	return;
}

/*
* runs a shell command and stops SNEEP if it fails (exit status != 0),
* otherwise e.g. a missing tool (bedtools, python3) leads to empty files and silently wrong results
*/
void BashCommand::run(string command){
	int status = system(command.c_str());
	if (status != 0){
		int exitCode = WIFEXITED(status) ? WEXITSTATUS(status) : status;
		cerr << "ERROR: command failed (exit code " << exitCode << "): " << command << endl;
		throw runtime_error("command failed: " + command);
	}
	return;
}

void BashCommand::sort(string options, string input, string output){
	// for chr positions options should be -k1,1 -k2,2n 
//	string command = "sort " + options + " -o " + input + "{,}";

	string command = "sort " + options  + " " +  input  + "> " + output;
	cout << command << endl;
	run(command);	
	return;
}
void BashCommand::uniq(string options, string input, string output){
	string command = "uniq " + options + " " + input + " > " + output;
	cout << command << endl;
	run(command);	
	return;
}

#endif/*CALLBASHCOMMAND_HPP*/
