#ifndef CALLBASHCOMMAND_HPP
#define CALLBASHCOMMAND_HPP

#include <cstdlib> // for system call
#include <string>
#include <iostream> 
#include <ostream>

using namespace std;

//defaul setting for bedtools path and genome
//string BEDTOOLS = "/MMCI/MS/DEEP-liver/work/Tools/bedtools/bedtools2/bin/";
//string BEDTOOLS = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/exportBedtools.sh";

//string GENOME = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/hg38.fa";
//string GENOME = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/allelSpecificHeart_Grote/Mus_musculus.GRCm38.dna.primary_assembly.fa";

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
	void callPythonScriptCheckActiveMotifs(string sourceDir, string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_id, double threshold);
	void callPythonScriptSplitPFMs(string sourceDir, string TransfacPFMs, string PFMsDir);
	void callPythonScriptSplitSEMs(string sourceDir, string TransfacPFMs, string PFMsDir);

	private:
	string path_bedtools;
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
//	cout << command << endl;
	system(command.c_str());
}

//options usually -name
void BashCommand::getFasta(string bed_file, string output, string options){

	string command = "bedtools getfasta " + options + " -fi " + genome_ + " -bed " + bed_file + " -fo " + output;
//	cout << command << endl;
	system(command.c_str());
}
//options usually -p (no error if existing and creating all parent dirs if necessary)
void BashCommand::mkdir(string dir, string options, bool remove){
	string command = "mkdir " + options + " " + dir; 
	system(command.c_str());
	if (remove == true){
		string command = "rm -r " + dir + "/*";
		system(command.c_str());
	}
}	

void BashCommand::callPythonScriptCheckActiveMotifs(string sourceDir, string activeTFs, string TransfacPFMs, string PFMsDir, string ensemble_name, double threshold){
	cout << "sourceDir: " << sourceDir<< " activeTFs: " << activeTFs << " PFMs: " << TransfacPFMs << "PFMdir: " << PFMsDir << " ensembl names: " << ensemble_name << " threshold: " << threshold << endl; 

	string command = "python3 ./" + sourceDir + "/seperatePFMsAndCheckActivity.py " + activeTFs + " "+ TransfacPFMs + " " +  PFMsDir + " " + ensemble_name + " " + to_string(threshold);
//	cout << "command: " << command << endl;

	system(command.c_str());
}

void BashCommand::callPythonScriptSplitPFMs(string sourceDir, string TransfacPFMs, string PFMsDir){
	string command = "python3 ./" + sourceDir + "/seperatePFMs.py " + TransfacPFMs + " " + PFMsDir;
	system(command.c_str());
}

void BashCommand::callPythonScriptSplitSEMs(string sourceDir, string TransfacPFMs, string PFMsDir){
	string command = "python3 ./" + sourceDir + "/seperateSEMs.py " + TransfacPFMs + " " + PFMsDir;
	cout << command << endl;
	system(command.c_str());
}

void BashCommand::rm(string dir){
	string command = "rm " + dir + "/*";
	system(command.c_str());
}
#endif/*CALLBASHCOMMAND_HPP*/
