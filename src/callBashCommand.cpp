#include "callBashCommand.hpp"


int main(){
	cout << "testtest" << endl;

	//create object
	BashCommand bc1;
	
	//call intersect
	string snp_file  = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/project_Luxembourg/test_snp.txt";
	string region_file =  "/MMCI/MS/EpiregDeep/work/TFtoMotifs/project_Luxembourg/test_1.txt";
	string output = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/project_Luxembourg/hihi2.fa";
	string option = "-wa -wb";
	
	bc1.intersect(snp_file, region_file, output, option);
	cout << "called intersection" << endl;

	//call getFasta
	output = "/MMCI/MS/EpiregDeep/work/TFtoMotifs/project_Luxembourg/hihi.fa";
	option = "-name";
	bc1.getFasta(region_file, output, option);
	cout << "called getFasta" <<endl;	
	return 0;
}





