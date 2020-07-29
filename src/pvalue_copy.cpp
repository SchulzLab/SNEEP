#include "pvalue_copy.hpp"
#include "Matrix_new.hpp"

string PATH_TO_TRANSITION_MATRIX = "../transition_matrix.txt";

int main(){

	ifstream transition_matrix_file(PATH_TO_TRANSITION_MATRIX); //open file from transition matrix
	Matrix<double> transition_matrix; 
	transition_matrix_file >> transition_matrix; //create transition matrix
	cout  << "transition matrix:\n" << transition_matrix << endl;
	cout << "hfkdshfs" << endl;
	pvalue p(3, 1);
	cout << "acc: " << p.getAcc() << endl;	
	cout << "thres: " << p.getThres() << endl;	
	cout << "Kmer size: " << p.getKmerSize() << endl;	
	cout << "---------------------------" << endl;

	cout << "----------------TEST1-------------------" << endl;

	// PWM	
	Matrix<double> PWM(4,3);
	PWM(1,1) = 4;
	PWM(1,2) = 1;
	PWM(1,3) = 2;
	PWM(2,1) = 3;
	PWM(2,2) = 2;
	PWM(2,3) = 2;
	PWM(3,1) = 1;
	PWM(3,2) = 4;
	PWM(3,3) = 3;
	PWM(4,1) = 2;
	PWM(4,2) = 1;
	PWM(4,3) = 2;
	cout << "PWM\n";
	cout <<PWM << endl;

	vector<double> freq(4, 0.25);
	for( auto& i: freq)
		cout << i << " ";
	cout << endl;
	
	vector<double> vec = p.calculatePvalues(PWM, freq, transition_matrix);
	cout << "thres: " << p.getThres() << endl;
	cout << "vec: " << endl;	
	for (auto& i: vec)
		cout << i << " ";
	cout << endl;
	cout << "-------------TEST2----------------" << endl;

	
	pvalue p2(3, 0.1);
	cout << "acc: " << p2.getAcc() << endl;	
	cout << "thres: " << p2.getThres() << endl;	
	cout << "Kmer size: " << p2.getKmerSize() << endl;	
	// PWM	
	Matrix<double> PWM2(4,3);
	PWM2(1,1) = 0.1;
	PWM2(1,2) = 0.9;
	PWM2(1,3) = 0.2;
	PWM2(2,1) = 0.2;
	PWM2(2,2) = 0.1;
	PWM2(2,3) = 0.1;
	PWM2(3,1) = 0.8;
	PWM2(3,2) = 0.1;
	PWM2(3,3) = 0.8;
	PWM2(4,1) = 0.4;
	PWM2(4,2) = 0.7;
	PWM2(4,3) = 0.9;
	cout <<PWM2 << endl;

	vector<double> freq2(4, 0.0);
	freq2[0]= 0.2;
	freq2[1]= 0.5;
	freq2[2]= 0.2;
	freq2[3]= 0.1;
	for( auto& i: freq2)
		cout << i << " ";
	cout << endl;

	vector<double> vec2 = p2.calculatePvalues(PWM2, freq2, transition_matrix);
	cout << "thres: " << p2.getThres() << endl;	
	cout << "vec2: " << endl;	
	for (auto& i: vec2)
		cout << i << " ";
	cout << endl;
	cout << "-------------------test3------------" << endl;

	pvalue p3(6,0.001);
	Matrix<double> PWM3;
	
	ifstream input("../random_sequence_examples/PWM_ARNT/ARNT.txt");	
	input>> PWM3;
	cout << PWM3 << endl;
	vector<double> vec3 = p3.calculatePvalues(PWM3,freq, transition_matrix);
	cout << "thres: " << p3.getThres() << endl;	
	cout << "vec3: " << endl;	
	for (auto& i: vec3)
		cout << i << " ";
	cout << endl;

	return 0; 
	cout << "test ATF3" << endl;

	pvalue p7(8,0.001);
	Matrix<double> PWM7;
	
	ifstream input4("../../ENCODE_data/PWMs_GRCh38_Gene_ID_0.001/ATF3.txt");	
	input4>> PWM7;
	cout << PWM7 << endl;
	vector<double> vec7 = p7.calculatePvalues(PWM7,freq, transition_matrix);
	cout << "thres: " << p7.getThres() << endl;	
	cout << "vec7: " << endl;	
	for (auto& i: vec7)
		cout << i << ", ";
	cout << endl;

	cout << "test REST" << endl;

	pvalue p6(21,0.001);
	Matrix<double> PWM6;
	
	ifstream input2("../../ENCODE_data/PWMs_GRCh38_Gene_ID_0.001/REST.txt");	
	input2>> PWM6;
	cout << PWM6 << endl;
	vector<double> vec6 = p6.calculatePvalues(PWM6,freq, transition_matrix);
	cout << "thres: " << p6.getThres() << endl;	
	cout << "vec6: " << endl;	
	for (auto& i: vec6)
		cout << i << ", ";
	cout << endl;



}
