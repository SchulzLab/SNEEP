#include "Matrix_new.hpp"

//column -> row
// notice that the matrix index starts by 0 but number of rows and number of columns starts by one
// zum testen 
int main(){

	cout << "HIHI" << endl;

	//create matrix with number of rows and number of columns
	Matrix<int> m_(2,2);
	cout << "........." << endl;
	cout << m_ << endl; 
	//read element on position  1,1 (index matrix 0,0)
	cout << m_(1,1)	<< endl;

	//write element
	m_(1,2) = 2;
	cout << m_(1,2)<< endl;
	cout << "......" << endl;
	// check << opertaor
	cout << m_ << endl;
	//check  >> operator
	cout << "check <<" << endl;
	Matrix<float> hihi;
	ifstream input("../small_example/PWMs/GFI1B_PWM.txt");	
	//ifstream input("test.txt");	
	input>> hihi;
	cout << hihi << endl;
	cout << hihi(1,1)<< endl;
	cout << hihi(1,2)<< endl;
	cout << hihi(1,3)<< endl;
	cout << hihi(2,1)<< endl;
	cout << hihi(2,2)<< endl;
	cout << hihi(2,3)<< endl;
	cout << hihi(3,1)<< endl;
	cout << hihi(3,2)<< endl;
	cout << hihi(3,3)<< endl;
	cout << hihi(4,1)<< endl;
	cout << hihi(4,2)<< endl;
	cout << hihi(4,3)<< endl; 

	//check copy constructur
	Matrix<float> c; 
	cout << "c: \n" << c << endl;
	c = hihi; 
	cout << "-----------------------__" << endl;
	cout << "c: \n" << c << endl;
	return 0;

}
