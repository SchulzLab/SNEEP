#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>

// really really important for getline !!!!
#include <fstream>
 
using namespace std;

// necessary for operator <<
//----------------------
template <class T>
class Matrix;

template<class T>
ostream& operator<< ( ostream& os,  Matrix<T>& m_ );
template<class T>
ifstream& operator>> ( ifstream& os,  Matrix<T>& m_ );
//----------------------

template <class T>
class Matrix{
	public:

	//constructors
	Matrix();
	Matrix(int num_row, int num_col);
	~Matrix<T>(); // deconstructor
	Matrix(const Matrix<T>&); // copyconstructor
	int nrow();
	int ncol();
	T& operator()(int row, int col);
	friend ostream& operator<< <>(ostream&,  Matrix<T>&);
	friend ifstream& operator>> <>(ifstream&,  Matrix<T>&);

	private: 
	//vector<T> matrix_;
	int num_row_;
	int num_col_;
	vector<T> matrix_;
};

//------template definition-------

//constructor
template<typename T> Matrix<T>::Matrix()
:num_row_(0), num_col_(0),matrix_(0, 0)
{
}
//constructor
template<typename T> Matrix<T>::Matrix(int num_row, int num_col)
: num_row_(num_row), num_col_(num_col), matrix_(num_row_*num_col_, 0)
{
}

//Deconstructor	
template<typename T> Matrix<T>::~Matrix<T>()
{
}

//copy-constructor
template<typename T> Matrix<T>::Matrix(const Matrix<T>& m_)
:num_row_{m_.num_row_},num_col_{m_.num_col_}//,matrix_{m_.matrix_}
{
//	cout << " copy constructor matrix" << endl; 
}

// returns numbers of rows
template<typename T> int Matrix<T>::nrow(){
	return (num_row_);
}

// returns number of columns
template<typename T> int Matrix<T>::ncol(){
	return (num_col_);
}

// operator() returns entry at position (x,y)
template<typename T> T& Matrix<T>::operator()(int row, int col){

	//check if index out of range
	if ((row > num_row_) or (col > num_col_))
		 throw range_error("index out of range!!!");
	return(matrix_[((row-1)*(num_col_)) + col-1]);	
}

//output a matrix on ostream
template<typename T> ostream& operator<<(ostream& os ,Matrix<T>& m_ ){
	int counter = 0;
	for(int i = 0; i < m_.num_row_* m_.num_col_; ++i){
		os << m_.matrix_[i]<< "   ";
		counter ++;
		if (counter == m_.num_col_){
			counter = 0;
			os << "\n";
		}
	}	
	return os;
} 

// read data from a file 
template<typename T> ifstream& operator>>(ifstream& in ,Matrix<T>& m_ ){

	string line = "";
	m_.num_row_ = 0;
	m_.num_col_ = 0;
	int i = 0;
	while(getline(in, line)){
		m_.num_row_++;
			
		i =  count(line.begin(), line.end(), '\t') + 1;
		if (m_.num_col_ != 0 and i != m_.num_col_){
			throw range_error(" Wrong matrix input format!!!");
		}else{
			m_.num_col_ = i;
		}
	}

	//reset pointer of in to begin of the input
	in.clear();
	in.seekg(0,in.beg);
	
	T help; 
	Matrix<T> m(m_.num_row_, m_.num_col_);
	m_ = m;
	for(int i = 0; i < m_.num_row_* m_.num_col_; ++i){
		in >> help;
		m_.matrix_[i] = help;
	}
	return in;
}
#endif/*MATRIX_HPP*/
