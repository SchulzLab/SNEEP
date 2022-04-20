#ifndef PVALUE_HPP
#define PVALUE_HPP

#include <math.h>  
#include <string>
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include "Matrix_new.hpp"

using namespace std;
//double E = 1.2e-16;
double E = 1.2e-16;
double E2 = 1.0e-12;



class pvalue{

	public: 
	
	//constructors
	pvalue();
	pvalue(int size_kmer, double acc);
	~pvalue(); //deconstructor
	pvalue(const pvalue&); //copy constructor

	//vector<double> calculatePvalues(Matrix<double>& PWM, Matrix<double>& freq);
	vector<double> calculatePvalues(Matrix<double>& PWM, vector<double>& freq, Matrix<double>& background);

	//getter
	double getAcc();
	double getThres();
	int getKmerSize();

	private:
	int size_kmer_;
	double acc_;
	double thres_;	
	double calculateMaxThres(Matrix<double>& PWM);
};

//constructor
pvalue::pvalue()
:size_kmer_(0), acc_(1), thres_(1)
{
}

//constructor
pvalue::pvalue(int size_kmer, double acc)
:size_kmer_(size_kmer), acc_(acc), thres_(0)
{
}

//deconstructor
pvalue::~pvalue()
{
//	cout << "deconstructor" << endl;
}

//copy constructor
pvalue::pvalue(const pvalue& p_)
:size_kmer_{p_.size_kmer_}, acc_{p_.acc_},thres_{p_.thres_}
{
}


vector<double> pvalue::calculatePvalues(Matrix<double>& PWM, vector<double>& freq, Matrix<double>& background){

	thres_ = calculateMaxThres(PWM);
	int helper = (round((1/(acc_))*thres_))+1;
	vector<double> vec(helper ,0.0);
	double t = 0.0;
	double min = 1;
	//Matrix<double> DP(size_kmer_+1, helper);
	Matrix<double> DP((size_kmer_ * 4)+1, helper); // für jeden buchstaben eigene reihe 
	int counter = 1;
	int pos_DP = 0;
	/// iterate over length of kmer -> all prefixe of the kmer		
	for (int i = 0; i <= size_kmer_; ++i){
		// case Q_-1(t)
		if (i == 0){
			//DP(i+1,helper) = 1.0;
			DP(pos_DP+1,helper) = 1.0;
			pos_DP++;
		//case Q_i(t)
		}else{
			//iterate over all possible thresholds
			for (int k = 1; k <=helper; k++){
				if (k == 1){
					//determine min entry
					// so ueberspringen wir die positionen in die wir eh null eintragen wuerden
					for( int n = 1; n <= 4; n++){
						if (PWM(n,i) < min)
							min = PWM(n,i);			
					}
				}
				t+=acc_;		
				if (t >= min){	
					//iterate over column i+1  of PWM
					for( int j = 1; j <=4; ++j){  
						if (i == 1){
							pos_DP = 1;
						}
						//cout << "pos_DP: " << pos_DP << endl;
						double pos = round(helper- ((t- PWM(j,i))*(1/acc_)));
						if (t - PWM(j,i) < -(acc_/2)){ 
							continue;
						}else{
							if (pos_DP == 1){
									
								DP(pos_DP+j, helper - counter) = DP(pos_DP +j, helper-counter) + (DP(pos_DP, pos)*freq[j-1]);
							}else{
								//rausfinden welcher buchstabe vorhergelesen wurde 
								if (abs(DP(pos_DP, pos)) > E){ //vorheriger buchstabe A
									DP(pos_DP+4, helper - counter) = DP(pos_DP+4, helper-counter) + (DP(pos_DP, pos)*background(j, 1));
								}
								if (abs(DP(pos_DP+1, pos)) > E){ //vorheriger buchstabe C
									DP(pos_DP+5, helper - counter) = DP(pos_DP+5, helper-counter) + (DP(pos_DP+1, pos)*background(j, 2));
								}
								if (abs(DP(pos_DP+2, pos)) > E){ //vorheriger buchstabe G
									DP(pos_DP+6, helper - counter) = DP(pos_DP+6, helper-counter) + (DP(pos_DP+2, pos)*background(j,3));
								}
								if (abs(DP(pos_DP+3, pos)) > E){ //vorheriger buchstabe T
									DP(pos_DP+7, helper - counter) = DP(pos_DP+7, helper-counter) + (DP(pos_DP+3, pos)*background(j,4));
								}
							}
						}						
					}
				}
				counter++;
			}
		if (pos_DP > 1)
			pos_DP+=4;
		if (pos_DP == 1)
			pos_DP++;
//		cout << "DP:\n" << DP << endl;
		min = 1;
		counter = 1;
		t = 0;
		}
	}
//	cout << "check" << endl;
	//----------------------------------
	//checks if algorithm works correct

	//da totale whrscheinlichkeit muessen die wahrscheinlichkeiten einer zeile aufaddiert eins sein
	// und hier müsste man immer vier aufeinanderfolgende zeilen aufaddieren
	vector<double> sum(size_kmer_);
	int counter_sum = 0;
	int h = 0;
	for(int j = 2; j <= DP.nrow(); j++){
		for(int k = 1; k <= DP.ncol(); k++){
			sum[counter_sum]+= DP(j,k);
		}
		h++;
		if (h == 4){
			counter_sum++;
			h = 0;
		}
	}
	bool check = false;
	for(auto& m: sum){
		if (abs((m-1.0)) > E2){
			cout << m << endl;
			cout << m  - 1<< endl;
			cout << "oh nein!!! da stimmt was nicht!!!!" << endl;
			check = true;
		}
	}
	if (check == true){
		for (auto& m: sum)
			cout << m << " ";
		cout << endl;
	}
	//-------------------------------

	// fuellen des vectors der cumulative wahrscheinlichkeit enthaelt, entspricht dem pvalue
	// und hier auch die letzen vier zeilen aufaddieren
	int last_row = (size_kmer_*4) +1;
	for(int i= 1; i<= helper; i++){
		if (i == 1){
			vec[0] = DP(last_row, i) + DP(last_row-1, i) + DP(last_row -2, i) + DP(last_row-3, i);
			
		}else{
			vec[i-1] = vec[i-2] + DP(last_row, i) + DP(last_row - 1, i) + DP(last_row - 2, i) + DP(last_row - 3, i);
		}
	}
	return vec;	
}

double pvalue::calculateMaxThres(Matrix<double>& PWM){

	double t = 0.0;
	double max = 0.0;
	for(int i = 1; i <= PWM.ncol(); i++){
		for(int j = 1; j <= 4; j++){
			if (PWM(j,i) > max){
				max = PWM(j,i);
			}	
		}
		t += max; 
		max = 0.0;
	}	
	return t;
}

//getter
double pvalue::getAcc(){
	return (acc_);	
}

double pvalue::getThres(){
	return (thres_);
}

int pvalue::getKmerSize(){
	return(size_kmer_);
}

#endif/*PVALUE_HPP*/
