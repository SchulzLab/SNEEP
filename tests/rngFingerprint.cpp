/*
* Prints random numbers drawn as in SNEEP's background sampling (std::mt19937 + std::uniform_int_distribution<int>).
* mt19937 is fully specified by the C++ standard, but the algorithm of uniform_int_distribution is up to the standard
* library (libc++, libstdc++, MSVC, and it may change between versions), so the same seed can give other random SNPs.
* runRegressionTests.sh uses the checksum of this output as a fingerprint of the toolchain: if it differs from the
* one stored with the reference, the sampling-dependent outputs are not compared.
*/
#include <iostream>
#include <random>

using namespace std;

int main(){
	int ranges[] = {1, 2, 6, 99, 1000, 65535, 1234567}; // upper bounds as used in SNEEP: alleles, tie-breaks, pool sizes
	for (int seed = 0; seed < 5; ++seed){
		for (int upper : ranges){
			mt19937 generator(seed * 1500 + upper);
			uniform_int_distribution<int> distribution(0, upper);
			for (int i = 0; i < 20; ++i){
				cout << distribution(generator) << ' ';
			}
			cout << '\n';
		}
	}
	return 0;
}
