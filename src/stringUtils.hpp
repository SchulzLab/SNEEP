#ifndef STRINGUTILS_HPP
#define STRINGUTILS_HPP

#include <string>
#include <stdexcept>

using namespace std;

/*
* returns the part of line before the first delim (or, if there is none, before the first '\n')
* and removes it together with the delimiter from line
*/
inline string getToken(string& line, char delim){
	size_t pos = line.find(delim);
	if (pos == string::npos){
		pos = line.find('\n');
	}
	if (pos == string::npos){
		throw invalid_argument("invalid file format:" + line);
	}
	string token = line.substr(0, pos);
	line.erase(0, pos + 1);
	return token;
}

#endif/*STRINGUTILS_HPP*/
