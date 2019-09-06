/*
 * inputParser.h
 *
 *  Created on: 25 Jun 2016
 *      Author: osboxes
 */

#ifndef INPUTPARSER_H_
#define INPUTPARSER_H_


#include <map>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cstdlib>

class InputParser {
public:
	std::string programName;
	InputParser(int argc, char *argv[], int startWith = 1) {
		if(startWith < 1 || startWith >= argc) std::runtime_error("InputParser: invalid startWith parameter.");
		programName = std::string(argv[0]);
		std::string opt;
		for(int i = startWith; i < argc; ++i) {
			std::string s(argv[i]);

			if(s.length() >= 2 && s[0]=='-' && !isNumber(s)) {
				opt = s.substr(1, s.length());
				inp[opt] = "";
				continue;
			}

			if(!opt.empty()) {
				inp[opt] = s;
				opt.clear();
			}
		}
	}

	std::istringstream& operator[](const std::string &key) {
		retVal.clear();
		try {
			retVal.str(inp.at(key));
		}catch(...) {
			throw std::runtime_error("InputParser: option \""+key+"\" is missing.");
		}
		return retVal;
	}

	bool exists(const std::string &key) {
		bool tmp = false;
		try {
			inp.at(key);
			tmp = true;
		}catch(...) {}
		return tmp;
	}

private:
	std::map<std::string,std::string> inp;
	std::istringstream retVal;

	bool isNumber(const std::string& s) {
	    char* p;
	    strtod(s.c_str(), &p);
	    return *p == 0;
	}
};


#endif /* INPUTPARSER_H_ */
