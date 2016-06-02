/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * File:   InputReader.cpp
 * Authors: Francis Ricci, David Krisiloff, Johannes M Dieterich
 *
 * \todo Add sanity checks on different inputs
 */

#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <algorithm>

#include "Keyword.h"
#include "StringUtils.h"
#include "InputReader.h"

using namespace std;

InputReader::InputReader() {
	/*
	 * This is where we define each of the input_keys in the TigerCI
	 * input file, their default values, nonlocal values, etc..
	 */

	// Basic molecule setup
	input_keys["NUMBER OF ELECTRONS"] = shared_ptr<Keyword>
	(new RequiredIntKeyword("NUMBER OF ELECTRONS"));
	input_keys["NUMBER OF ORBITALS"] = shared_ptr<Keyword>
	(new RequiredIntKeyword("NUMBER OF ORBITALS"));
	input_keys["NUMBER OF FROZEN ORBITALS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("NUMBER OF FROZEN ORBITALS", 0));
	input_keys["NUMBER OF INACTIVE ORBITALS"] = shared_ptr<Keyword>
	(new RequiredIntKeyword("NUMBER OF INACTIVE ORBITALS"));
	input_keys["NUMBER OF ACTIVE ORBITALS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("NUMBER OF ACTIVE ORBITALS",0));
	input_keys["NUMBER OF BASIS FUNCTIONS"] = shared_ptr<Keyword>
	(new RequiredIntKeyword("NUMBER OF BASIS FUNCTIONS"));
	input_keys["SPIN MULTIPLICITY"] = shared_ptr<Keyword>
	(new RequiredIntKeyword("SPIN MULTIPLICITY"));
	input_keys["NUMBER OF REFERENCES"] = shared_ptr<Keyword>
	(new RequiredIntKeyword("NUMBER OF REFERENCES"));

	// Tolerance and thresholds
	input_keys["ENERGY TOLERANCE"] = shared_ptr<Keyword>
	(new DoubleKeyword("ENERGY TOLERANCE", 1.0e-6));
	input_keys["RESIDUAL TOLERANCE"] = shared_ptr<Keyword>
	(new DoubleKeyword("RESIDUAL TOLERANCE", 1.0e-6));
	input_keys["INTEGRAL THRESHOLD"] = shared_ptr<Keyword>
	(new DoubleKeyword("INTEGRAL THRESHOLD", 1.0e-12));
	input_keys["AO INTEGRAL THRESHOLD"] = shared_ptr<Keyword>
	(new DoubleKeyword("AO INTEGRAL THRSHOLD", 1.0e-12));
	input_keys["OCCUPATION THRESHOLD"] = shared_ptr<Keyword>
	(new DoubleKeyword("OCCUPATION THRESHOLD", 0.8));
	input_keys["VIRTUAL OCCUPATION THRESHOLD"] = shared_ptr<Keyword>
	(new DoubleKeyword("VIRTUAL OCCUPATION THRESHOLD",0.8));
	input_keys["CD THRESHOLD"] = shared_ptr<Keyword>
	(new DoubleKeyword("CD THRESHOLD", 1.0e-12));

	// Allocation-related variables
	input_keys["CD MEMORY"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("CD MEMORY", 1000));
	input_keys["INT MEMORY"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("INT MEMORY", 500));
	input_keys["AIJK MEMORY"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("AIJK MEMORY", 10));
	input_keys["DAVIDSON HARD DRIVE SPACE"] = shared_ptr<Keyword>
	(new DoubleKeyword("DAVIDSON HARD DRIVE SPACE",30.0));
	input_keys["STRIDE SIZE SIGMA"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("STRIDE SIZE SIGMA", 32));
	input_keys["STRIDE SIZE CI"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("STRIDE SIZE CI", 32));
	input_keys["IOBUFFBLOCKINTS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("IOBUFFBLOCKINTS", 10));
	input_keys["IOBUFFMAXMEMINTS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("IOBUFFMAXMEMINTS", 500));
	input_keys["IOBUFFMAXMEMCD"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("IOBUFFMAXMEMCD", 500));
	input_keys["IOBUFFMAXMEMSEGS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("IOBUFFMAXMEMSEGS", 50));
	input_keys["IOBUFFBLOCKSEGS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("IOBUFFBLOCKSEGS", 10));

	// Geometric parameters
	input_keys["WP DEFAULT RADIUS"] = shared_ptr<Keyword>
	(new DoubleKeyword("WP DEFAULT RADIUS", 2.65));
	input_keys["TOV OCCUPIED DEFAULT RADIUS"] = shared_ptr<Keyword>
	(new DoubleKeyword("TOV OCCUPIED DEFAULT RADIUS", 2.65));
	input_keys["TOV VIRTUAL DEFAULT RADIUS"] = shared_ptr<Keyword>
	(new DoubleKeyword("TOV VIRTUAL DEFAULT RADIUS", 2.65));
	input_keys["WP RADIUS MULTIPLIER"] = shared_ptr<Keyword>
	(new DoubleKeyword("WP RADIUS MULTIPLIER", 1.7));
	input_keys["TOV OCCUPIED RADIUS MULTIPLIER"] = shared_ptr<Keyword>
	(new DoubleKeyword("TOV OCCUPIED RADIUS MULTIPLIER", 1.95));
	input_keys["TOV VIRTUAL RADIUS MULTIPLIER"] = shared_ptr<Keyword>
	(new DoubleKeyword("TOV VIRTUAL RADIUS MULTIPLIER", 2.0));
	input_keys["CYLINDER RADIUS"] = shared_ptr<Keyword>
	(new DoubleKeyword("CYLINDER RADIUS", 2.0));
	input_keys["TOV CYLINDER RADIUS"] = shared_ptr<Keyword>
	(new DoubleKeyword("TOV CYLINDER RADIUS", 2.0));

	// Flags
	input_keys["VALENCE CI FLAG"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("VALENCE CI FLAG", 0));
	input_keys["RESTART FLAG"] = shared_ptr<Keyword>
	(new DefaultBoolKeyword("RESTART FLAG", false));
	input_keys["NATURAL ORBITAL FLAG"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("NATURAL ORBITAL FLAG", 0));
	input_keys["REFERENCE CI FLAG"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("REFERENCE CI FLAG", 0));
	input_keys["ACPF FLAG"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("ACPF FLAG", 0));
	input_keys["NONLOCAL"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("NONLOCAL", 0));
	input_keys["CD RESTART"] = shared_ptr<Keyword>
	(new DefaultBoolKeyword("CD RESTART", false));
	input_keys["PEDERSEN CD"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("PEDERSEN CD"));
	input_keys["DENSITY FITTING"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("DENSITY FITTING"));
	input_keys["CARTESIAN"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("CARTESIAN"));
	input_keys["SKIP INACTIVE PM"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("SKIP INACTIVE PM"));
	input_keys["SKIP LOVOS"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("SKIP LOVOS"));

	// Flags which take no arguments
	input_keys["INTEGRAL DIRECT MODE"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("INTEGRAL DIRECT MODE"));
	input_keys["FULLY INTEGRAL DIRECT MODE"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("FULLY INTEGRAL DIRECT MODE"));
        input_keys["NO DIRECT FOUR INTERNAL"] = shared_ptr<Keyword>
	(new DefaultBoolKeyword("NO DIRECT FOUR INTERNAL",false));
	input_keys["LOW MEMORY DIRECT"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("LOW MEMORY DIRECT"));
	input_keys["SUPER LOW MEMORY DIRECT"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("SUPER LOW MEMORY DIRECT"));
	input_keys["CD VECTORS IN MEMORY"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("CD VECTORS IN MEMORY"));
	input_keys["SPHERE PRINT"] = shared_ptr<Keyword>
	(new DefaultBoolKeyword("SPHERE PRINT",0));
	input_keys["NO ACPF ROOT FOLLOW"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("NO ACPF ROOT FOLLOW"));
	input_keys["ROOT FOLLOW HELP"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("ROOT FOLLOW HELP"));
	input_keys["DONTSTOREINTS"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("DONTSTOREINTS"));
	input_keys["DONTSTORECDVECS"] = shared_ptr<Keyword>
	(new FlagBoolKeyword("DONTSTORECDVECS"));
    input_keys["EMBEDDING"] = shared_ptr<Keyword>
    (new FlagBoolKeyword("EMBEDDING"));

	// Other parameters
	input_keys["NUM ROOTS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("NUM ROOTS", 1));
	input_keys["CUSTOM G VAL"] = shared_ptr<Keyword>
	(new DoubleKeyword("CUSTOM G VAL", 1.0));
	input_keys["NUM THREADS"] = shared_ptr<Keyword>
	(new DefaultIntKeyword("NUM THREADS", 1));

	//Vector parameters
	input_keys["REFERENCE OCCUPATIONS"] = shared_ptr<Keyword>
	(new RefVecsKeyword("REFERENCE OCCUPATIONS"));
	input_keys["CYLINDERS"] = shared_ptr<Keyword>
	(new CylVecsKeyword("CYLINDERS"));
    input_keys["GRID DIMENSIONS"] = shared_ptr<Keyword>
    (new IntVecsKeyword("GRID DIMENSIONS"));
    input_keys["LATTICE VECTORS"] = shared_ptr<Keyword>
    (new DoubleVecsKeyword("LATTICE VECTORS"));
    input_keys["SHIFT VECTOR"] = shared_ptr<Keyword>
    (new DefaultDoubleVecsKeyword("SHIFT VECTOR",{{0.0,0.0,0.0}}));

	// Files and directories
	input_keys["SCRATCH DIRECTORY"] = shared_ptr<Keyword>
	(new RequiredStringKeyword("SCRATCH DIRECTORY"));

	// Stand-alone related variables
        input_keys["XYZ FILE"] = shared_ptr<Keyword>
	(new DefaultStringKeyword("XYZ FILE", "")); // not used anymore, deprecated!
	input_keys["ORB FILE"] = shared_ptr<Keyword>
	(new RequiredStringKeyword("ORB FILE"));
	input_keys["BASIS SET"] = shared_ptr<Keyword>
	(new RequiredStringKeyword("BASIS SET"));
	input_keys["AUXILIARY BASIS SET"] = shared_ptr<Keyword>
	(new DefaultStringKeyword("AUXILIARY BASIS SET", ""));
    input_keys["EMBEDDING POTENTIAL FILE"] = shared_ptr<Keyword>
    (new DefaultStringKeyword("EMBEDDING POTENTIAL FILE","embpot.dat"));
}

// Read and store input parameters from string file_name
void InputReader::read_tiguar_input(string file_name) {
	ifstream input;
	input.open(file_name);

	string line;
	shared_ptr<Keyword> key = nullptr;
	unordered_map<string, shared_ptr < Keyword >> ::iterator loc_finder;

	while (getline(input, line)) {
		line = trim(line);
		if (line == "END OF INPUT")
			break;
		loc_finder = input_keys.find(line);
		if (loc_finder != input_keys.end()) {
			// this is a keyword
			key = loc_finder->second;
			// check that we haven't set it yet
			if (key->isValueSet())
				throw runtime_error("Duplicate keyword: " + key->getName());
			// read in the value(s)
			while (!key->isValueSet()){
				key->parse(input);
			}
		}
		else
			throw runtime_error("Unknown keyword: " + line);
	}

	// check that we haven't missed any required keywords
	for (auto key_val : input_keys) {
		key = key_val.second;
		if (key->isRequired() && !key->isValueSet())
			throw runtime_error("The required keyword " + key_val.first +
					" was not found in the input file");
	}

	// check references
	int num_ref = getIntKeyword("NUMBER OF REFERENCES");
	int num_electrons = getIntKeyword("NUMBER OF ELECTRONS");
	int num_inactive = getIntKeyword("NUMBER OF INACTIVE ORBITALS");
	int num_frozen = getIntKeyword("NUMBER OF FROZEN ORBITALS");
	shared_ptr<RefVecsKeyword> references = dynamic_pointer_cast<RefVecsKeyword>
	(input_keys.at("REFERENCE OCCUPATIONS"));
	references->RefVecsCheck(num_ref, num_electrons, num_inactive, num_frozen);

    // check embedding data
    bool embedding = getBoolKeyword("EMBEDDING");
    if (embedding == true) {
        for (auto key_val : input_keys) {
            key = key_val.second;
            if (key->getName() == "GRID DIMENSIONS" && !key->isValueSet()) 
                throw invalid_argument("The keyword " + key_val.first + 
                    " required for embedding calculations was not found in the input file");
            else if (key->getName() == "LATTICE VECTORS" && !key->isValueSet())
                throw invalid_argument("The keyword " + key_val.first + 
                      " required for embedding calculations was not found in the input file");
        }
    }
}

// print out all keywords and their values
void InputReader::output_keywords() {
	
    cout << endl;
    cout << "********************************************************************************" << endl;
    cout << "**********                                                            **********" << endl;
    cout << "                        CONFIGURATION OF THIS CALCULATION                       " << endl;
    cout << endl;
    cout << setw(40) << left << "Keyword" << "Value" << endl;
    cout << setw(40) << left << "-------" << "-----" << endl;


    // let's make the output as easy to parse as possible
    // by writing out the keywords in alphabetical order
    vector<string> keys;
    for (auto key_val : input_keys)
        keys.push_back(key_val.first);
    sort(keys.begin(), keys.end());

    for (string key : keys) {
        // keyword name
        cout << setw(40) << left << key;

        // keyword value
        auto val = input_keys[key];
        // check if int keyword
        shared_ptr<IntKeyword> ik = dynamic_pointer_cast<IntKeyword> (val);
        if (ik) {
            cout << ik->getIntKeyword();
        }

        // check if string keyword
        shared_ptr<StringKeyword> sk = dynamic_pointer_cast<StringKeyword> (val);
        if (sk) {
            cout << sk->getStringKeyword();
        }

        // check if bool keyword
        shared_ptr<BoolKeyword> bk = dynamic_pointer_cast<BoolKeyword> (val);
        if (bk) {
            cout << bk->getBoolKeyword();
        }

        // check if double keyword
        shared_ptr<DoubleKeyword> dk = dynamic_pointer_cast<DoubleKeyword> (val);
        if (dk) {
            cout << dk->getDoubleKeyword();
        }

        // check if ref keyword
        shared_ptr<RefVecsKeyword> rk =
                dynamic_pointer_cast<RefVecsKeyword> (val);
        if (rk) {
            vector<vector<int>> refs = rk->getRefVecsKeyword();
            bool init = false; // special formatting for first line
            for (vector<int> vec : refs) {
                if (!init)
                    init = true;
                else
                    cout << "                                        ";
                for (int i : vec) {
                    cout << i;
                }
                cout << endl;
            }
            continue;
        }

        // check if cyl keyword
        shared_ptr<CylVecsKeyword> ck =
                dynamic_pointer_cast<CylVecsKeyword> (val);
        if (ck) {
            vector<Chain> chains = ck->getCylVecsKeyword();
            if (!(ck->isValueSet())) {
                cout << "Not defined" << endl;
                continue;
            }
            bool init = false; //special formatting for first line
            for (Chain cyls : chains) {
                if (!init)
                    init = true;
                else
                    cout << "                                        ";
                cyls.print();
                cout << endl;
            }
            continue;
        }

        // check if IntVecs keyword
        shared_ptr<IntVecsKeyword> iv =
                dynamic_pointer_cast<IntVecsKeyword> (val);
        if (iv) {
            vector<vector<int>> vecs = iv->getIntVecsKeyword();
            if (!(iv->isValueSet())) {
                cout << "Not defined" << endl;
                continue;
            }
            bool init = false; // special formatting for first line
            for (vector<int> vec : vecs) {
                if (!init)
                    init = true;
                else
                    cout << "                                        ";
                for (int i : vec) {
                    cout << i << "  ";
                }
                cout << endl;
            }
            continue;
        }

        // check if DoubleVecs keyword
        shared_ptr<DoubleVecsKeyword> dv =
                dynamic_pointer_cast<DoubleVecsKeyword> (val);
        if (dv) {
            vector<vector<double>> vecs = dv->getDoubleVecsKeyword();
            if (!(dv->isValueSet())) {
                cout << "Not defined";
                if (dv->hasDefault()) {
                    cout << ", use default:" << endl;
                    cout << "                                        ";
                } else {
                    cout << endl;
                    continue;
                }
            }
            bool init = false; // special formatting for first line
            for (vector<double> vec : vecs) {
                if (!init)
                    init = true;
                else
                    cout << "                                        ";
                for (double x : vec) {
                    cout << x << "  ";
                }
                cout << endl;
            }
            continue;
        }
        cout << endl;
    }
    
    cout << "**********                                                            **********" << endl;
    cout << "********************************************************************************" << endl;
}

// return value for string keyword given by name
string InputReader::getStringKeyword(string name){
	shared_ptr<Keyword> k = input_keys.at(name);
	shared_ptr<StringKeyword> sk =
			dynamic_pointer_cast<StringKeyword>(k);
	return sk->getStringKeyword();
}

// return value for int keyword given by name
int InputReader::getIntKeyword(string name){
	shared_ptr<Keyword> k = input_keys.at(name);
	shared_ptr<IntKeyword> ik = dynamic_pointer_cast<IntKeyword>(k);
	return ik->getIntKeyword();
}

// return value for bool keyword given by name
bool InputReader::getBoolKeyword(string name){
	shared_ptr<Keyword> k = input_keys.at(name);
	shared_ptr<BoolKeyword> bk = dynamic_pointer_cast<BoolKeyword>(k);
	return bk->getBoolKeyword();
}

// return value for double keyword given by name
double InputReader::getDoubleKeyword(string name){
	shared_ptr<Keyword> k = input_keys.at(name);
	shared_ptr<DoubleKeyword> dk =
			dynamic_pointer_cast<DoubleKeyword>(k);
	return dk->getDoubleKeyword();
}

// return value for ref keyword given by name
vector<vector<int>> InputReader::getRefVecsKeyword(string name){
	shared_ptr<Keyword> k = input_keys.at(name);
	shared_ptr<RefVecsKeyword> rk = dynamic_pointer_cast<RefVecsKeyword>(k);
	return rk->getRefVecsKeyword();
}


// return value for cyl keyword given by name
vector<Chain> InputReader::getCylVecsKeyword(string name){
	shared_ptr<Keyword> k = input_keys.at(name);
	shared_ptr<CylVecsKeyword> ck = dynamic_pointer_cast<CylVecsKeyword>(k);
	return ck->getCylVecsKeyword();
}

// return value for IntVecs keyword given by name
vector<vector<int>> InputReader::getIntVecsKeyword(string name){
  shared_ptr<Keyword> k = input_keys.at(name);
  shared_ptr<IntVecsKeyword> iv = dynamic_pointer_cast<IntVecsKeyword>(k);
  return iv->getIntVecsKeyword();
}

// return value for DoubleVecs keyword given by name
vector<vector<double> > InputReader::getDoubleVecsKeyword(string name){
    shared_ptr<Keyword> k = input_keys.at(name);
    shared_ptr<DoubleVecsKeyword> dv = dynamic_pointer_cast<DoubleVecsKeyword>(k);
  return dv->getDoubleVecsKeyword();
}

// return a vector containing all keyword names
vector<string> InputReader::getNames() {
	vector<string> names;
	for (pair<string, shared_ptr<Keyword>> key : input_keys)
		names.push_back(key.first);
	return names;
}
