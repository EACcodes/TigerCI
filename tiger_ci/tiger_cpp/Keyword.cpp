/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * File:   Keyword.cpp
 * Authors: Francis Ricci and David Krisiloff
 *
 */

#include <string>
#include <sstream>
#include <iostream>
#include <exception>
#include <iomanip>
#include <fstream>
#include <stdexcept>

#include "Keyword.h"
#include "StringUtils.h"

using namespace std;

// Int keywords consist of a keyword plus an int value
// Must provide name, default value, and required status on construction
IntKeyword::IntKeyword(string name, int default_val, bool required):
  Keyword(name, required), _val(default_val)  {}

// return the value of this int keyword
int IntKeyword::getIntKeyword(){
    return _val;
}

// Specific parsing for int keywords
void IntKeyword::parse(ifstream& input){
    string s;
    getline(input, s);
    s = trim(s);
   
#ifdef __FreeBSD__
    int value = -1;
    stringstream ss(s);
    ss >> value;
#else
    int value = stoi(s);
#endif
    if (value<0)
        throw runtime_error("Detected a non-positive integer for keyword " +
                _name + '\n' + "Value read in was = " + s);
    _val = value;
    _finished_reading = true;
}

// constructor requires the name of the keyword and its default value
DefaultIntKeyword::DefaultIntKeyword(string name, int default_value):
  IntKeyword(name, default_value, false) {}

// constructor requires only the keyword name
RequiredIntKeyword::RequiredIntKeyword(string name):
  IntKeyword(name, -1, true) {}

// constructor requires the name of the keyword
StringKeyword::StringKeyword(string name, string default_val, bool required):
  Keyword(name, required), _val(default_val) {}

// specific parsing for string keywords
void StringKeyword::parse(ifstream& input){
    string s;
    getline(input, s);
    s = trim(s);
    
    _val = s;
    _finished_reading = true;
}

// return the value of the string keyword
string& StringKeyword::getStringKeyword(){
    string& value = _val;
    return value;
}     

DefaultStringKeyword::DefaultStringKeyword(string name, string default_val):
  StringKeyword(name, default_val, false) {}

RequiredStringKeyword::RequiredStringKeyword(string name):
  StringKeyword(name, "", true) {}

// constructor requires the name of the keyword, its default value, and
// whether or not it is a required input
BoolKeyword::BoolKeyword(string name, bool default_val, bool required):
  Keyword(name, required), _val(default_val) {}

// constructor requires keyword name
RequiredBoolKeyword::RequiredBoolKeyword(string name):
  BoolKeyword(name, 0, true) {}

// constructor requires keyword name and default value
DefaultBoolKeyword::DefaultBoolKeyword(string name, bool default_val):
  BoolKeyword(name, default_val, false) {}

// constructor requires keyword name
FlagBoolKeyword::FlagBoolKeyword(string name):
  BoolKeyword(name, false, false) {}

// specific parsing for boolean keywords
void BoolKeyword::parse(ifstream& input){
    string s;
    getline(input, s);
    s = trim(s);
    
    /* Boolean keywords are inputed as 0 or 1*/
#ifdef __FreeBSD__
    int value = -1;
    stringstream ss(s);
    ss >> value;
#else
    int value = stoi(s);
#endif
    if (!((value==0)|| (value==1)))
        throw runtime_error("Detected an input of not 0 or 1 for a boolean"
                + string(" keyword - ") + _name);
    if (value==1){
        _val = true;
        
   }
   else{
       _val = false;
   }
   _finished_reading = true;   
}

// specific parsing for flag bool keywords
void FlagBoolKeyword::parse(ifstream& input) {
    // Boolean keywords given as flags have no arguments!
    char c = input.peek();
    if (isdigit(c))
        throw runtime_error("Detected an argument for a flag taking no" +
                string(" inputs - ") + _name);
    _val = true;
    _finished_reading = true;
}
       
// return value of the boolean keyword
bool BoolKeyword::getBoolKeyword(){
    return _val;
}
      
// constructor requires keyword name and default value
DoubleKeyword::DoubleKeyword(string name, double default_val):
  Keyword(name, false), _val(default_val) {}

// return the value of the double keyword
double DoubleKeyword::getDoubleKeyword(){
    return _val;
}

// specific parsing for double keywords
void DoubleKeyword::parse(ifstream& input){
    string s;
    getline(input, s);
    s = trim(s);
    stringstream ss(s);
    ss >> _val;
    _finished_reading = true;
    
}

// constructor requires keyword name
RefVecsKeyword::RefVecsKeyword(string key):
  Keyword(key, true) {}
        
// returns a reference to the vector of reference vectors
vector<vector<int>>& RefVecsKeyword::getRefVecsKeyword(){
    vector<vector<int>>& res = _val;
    return res;
}
      
// specific parsing for reference keywords
void RefVecsKeyword::parse(ifstream& input){
    char c = input.peek();

    // if there are no more references to read, we are finished
    if (!isdigit(c)){
        if (_val.size() == 0)
            throw runtime_error("Must provide vector arguments for " + _name);
        else{
            _finished_reading = true;
            return;
        }
    }
    
    string s;
    getline(input, s);
    s = trim(s);
    
    vector<int> ref;
    
    // extract one int at a time from the reference
    while(s.length() != 0){
        int i = static_cast<int>(s[0]-'0');
        ref.push_back(i);
        s = s.substr(1);
        while (isspace(s[0]))
            s = s.substr(1);
    }
    // add reference vector to the vector of vectors
    _val.push_back(ref);
}

// run a check for reference consistency, given several other input parameters
void RefVecsKeyword::RefVecsCheck(int num_ref, int num_electrons,
                                  int num_inactive, int num_frozen){
    // global checks
    if (num_ref == 0)
        throw runtime_error("Number of references must not be zero.");
    
    else if (num_electrons <= 0)
        throw runtime_error("Number of electrons must be greater than zero.");
    
    else if (num_inactive < num_frozen)
        throw runtime_error("Number of frozen orbitals should be less than " +
                string("the number of inactive orbitals."));
    
    // iterate through all references
    int ref_num {0};
    for (vector<int> vec : _val){
        int orb_num {0};
        int elec_num {0};
        for (int i : vec){
            if ((orb_num < num_inactive) && (i != 2))
                throw runtime_error("Inactive orbitals must be doubly occupied.");
            if (i < 0 || i > 2){
                ostringstream oss;
                oss << "Invalid reference occupation - " << i;
                throw runtime_error(oss.str());
            }
            orb_num++;
            elec_num += i;
        }
        ref_num++;
        if (elec_num != num_electrons){
            ostringstream oss;
            oss << "Incorrect number of electrons in reference " << num_ref;
            throw runtime_error(oss.str());
        }
    }
    if (ref_num != num_ref)
        throw runtime_error("Incorrect number of reference vectors provided");
}

// constructor requires keyword name
CylVecsKeyword::CylVecsKeyword(string key):
  Keyword(key, false) {}
        
// returns a reference to the vector of cylinder coordinate vectors
vector<Chain>& CylVecsKeyword::getCylVecsKeyword(){
    vector<Chain>& res = _val;
    return res;
}
     
// specific parsing for cylinder coordinates
void CylVecsKeyword::parse(ifstream& input){
    char c = input.peek();

    if (!isdigit(c)){
        if (_val.size() == 0)
            throw runtime_error("Must provide vector arguments for " + _name);
        else{
            _finished_reading = true;
            return;
        }
    }
    
    string s;
    getline(input, s);
    s = trim(s);
    
    istringstream line_input(s + ',');
    
    // split out into individual cylinders
    string sub;
    getline(line_input, sub, ':');
#ifdef __FreeBSD__
    int orbital = -1;
    stringstream ss(sub);
    ss >> orbital;
#else
    int orbital = stoi(sub);
#endif
    Chain orb_chain(orbital-1);
    
    getline(line_input, sub, ',');
    while (!sub.empty()){
#ifdef __FreeBSD__
        sub = trim(sub);
        int first;
        int second;
        stringstream ss(sub);
        ss >> first;
        ss >> second;
#else
        string::size_type sz;
        int first = stoi(sub, &sz);
        int second = stoi(sub.substr(sz));
#endif
        
        Cylinder cyl = Cylinder(first, second);
        orb_chain.addCylinder(cyl);
        
        getline(line_input, sub, ',');
    }
    
    _val.push_back(orb_chain);
}

// constructor requires keyword name
IntVecsKeyword::IntVecsKeyword(string key):
    Keyword(key, false) {}

// returns a reference to the vector of integer vectors
vector<vector<int>>& IntVecsKeyword::getIntVecsKeyword(){
    vector<vector<int>>& res = _val;
    return res;
}

// specific parsing for integer vectors
void IntVecsKeyword::parse(ifstream& input){
    char c = input.peek();

    // if there are no more vectors to read, we are finished
    if (!isdigit(c)){
        if (_val.size() == 0)
            throw invalid_argument("Must provide vector arguments for " + _name);
        else{
            _finished_reading = true;
            return;
        }
    }

    string s;
    getline(input, s);
    s = trim(s);

    vector<int> vec;
    // split out into individual integers
    istringstream is( s );
    int n;
    while (is >> n) {
        vec.push_back(n);
    }

    _val.push_back(vec);
}

// constructor requires keyword name
DoubleVecsKeyword::DoubleVecsKeyword(string key):
    Keyword(key, false), _default(false) {}

// constructor requires keyword name
DoubleVecsKeyword::DoubleVecsKeyword(string key, vector<vector<double>> default_val):
    Keyword(key, false),_val(default_val), _default(true)  {}

// returns a reference to the matrix of doubles
vector<vector<double>>& DoubleVecsKeyword::getDoubleVecsKeyword(){
    vector<vector<double>>& vecs = _val;
    return vecs;
}

// specific parsing for double matrices
void DoubleVecsKeyword::parse(ifstream& input){
    char c = input.peek();

    // account for possibility of negative values -> peek second character if line starts with '-'
    if (c == '-') {
        input.get(c);
        char c1 = input.peek();
        input.putback(c);
        c=c1;
    }

    // if there are no more vectors to read, we are finished
    if (!isdigit(c)){
        if (_val.size() == 0)
            throw invalid_argument("Must provide vector arguments for " + _name);
        else{
            _finished_reading = true;
            return;
        }
    }

    string s;
    getline(input, s);
    s = trim(s);

    vector<double> vec;

    istringstream is( s );
    double n;
    while (is >> n) {
        vec.push_back(n);
    }

    _val.push_back(vec);
}

// Constructore when default value exists
DefaultDoubleVecsKeyword::DefaultDoubleVecsKeyword(string key, vector<vector<double>> default_val):
    DoubleVecsKeyword(key, default_val) {}
