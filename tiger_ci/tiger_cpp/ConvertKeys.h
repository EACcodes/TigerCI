/* 
 * File:   ConvertKeys.h
 * Author: Francis Ricci
 *
 * Created on June 10, 2014, 11:00 AM
 */

#ifndef CONVERTKEYS_H
#define	CONVERTKEYS_H

#include <string>

#include "InputReader.h"

using namespace std;

// external and without C++ "name mangling" to allow for calling by Fortran
extern "C"
{
    // get int val for global variable c_name
    void c_get_global_int(const char *c_name, int64_t& val);
    
    // get double val for global variable c_name
    void c_get_global_double(const char *c_name, double& val);
    
    // get bool val for global variable c_name
    void c_get_global_bool(const char *c_name, int64_t& val);
    
    // get string val for global variable c_name
    void c_get_global_string(const char *c_name, char *val, int64_t& length);
    
    // get int array val for vector num of global variable c_name
    void c_get_global_int_vec(const char *c_name, int64_t& num, int64_t *val);
    
    // get number of cylinders for a given orbital
    void c_num_cylinders(int64_t& orb, int64_t& num_cyl);
    
    // get cylinders on a given orbital
    void c_get_cylinders(int64_t& orb, int64_t *val);
}

class ConvertKeys {
public:
    // default constructor
    ConvertKeys(){};
    
    // convert the keys in key_map into variable name-value pairs,
    // store in an unordered map
    void keys_to_globals(InputReader& key_map);
    string get_my_string(string name) const;
    int get_global_int(string name) const;
    bool get_global_bool(string name) const;
    string get_global_string(string name) const;
    double get_global_double(string name) const;
    vector<vector<int>> get_global_int_vecs(string name) const;
    vector<vector<double>> get_global_double_vecs(string name) const;

    // function to add global variables that old TigerCI needs
    void add_global_int(string name, int val);
    void add_global_bool(string name, bool val);
    
private:
    // setup nonlocal parameters
    void setup_nonlocal();
    // fill all references with zeros
    void fill_references();
    // check input parameters
    void check_input();
    // check for unique references
    void unique_ref_check(vector<vector<int>>& refs);
    unordered_map<string, string> my_strings;
};
#endif	/* CONVERTKEYS_H */
