/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * File:   ConvertKeys.cpp
 * Author: Francis Ricci and Johannes M Dieterich
 */

#include <string.h>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <map>
#include <memory>
#include <stdexcept>

#include "ConvertKeys.h"
#include "Keyword.h"
#include "InputReader.h"
#include "Cylinder.h"

// Unordered maps for the storage of name-value pairs.
// Must be global to allow for calling by "C" routines from Fortran.
static unordered_map<string, int> global_ints;
static unordered_map<string, double> global_doubles;
static unordered_map<string, bool> global_bools;
static unordered_map<string, string> global_strings;
static unordered_map<string, vector<vector<int>>> global_int_vecs;
static unordered_map<string, vector<vector<double>>> global_double_vecs;
static map<int, Chain> cyl_endpoints;

using namespace std;

// convert the keys in key_map into variable name-value pairs,
// store in an unordered map
void ConvertKeys::keys_to_globals(InputReader& key_map) {
    bool fullyInt = false;
    vector<string> names = key_map.getNames();
    for (string name : names){
        if (name == "NUMBER OF ORBITALS"){
            int num_orb = key_map.getIntKeyword(name);
            global_ints["num_orbitals"] = num_orb;
            global_ints["num_orbitalsC2"] = num_orb * (num_orb + 1) / 2;
        }
        else if (name == "ENERGY TOLERANCE")
            global_doubles["energy_tol"] = key_map.getDoubleKeyword(name);
        else if (name == "RESIDUAL TOLERANCE")
            global_doubles["norm_tol"] = key_map.getDoubleKeyword(name);
        else if (name == "VALENCE CI FLAG")
            global_ints["valence_ci_flag"] = key_map.getIntKeyword(name);
        else if (name == "RESTART FLAG"){
            global_bools["restart_flag"] = key_map.getBoolKeyword(name);
            if (global_bools["restart_flag"])
                throw runtime_error("The restart flag is currently not implemented!");
        }
        else if (name == "PEDERSEN CD")
            global_bools["PEDERSEN CD"] = key_map.getBoolKeyword(name);
        else if (name == "INTEGRAL DIRECT MODE"){
            if (!fullyInt)
                global_bools["integralDirect"] = key_map.getBoolKeyword(name);
        }
        else if (name == "FULLY INTEGRAL DIRECT MODE"){
            bool flag = key_map.getBoolKeyword(name);
            global_bools["fullyIntegralDirect"] = flag;
            if (flag){
                global_bools["integralDirect"] = true;
                global_bools["cdVecsInMemory"] = true;
                fullyInt = true;
            }
        }
        else if (name == "NO DIRECT FOUR INTERNAL")
            global_bools["directFourInternal"] = !key_map.getBoolKeyword(name);
        else if (name == "LOW MEMORY DIRECT")
            global_bools["directLowMem"] = key_map.getBoolKeyword(name);
        else if (name == "SUPER LOW MEMORY DIRECT")
            global_bools["directSuperLowMem"] = key_map.getBoolKeyword(name);
        else if (name == "CD VECTORS IN MEMORY"){
            if (!fullyInt)
                global_bools["cdVecsInMemory"]= key_map.getBoolKeyword(name);
        }
        else if (name == "DENSITY FITTING")
            global_bools["DENSITY FITTING"] = key_map.getBoolKeyword(name);
        else if (name == "CARTESIAN")
            global_bools["spherical"] = !key_map.getBoolKeyword(name);
        else if (name == "NUMBER OF BASIS FUNCTIONS"){
            int num_bas = key_map.getIntKeyword(name);
            global_ints["number_bas"] = num_bas;
            global_ints["number_basC2"] = num_bas * (num_bas + 1) / 2;
        }
        else if (name == "NUMBER OF FROZEN ORBITALS")
            global_ints["num_frozen"] = key_map.getIntKeyword(name);
        else if (name == "NUMBER OF INACTIVE ORBITALS")
            global_ints["num_inactive"] = key_map.getIntKeyword(name);
        else if (name == "NUMBER OF ACTIVE ORBITALS")
            global_ints["num_active"] = key_map.getIntKeyword(name);
        else if (name == "SPIN MULTIPLICITY")
            global_ints["spinM"] = key_map.getIntKeyword(name);
        else if (name == "NUMBER OF ELECTRONS")
            global_ints["num_elec"] = key_map.getIntKeyword(name);
        else if (name == "NUMBER OF REFERENCES")
            global_ints["num_ref"] = key_map.getIntKeyword(name);
        else if (name == "NATURAL ORBITAL FLAG")
            global_ints["nat_orb_flag"] = key_map.getIntKeyword(name);
        else if (name == "SCRATCH DIRECTORY")
            global_strings["scratch_directory"] = key_map.getStringKeyword(name) + '/';
        else if (name == "INTEGRAL THRESHOLD")
            global_doubles["integral_threshold"] = key_map.getDoubleKeyword(name);
        else if (name == "AO INTEGRAL THRESHOLD")
            global_doubles["ao_integral_threshold"] = key_map.getDoubleKeyword(name);
        else if (name == "SKIP INACTIVE PM")
        	global_bools["localize_inactives"] = !(key_map.getBoolKeyword(name));
        else if (name == "SKIP LOVOS")
        	global_bools["skip_lovos"] = key_map.getBoolKeyword(name);
        else if (name == "OCCUPATION THRESHOLD")
            global_doubles["internal_threshold"] = key_map.getDoubleKeyword(name);
        else if (name == "WP DEFAULT RADIUS")
            global_doubles["wp_default_radius"] = key_map.getDoubleKeyword(name);
        else if (name == "TOV OCCUPIED DEFAULT RADIUS")
            global_doubles["tov_occupied_default_radius"] = key_map.getDoubleKeyword(name);
        else if (name == "TOV OCCUPIED RADIUS MULTIPLIER")
            global_doubles["tov_occupied_multiplier"] = key_map.getDoubleKeyword(name);
        else if (name == "WP RADIUS MULTIPLIER")
            global_doubles["wp_multiplier"] = key_map.getDoubleKeyword(name);
        else if (name == "VIRTUAL OCCUPATION THRESHOLD")
            global_doubles["virtual_threshold"] = key_map.getDoubleKeyword(name);
        else if (name == "TOV VIRTUAL DEFAULT RADIUS")
            global_doubles["tov_virtual_default_radius"] = key_map.getDoubleKeyword(name);
        else if (name == "SPHERE PRINT")
            global_bools["sphereprint"] = key_map.getBoolKeyword(name);
        else if (name == "TOV VIRTUAL RADIUS MULTIPLIER")
            global_doubles["tov_virtual_multiplier"] = key_map.getDoubleKeyword(name);
        else if (name == "TOV CYLINDER RADIUS")
            global_doubles["tov_cylinder_radius"] = key_map.getDoubleKeyword(name);
        else if (name == "REFERENCE CI FLAG")
            global_ints["reference_ci_flag"] = key_map.getIntKeyword(name);
        else if (name == "CYLINDER RADIUS")
            global_doubles["wp_cylinder_radius"] = key_map.getDoubleKeyword(name);
        else if (name == "ACPF FLAG")
            global_ints["acpf_flag"] = key_map.getIntKeyword(name);
        else if (name == "CD THRESHOLD")
            global_doubles["cd_thresh"] = key_map.getDoubleKeyword(name);
        else if (name == "CD MEMORY")
            global_ints["max_mem"] =  key_map.getIntKeyword(name);
        else if (name == "INT MEMORY")
            global_ints["max_mem_ints"] = key_map.getIntKeyword(name);
        else if (name == "AIJK MEMORY")
            global_ints["for_buf_maxmemAIJK"] = key_map.getIntKeyword(name);
        else if (name == "NONLOCAL")
            global_ints["nonlocal_flag"] = key_map.getIntKeyword(name);
        else if (name == "NUM ROOTS")
            global_ints["num_roots"] = key_map.getIntKeyword(name);
        else if (name == "CD RESTART")
            global_bools["restart_from_old_CD"] = key_map.getBoolKeyword(name);
        else if (name == "DAVIDSON HARD DRIVE SPACE")
            global_doubles["ASSUMED_DAV_DISK_SPACE"] = key_map.getDoubleKeyword(name);
        else if (name == "NO ACPF ROOT FOLLOW")
            global_bools["acpf_root_follow"] = !(key_map.getBoolKeyword(name));
        else if (name == "ROOT FOLLOW HELP")
            global_bools["ACPF_ROOT_FOLLOW_HELP"] = key_map.getBoolKeyword(name);
        else if (name == "CUSTOM G VAL")
            global_doubles["custom_g_val"] = key_map.getDoubleKeyword(name);
        else if (name == "NUM THREADS")
            global_ints["numThreads"] = key_map.getIntKeyword(name);
        else if (name == "STRIDE SIZE SIGMA")
            global_ints["blockSizeSigma"] = key_map.getIntKeyword(name);
        else if (name == "STRIDE SIZE CI")
            global_ints["blockSizeCI"] = key_map.getIntKeyword(name);
        else if (name == "IOBUFFBLOCKINTS")
            global_ints["for_buf_blocksizeInts"] = key_map.getIntKeyword(name);
        else if (name == "IOBUFFMAXMEMINTS")
            global_ints["for_buf_maxmemInts"] = key_map.getIntKeyword(name);
        else if (name == "DONTSTOREINTS")
            global_bools["for_buf_storeIntegrals"] = !(key_map.getBoolKeyword(name));
        else if (name == "IOBUFFMAXMEMCD")
            global_ints["for_buf_maxmemCD"] = key_map.getIntKeyword(name);
        else if (name == "DONTSTORECDVECS")
            global_bools["for_buf_storeCDVecs"] = !(key_map.getBoolKeyword(name));
        else if (name == "IOBUFFMAXMEMSEGS")
            global_ints["for_buf_maxmemSegs"] = key_map.getIntKeyword(name);
        else if (name == "IOBUFFBLOCKSEGS")
            global_ints["for_buf_blocksizeSegs"] = key_map.getIntKeyword(name);
        else if (name == "REFERENCE OCCUPATIONS"){
            global_int_vecs["references"] = key_map.getRefVecsKeyword(name);
        }
        else if (name == "CYLINDERS"){
            vector<Chain> cyl_vec = key_map.getCylVecsKeyword(name);
            if (cyl_vec.size() != 0){
                for (Chain chain : cyl_vec){
                    pair<int, Chain> p (chain.getOrb(), chain);
                    cyl_endpoints.insert(p);
                }
                cout << "manual cylinders:" << endl;
                for (pair<int, Chain> p : cyl_endpoints)
                    p.second.print();
            }
            else
                global_bools["use_cyl_endpoints"] = false;
        }
        else if (name == "XYZ FILE"){
            std::cerr << "IGNORING XYZ FILE SETTING " << key_map.getStringKeyword(name) << std::endl;
        }
        else if (name == "ORB FILE"){
            my_strings["orb_file"] = key_map.getStringKeyword(name);
        }
        else if (name == "BASIS SET"){
            my_strings["basis_set"] = key_map.getStringKeyword(name);
        }
        else if (name == "AUXILIARY BASIS SET"){
            my_strings["auxiliary basis set"] = key_map.getStringKeyword(name);
        }
        else if (name == "EMBEDDING"){
            global_bools["embedding"] = key_map.getBoolKeyword(name);
        }
        else if (name == "GRID DIMENSIONS"){
            global_int_vecs["grid dimensions"] = key_map.getIntVecsKeyword(name);
        }
        else if (name == "LATTICE VECTORS"){
            global_double_vecs["lattice vectors"] = key_map.getDoubleVecsKeyword(name);
        }
        else if (name == "SHIFT VECTOR"){
            global_double_vecs["shift vector"] = key_map.getDoubleVecsKeyword(name);
        }
        else if (name == "EMBEDDING POTENTIAL FILE"){
            global_strings["embedding potential file"] = key_map.getStringKeyword(name);
        }

        else
            throw runtime_error("Tried to read variable " + name +
                    " but it cannot be found");
    }
    // take care of a few reference things
    global_ints["num_internal"] = global_ints.at("num_inactive")
            + global_ints.at("num_active");
    if (global_ints["nonlocal_flag"])
      setup_nonlocal();
    fill_references();
    check_input();
}

// get int val for global variable c_name
void c_get_global_int(const char *c_name, int64_t& val){
    string name(c_name);
    val = global_ints.at(name);
}

// get double val for global variable c_name
void c_get_global_double(const char *c_name, double& val) {
    string name(c_name);
    val = global_doubles.at(name);
}

// get bool val for global variable c_name
void c_get_global_bool(const char *c_name, int64_t& val) {
    string name(c_name);
    bool bool_val = global_bools.at(name);
    if (bool_val)
        val = 1;
    else
        val = 0;
}

double ConvertKeys::get_global_double(string name) const{
    return global_doubles.at(name);
}

int ConvertKeys::get_global_int(string name) const{
    return global_ints.at(name);
}

bool ConvertKeys::get_global_bool(string name) const{
    return global_bools.at(name);
}

string ConvertKeys::get_global_string(string name) const{
    return global_strings.at(name);
}

vector<vector<int>> ConvertKeys::get_global_int_vecs(string name) const{
    return global_int_vecs.at(name);
}

vector<vector<double>> ConvertKeys::get_global_double_vecs(string name) const{
    return global_double_vecs.at(name);
}

void ConvertKeys::add_global_int(string name, int val){
    global_ints[name] = val;
}

void ConvertKeys::add_global_bool(string name, bool val){
    global_bools[name]=val;
}

// get string val for global variable c_name
void c_get_global_string(const char *c_name, char *val, int64_t& length){
    string name(c_name);
    string s = global_strings.at(name);
    strcpy(val, s.c_str());
    length = s.length();
}

// get int array val for vector num of global variable c_name
void c_get_global_int_vec(const char *c_name, int64_t& num, int64_t *val) {
    string name(c_name);
    vector<vector<int>> key = global_int_vecs.at(name);
    vector<int> vec = key[num];
    
    int i = 0;
    
    for (int64_t int_val : vec)
        val[i++] = int_val;
}

// get cylinders on a given orbital
void c_get_cylinders(int64_t& orb, int64_t *val){
    Chain chain = cyl_endpoints.at(orb);
    vector<Cylinder> cyls = chain.getCylinders();
    
    int i = 0;
    for (Cylinder cyl : cyls){
        val[i++] = cyl.firstEndpoint();
        val[i++] = cyl.secondEndpoint();
    }
}

// get number of cylinders for a given orbital
void c_num_cylinders(int64_t& orb, int64_t& num_cyl){
    Chain chain;
    try{
        chain = cyl_endpoints.at(orb);
    }
    catch(const out_of_range& oor){
        num_cyl = 0;
        return;
    }
    num_cyl = chain.numCylinders();
}

// get string val for global variable c_name
string ConvertKeys::get_my_string(string name) const{
    return my_strings.at(name);
}

// fill all references with zeros
void ConvertKeys::fill_references(){
    size_t num_internal = global_ints.at("num_internal");
    size_t num_orbitals = global_ints.at("num_orbitals");
    vector<vector<int>>& refs = global_int_vecs.at("references");

    for (vector<int>& ref : refs){
      if (num_internal != ref.size()){
	throw runtime_error("The length of each reference must be "
          + string("equal to the sum of active and inactive orbitals."));
      }
      for (size_t i = num_internal; i < num_orbitals; i++)
            ref.push_back(0);
    }
}

// do some input parameter checking
void ConvertKeys::check_input(){

  // check number of threads
  int numThreads = global_ints.at("numThreads");
#if !defined(_OPENMP)
  if (numThreads > 1)
    throw runtime_error("ERROR: This is a serial code and you want more than one thread!");
#endif

  // some ACPF checks
  int acpf_flag = global_ints.at("acpf_flag");
  int reference_ci_flag = global_ints.at("reference_ci_flag");
  bool restart_flag = global_bools.at("restart_flag");
  bool acpf_root_follow = global_bools.at("acpf_root_follow");
  int num_roots = global_ints.at("num_roots");

  if (acpf_flag && reference_ci_flag){
    stringstream err("WHAT? You want me to run a acpf calculation\n");
    err << "with only references? That doesn't make any\n";
    err << "sense. Please fix the input file!\n";
    throw runtime_error(err.str());
  }

  if (acpf_flag && restart_flag){
    stringstream err("You have asked to restart an ACPF calculation.\n");
    err << "Honestly not sure if this will work (not sure if the reference energy" ;
    err << " will be calculated correctly \n";
    throw runtime_error(err.str());
  }
  
  if (acpf_flag && acpf_root_follow && (num_roots > 1)){
    cout << endl;
    cout << "*******************************************************************" << endl;
    cout << "* You want me to run an ACPF calculation with root following      *" << endl;
    cout << "* on more than 1 root. Be careful, this hasn't really been tested.*" << endl;
    cout << "*******************************************************************" << endl;
    cout << endl;
  }

  // check AO threshold
  double ao_integral_threshold = global_doubles.at("ao_integral_threshold");
  double cd_thresh = global_doubles.at("cd_thresh");

  if (ao_integral_threshold > cd_thresh){
    stringstream err("TigerCI does not have stable performance when the AO Integral");
    err << "threshold is looser than the Cholesky decomposition threshold.";
    err << "Please change your input parameters.";
    throw runtime_error(err.str());
  }

  // check the references
  vector<vector<int>>& refs = global_int_vecs.at("references");
  unique_ref_check(refs);

  // check for overly small systems
  int num_internal = global_ints.at("num_internal");
  int num_orbitals = global_ints.at("num_orbitals");
  int num_external = num_orbitals - num_internal;
  int num_elec = global_ints.at("num_elec");

  if (num_internal < 2)
    throw runtime_error("Sorry, this code won't handle calculations with fewer than 2 internal orbitals");
  if (num_external < 4)
    throw runtime_error("Sorry, this code won't handle calculations with fewer than 4 external orbitals");
  if (num_elec <= 2)
    throw runtime_error("Sorry, this code won't handle calculations with fewer than 3 electrons");
}

//> \brief Checks to make sure that each reference is unique
//> \param num_ref    The number of references
//> \param references The occupation(internal) of each reference
//> 
//> You wouldn't think this is a big deal. Thanks to the graphical methods a calculation with 2 identical references
//> isn't more expensive. But it does break a fundamental assumption we use that each path we get out of the treesearch
//> represents a unique internal occupation (somehow)
void ConvertKeys::unique_ref_check(vector<vector<int>>& refs){
  int num_refs = refs.size();
  for (int i = 0; i < num_refs; i++){
    for (int j = i+1; j < num_refs; j++){
      vector<int> refA = refs.at(i);
      vector<int> refB = refs.at(j);
      if (equal(refA.begin(), refA.end(), refB.begin())){
	stringstream err;
	err << "FATAL INPUT ERROR: References " << i << " and " << j << " are identical!";
	throw runtime_error(err.str());
      }
    }
  }
}

// Here we reset parameters if the nonlocal_flag is active
void ConvertKeys::setup_nonlocal(){
  cout << endl;
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  cout << "WARNING NONLOCAL FLAG DETECTED .... " << endl;
  cout << "THE FOLLOWING PARAMETERS ARE BEING  " << endl;
  cout << "READJUSTED TO PERFORM A NONLOCAL    " << endl;
  cout << "CALCULATION. USER INPUTS MAY BE     " << endl;
  cout << "OVERWRITTEN                         " << endl;
  cout << endl;
  cout << "INTEGRAL THRESHOLD SET TO                    0.0" << endl;
  cout << "OCCUPATION THRESHOLD SET TO                 .999" << endl;
  cout << "PAO OCCUPATION THRESHOLD SET TO             .999" << endl;
  cout << "DEFAULT RADIUS SET TO                10000000.D0" << endl;
  cout << "DEFAULT RADIUS MULTIPLIER SET TO     10000000.D0" << endl;
  cout << "MO DEFAULT RADIUS SET TO             10000000.D0" << endl;
  cout << "MO DEFAULT RADIUS MULTIPLIER SET TO  10000000.D0" << endl;
  cout << "TOV DEFAULT RADIUS SET TO            10000000.D0" << endl;
  cout << "TOV RADIUS MULTIPLIER SET TO         10000000.D0" << endl;
  cout << "PAO RADIUS SET TO                    10000000.D0" << endl;
  cout << "PAO RADIUS MULTIPLIER SET TO         10000000.D0" << endl;
  cout << "CYLINDER RADIUS SET TO               10000000.D0" << endl;
  cout << "TOV CYLINDER RADIUS SET TO           10000000.D0" << endl;
  cout << "LargeCholeskyVecDistanceCutOff is true          " << endl;
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
  cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;

  global_doubles.at("integral_threshold") = 0.0;
  global_doubles.at("internal_threshold") = 0.999E0;
  global_doubles.at("virtual_threshold")  = 0.999E0;
  global_doubles.at("wp_default_radius")  = 10000000.E0;
  global_doubles.at("tov_occupied_default_radius") = 10000000.E0;
  global_doubles.at("tov_occupied_multiplier") = 10000000.E0;
  global_doubles.at("wp_multiplier")  = 10000000.E0;
  global_doubles.at("tov_virtual_default_radius") = 10000000.E0;
  global_doubles.at("tov_virtual_multiplier") = 10000000.E0;
  global_doubles.at("tov_cylinder_radius")   = 10000000.E0;
  global_doubles.at("wp_cylinder_radius")    = 10000000.E0;
  global_doubles.at("ao_integral_threshold") = 1.0E-15;
}
