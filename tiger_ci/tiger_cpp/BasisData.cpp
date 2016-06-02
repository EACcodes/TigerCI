/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * @file BasisData.cpp
 *
 *  Created on: Jan 15, 2015
 *  Updated on: Jan 7, 2016
 *      Author: Francis Ricci, Johannes M Dieterich
 */

#include <tiger_cpp/BasisData.h>
#include <tiger_cpp/EmbeddingData.h>
#include <tiger_cpp/OrbitalParsers.h>

// pointer to an BasisData object (stored by main.cpp), which allows the
// Fortran interface functions to access the class objects (a workaround for
// Fortran not having access to C++ objects)
static BasisData *data;
static int64_t nCho;
static bool init = false;

/* C/Fortran interface functions */

void init_bas(BasisData *bas_data){
	if (init)
		return;
	data = bas_data;
	init = true;
}

// get nuclear repulsion energy, returning value into nuclear_energy
void c_get_nuclear_repulsion(double& nuclear_energy){
    nuclear_energy = data->get_nuclear_repulsion();
}

// get number of atoms
void c_get_num_atoms(int64_t& num_atoms){
    num_atoms = (int64_t) data->get_nAtoms();
}

// get number of shells
void c_nr_shells(int64_t& nshell){
    nshell = (int64_t) data->get_nShell();
}

// get nuclear coordinates
void c_get_coordinates(int64_t& num_atoms, double *coordinates){
	if (num_atoms != (int64_t) data->get_nAtoms())
		throw runtime_error("Number of atoms does not match basis set");

	vector<tiger::coords_t> coords = data->get_coordinates();

	int idx = 0;
	for (int i = 0; i < num_atoms; i++){
		tiger::coords_t coord = coords[i];
		coordinates[idx++] = coord.x;
		coordinates[idx++] = coord.y;
		coordinates[idx++] = coord.z;
	}
}

// data about location and number of basis functions
void c_basis_func_data(int64_t *nFuncInShell, int64_t& number_bas, int64_t& nShell){
    if (number_bas != (int64_t) data->get_nBas())
        throw runtime_error("Mismatch in number of basis functions");

    nShell = (int64_t) data->get_nShell();
    vector<int> nFunc = data->get_nFuncInShell();

    for (int i = 0; i < nShell; i++){
        nFuncInShell[i] = (int64_t) nFunc.at(i);
    }
}

// generate mapping of each basis function to its atomic center
void c_get_basis_map(int64_t& number_bas, int64_t *basismap){
    if (number_bas != (int64_t) data->get_nBas())
        throw runtime_error("Mismatch in number of basis functions");

    vector<int> map = data->get_basis_map();
    for (int i = 0; i < number_bas; i++)
        basismap[i] = (int64_t) map[i];
}

/* Integral stuff */

void set_cho(int n){
	nCho = (int64_t) n;
}

void c_get_cholesky_data(int64_t& n){
	n = nCho;
}

// get two electron integrals for shells ijkl, returning values into TInt
void c_eval_ijkl(int64_t& i, int64_t& j, int64_t& k, int64_t& l,
		double *TInt, int64_t& nTInt){

	vector<double> eris = data->eval_ijkl(i, j, k, l);

	int64_t num_eris = eris.size();

	if (nTInt < num_eris)
		throw length_error("Insufficient memory allocated for two-electron integrals");

	//convert vector to array
	for (int m = 0; m < num_eris; m++){
		TInt[m] = eris[m];
	}
}

// get one electron integrals, returning values into one_ints_AO
void c_get_one_ints(double *one_ints_AO, int64_t& number_bas, int64_t& nao){
	//calculate one electron integrals
	arma::mat energy_ints = data->get_one_ints();

	if ((energy_ints.n_rows != number_bas) || (energy_ints.n_cols != number_bas))
		throw runtime_error("Number of one-electron integrals does not match "
				+ string("number of basis functions"));

	// convert upper-triangular matrix of integrals to one-dimensional array
	int64_t ao_count = 0;
	for (int64_t i = 0; i < number_bas; i++){
		for (int64_t j = 0; j <= i; j++){
			one_ints_AO[ao_count++] = energy_ints(j,i);

			if (ao_count > nao){
				throw length_error("Insufficient memory allocated for " +
						string("one-electron integrals"));
			}
		}
	}
}

// get overlap integrals, returning values into AOoverlap
void c_get_overlap(int64_t& nBas, double *AOoverlap){
	//calculate overlap integrals
	arma::mat overlap_ints = data->get_overlap();

	if ((overlap_ints.n_rows != nBas) || (overlap_ints.n_cols != nBas))
		throw runtime_error("Number of one-electron integrals does not match"
				+ string(" number of basis functions"));

	double *ptr = AOoverlap;

	//convert integral matrix to two-dimensional array
	for (int64_t i = 0; i < nBas; i++){
		for (int64_t j = 0; j < nBas; j++){
			*ptr = overlap_ints(j,i); //[j][i] for Fortran arrays
			ptr++;
		}
	}
}

BasisData::~BasisData(){}

ErkaleBasis::ErkaleBasis(ConvertKeys& global_vars, tiger::basis_type type) {
	init_libint_base();

    _globals = &global_vars;
    _type = type;

    basis = create_basis();

    data = new ErkaleData(basis);
    ints = new ErkaleIntegral(basis);
	if (_globals->get_global_bool("DENSITY FITTING")){
		auxbasis = create_aux();
		df = new ErkaleDF(basis, auxbasis);
	}

    init_bas(this);
}

ErkaleBasis::~ErkaleBasis(){
	delete data;
	delete ints;
	if (_globals->get_global_bool("DENSITY FITTING")){
		delete df;
	}
}

std::shared_ptr<BasisSet>  ErkaleBasis::create_aux(){
	tiger::basis_type old_type = _type;
	_type = tiger::aux;
	std::shared_ptr<BasisSet> aux = create_basis();
	_type = old_type;
	return aux;
}

std::shared_ptr<BasisSet> ErkaleBasis::create_basis(){
    Settings erkale_set;
    BasisSetLibrary basis_lib;
    std::shared_ptr<BasisSet> erk_basis(new BasisSet());

    bool spherical = _globals->get_global_bool("spherical");

    erkale_set.add_bool("UseLM", "", spherical); // don't use spherical harmonics
    erkale_set.add_string("Decontract", "", "");
    erkale_set.add_bool("BasisRotate", "", false);
    erkale_set.add_double("BasisCutoff", "", 1.0e-9);

    // parse the nuclear positions using the molden file
    MoldenOrbitalParser parser = MoldenOrbitalParser(spherical);
    std::vector<atom_t> atoms = parser.readCoordinates(_globals->get_my_string("orb_file"));
    
    string basis_set;
    switch (_type){
        case tiger::full: basis_set = _globals->get_my_string("basis_set");
            break;
        case tiger::minimal: basis_set = "STO-3G";
            break;
        case tiger::aux: basis_set = _globals->get_my_string("auxiliary basis set");
            break;
        default: throw runtime_error("Invalid basis type");
    }

    basis_lib.load_gaussian94(basis_set, false);

    construct_basis(*erk_basis, atoms, basis_lib, erkale_set);
    erk_basis->finalize(false,false);
    return erk_basis;
}

BasisData* ErkaleBasis::duplicate(){
	BasisData *new_bas = new ErkaleBasis(*this);
	return new_bas;
}

ErkaleBasis::ErkaleBasis(ErkaleBasis& bas){
	_globals = bas._globals;
	_type = bas._type;
    basis = bas.basis;

    data = new ErkaleData(basis);
    ints = new ErkaleIntegral(basis);
    if (_globals->get_global_bool("DENSITY FITTING")){
    	auxbasis = create_aux();
    	df = new ErkaleDF(basis, auxbasis);
    }
}

double ErkaleBasis::get_nuclear_repulsion(){
	return data->get_nuclear_repulsion();
}

std::vector<tiger::coords_t> ErkaleBasis::get_coordinates(){
	return data->get_coordinates();
}

std::vector<tiger::atom_t> ErkaleBasis::get_atoms(){
	return data->get_atoms();
}

std::vector<int> ErkaleBasis::get_nFuncInShell(){
	return data->get_nFuncInShell();
}

int ErkaleBasis::get_nBas(){
	return data->get_nBas();
}

int ErkaleBasis::get_nShell(){
	return data->get_nShell();
}

int ErkaleBasis::get_nAtoms(){
	return data->get_nAtoms();
}

std::vector<int> ErkaleBasis::get_basis_map(){
	return data->get_basis_map();
}

arma::vec ErkaleBasis::eval_func(double x, double y, double z){
    return data->eval_func(x, y, z);
}

std::vector<double> ErkaleBasis::eval_ijkl(int i, int j, int k, int l){
	return ints->eval_ijkl(i, j, k, l);
}

arma::mat ErkaleBasis::get_one_ints(){
    arma::mat one_ints = ints->get_one_ints();
    if (_globals->get_global_bool("embedding")) {
		std::cout << "***********************" << std::endl;
		std::cout << "*                      " << std::endl;
		std::cout << "* Embedding Integrals  " << std::endl;
		std::cout << "*                      " << std::endl;
		std::cout << "***********************" << std::endl;
		std::cout << endl;

        EmbeddingData *embpot = new EmbeddingData(*_globals, duplicate());
		one_ints += embpot->compute_emb_one_ints();
	}
	return one_ints;
}

arma::mat ErkaleBasis::get_overlap(){
	return ints->get_overlap();
}

/* LOVO stuff */

arma::mat ErkaleBasis::get_overlap(BasisData *bas){
	ErkaleBasis *erk_bas = dynamic_cast<ErkaleBasis*>(bas);
	return ints->get_overlap(erk_bas->ints);
}

// separate shells by angular momentum
vector<vector<size_t>> ErkaleBasis::separate(){
    return data->separate();
}

vector<size_t> ErkaleBasis::funcs(){
	return data->funcs();
}

/* DF stuff */

void ErkaleBasis::compute_three_center(arma::mat& integrals, PairIndex& addr){
	return df->compute_three_center_integrals(integrals, addr);
}

void ErkaleBasis::compute_two_center(arma::mat& ab){
	return df->compute_two_center_integrals(ab);
}

int ErkaleBasis::get_nAux(){
	return df->get_nAux();
}
