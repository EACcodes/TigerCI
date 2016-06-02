/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * @file BasisData.h
 *
 *  Created on: Jan 15, 2015
 *      Author: Francis Ricci
 */

#ifndef TIGER_CPP_BASIS_H_
#define TIGER_CPP_BASIS_H_

#include <erkale/settings.h>
#include <erkale/basislibrary.h>
#include <erkale/basis.h>

#include <tiger_cpp/ConvertKeys.h>
#include <tiger_cpp/ErkaleData.h>
#include <tiger_cpp/ErkaleIntegral.h>
#include <tiger_cpp/ErkaleDF.h>
#include <tiger_cpp/PairIndex.h>
#include <tiger_cpp/TigerStructs.h>

#include <armadillo>
#include <vector>
#include <memory>

/**
 * @class BasisData
 * @brief An abstract class, presenting an interface for basis set and integral data
 */

class BasisData {
public:
	//copy constructor
	virtual BasisData* duplicate() = 0;

	//basis data
	virtual double get_nuclear_repulsion() = 0;
	virtual std::vector<tiger::coords_t> get_coordinates() = 0;
	virtual std::vector<tiger::atom_t> get_atoms() = 0;
	virtual std::vector<int> get_nFuncInShell() = 0;
	virtual int get_nBas() = 0;
	virtual int get_nShell() = 0;
	virtual int get_nAtoms() = 0;
	virtual std::vector<int> get_basis_map() = 0; //maps basis functions to atom number
	virtual arma::vec eval_func(double x, double y, double z) = 0; // evaluates the basis functions at point (x,y,z)

	//integrals
	virtual std::vector<double> eval_ijkl(int i, int j, int k, int l) = 0;
	virtual arma::mat get_one_ints() = 0;
	virtual arma::mat get_overlap() = 0;

	//LOVOs
	virtual arma::mat get_overlap(BasisData *bas) = 0;
	virtual vector<vector<size_t>> separate() = 0; //groups basis funcs by atom and then am
	virtual vector<size_t> funcs() = 0; //groups basis funcs by atom

	//DF
	virtual void compute_three_center(arma::mat& integrals, PairIndex& addr) = 0;
	virtual void compute_two_center(arma::mat& ab) = 0;
	virtual int get_nAux() = 0;

	virtual ~BasisData();

	ConvertKeys *_globals;

protected:
	tiger::basis_type _type;
};

/**
 * @class ErkaleBasis
 * @brief Uses Erkale to implement the BasisData interface for basis set and integral data
 */

class ErkaleBasis : public BasisData{
public:
	ErkaleBasis(ConvertKeys& global_vars, tiger::basis_type type);
	~ErkaleBasis();

	//copy constructor (for thread-safe integrals)
	BasisData* duplicate();
	ErkaleBasis(ErkaleBasis& bas);

	//basis data
	double get_nuclear_repulsion();
	std::vector<tiger::coords_t> get_coordinates();
	std::vector<tiger::atom_t> get_atoms();
	std::vector<int> get_nFuncInShell();
	int get_nBas();
	int get_nShell();
	int get_nAtoms();
	std::vector<int> get_basis_map();
	arma::vec eval_func(double x, double y, double z);

	//integrals
	std::vector<double> eval_ijkl(int i, int j, int k, int l);
	arma::mat get_one_ints();
	arma::mat get_overlap();

	//LOVOs
	arma::mat get_overlap(BasisData *bas);
	vector<vector<size_t>> separate();
	vector<size_t> funcs();

	//DF
	void compute_three_center(arma::mat& integrals, PairIndex& addr);
	void compute_two_center(arma::mat& ab);
	int get_nAux();

protected:
	std::shared_ptr<BasisSet> create_aux();
	std::shared_ptr<BasisSet> create_basis();

	std::shared_ptr<BasisSet> basis;
	std::shared_ptr<BasisSet> auxbasis;

	ErkaleData *data;
	ErkaleIntegral *ints;
	ErkaleDF *df;
};

// external and without C++ "name mangling" to allow for calling by Fortran
extern "C"
{
/* Basis set data */

void init_bas(BasisData *data);

// get nuclear repulsion energy, returning value into nuclear_energy
void c_get_nuclear_repulsion(double& nuclear_energy);

// get number of atoms
void c_get_num_atoms(int64_t& num_atoms);

// get number of shells
void c_nr_shells(int64_t& nshell);

// get nuclear coordinates
void c_get_coordinates(int64_t& num_atoms, double *coordinates);

// data about location and number of basis functions
void c_basis_func_data(int64_t *nFuncInShell, int64_t& number_bas,
		int64_t& nShell);

// generate mapping of each basis function to its atomic center
void c_get_basis_map(int64_t& number_bas, int64_t *basisfunc2atom);

/* Integral information */

// initialize integrals (for calling from the Fortran)
void init_ints(ErkaleIntegral *ints);

void set_cho(int nCho);
void c_get_cholesky_data(int64_t& nCho);

// get two electron integrals for shells ijkl, returning values into TInt
void c_eval_ijkl(int64_t& iShell, int64_t& jShell, int64_t& kShell,
		int64_t& lShell, double *TInt, int64_t& nTInt);

// get one electron integrals, returning values into one_ints_AO
void c_get_one_ints(double *one_ints_AO, int64_t& number_bas, int64_t& nao);

// get overlap integrals, returning values into AOoverlap
void c_get_overlap(int64_t& nBas, double *AOoverlap);
}

#endif /* TIGER_CPP_BASIS_H_ */
