/*
 * @file EmbeddingData.h
 *
 *      Author: Caroline M. Krauter
 */

#ifndef TIGER_CPP_EMBEDDING_H_
#define TIGER_CPP_EMBEDDING_H_

#include <tiger_cpp/ConvertKeys.h>
#include <tiger_cpp/BasisData.h>

#include <armadillo>
#include <vector>

/**
 * @class EmbeddingData
 * @brief Contains the implementation of embedding data 
 */

class EmbeddingData {
public:
	// constructor
	EmbeddingData(ConvertKeys& global_vars, BasisData *bas);
	~EmbeddingData();

	// Compute the one-electron integrals over the embedding potential
	arma::mat compute_emb_one_ints();

private:
	//generate embedding data
    double compute_cell_volume(std::vector<std::vector<double>>& lattice_vecs);
    int compute_nr_grid_pts(std::vector<int>& grid_dimension);
    std::vector<double> read_embedding_potential(int nr_grid_points);
	std::vector<std::vector<double>> convert_to_bohr(std::vector<std::vector<double>>& vecs);

    //data
    BasisData *_basis;
	ConvertKeys *_globals;

	// grid information
    std::vector<int> _grid_dimension;
    std::vector<std::vector<double>> _shift_vec;
    std::vector<std::vector<double>> _lattice_vecs;
	int _nr_grid_points;
	double _cell_volume;

	// embedding potential
	std::vector<double> _embedding_potential;

	// Use this conversion factor for the moment for better comparability with Molcas
	//
	// from embedding implementation in Molcas. Originally from (checked):
	// *     Conversion factor angstrom to bohr from the IUPAC
	// *     publication
	// *     .529177249(24) angstrom / bohr
	// *     "Quantities, Units and Symbols in Physical Chemistry"
	// *     I. Mills, T. Cvitas, K. Homann, N. Kallay and
	// *     K. Kuchitsu, Blackwell Scientific Publications,
	// *     Oxford, 1988.
	const double _bohr_to_angstr = 1./0.529177249;

};

#endif /* TIGER_CPP_EMBEDDING_H_ */
