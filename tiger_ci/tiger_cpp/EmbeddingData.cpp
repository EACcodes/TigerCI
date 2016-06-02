/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * @file EmbeddingData.cpp
 *
 *      Author: Caroline M. Krauter
 */

#include <tiger_cpp/EmbeddingData.h>
#include <tiger_cpp/Timer.h>

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <fstream>
#include <iostream>


EmbeddingData::EmbeddingData(ConvertKeys& global_vars, BasisData *bas) {
	_globals = &global_vars;
	_basis = bas;

	Timer embedding_info;

	// get information about cell from input file
	std::vector<std::vector<int>> vecs = _globals->get_global_int_vecs("grid dimensions");
	_grid_dimension = vecs[0];	
	_shift_vec = _globals->get_global_double_vecs("shift vector");
	_lattice_vecs = _globals->get_global_double_vecs("lattice vectors");

	// convert to bohr
	_shift_vec = convert_to_bohr(_shift_vec);
	_lattice_vecs = convert_to_bohr(_lattice_vecs);

	// calculate number of grid points and cell volume
	_nr_grid_points = compute_nr_grid_pts(_grid_dimension);
	_cell_volume = compute_cell_volume(_lattice_vecs);

	// read embedding potential
	_embedding_potential = read_embedding_potential(_nr_grid_points);
	// Debug output
	std::cout << "Size of embpot: " << _embedding_potential.size() << "\nNr grid points: " << _nr_grid_points << "\nCapacity: " << _embedding_potential.capacity() << "\n";

	embedding_info.print_elapsed_time("  Read grid data and embedding potential");
}

arma::mat EmbeddingData::compute_emb_one_ints() {

	Timer embedding_integrals;

	// determine incremental lattice vectors in bohr
	std::vector<std::vector<double>> incr_lattice_vecs (3, std::vector<double>(3));
	for (size_t i=0; i<3; i++){
        for (size_t j =0; j<3; j++){
            incr_lattice_vecs[i][j] = _lattice_vecs[i][j]/_grid_dimension[i];
        }
    }

	// determine volume element for integration
	double volume_element = _cell_volume / _nr_grid_points;

	// Debug output
	std::cout << "nr_grid_points " << _nr_grid_points << std::endl;
	std::cout << std::setprecision (15) << "cell_volume " << _cell_volume << std::endl;
	std::cout << std::setprecision (15) << "volume_elemet " << volume_element << std::endl;
	std::cout << "incremental lattice vectors in bohr: " << std::endl;
	for (size_t i=0; i<incr_lattice_vecs.size();i++){
        for (size_t j =0; j<incr_lattice_vecs[i].size();j++){
            std::cout << incr_lattice_vecs[i][j] << " ";
        } 
		std::cout << std::endl;
    }

	// Initialize matrix for one electron integrals
	int NBas =  _basis->get_nBas(); // Number of basis functions
	arma::mat ints(NBas,NBas);
	ints.zeros();

	// loop over grid points
	size_t i = 0; // index to determine the index of a grid point in V_emb
	for(size_t i_z=0; i_z<_grid_dimension[2]; i_z++) {
		for (size_t i_y=0; i_y<_grid_dimension[1]; i_y++) {
			for (size_t i_x=0; i_x<_grid_dimension[0]; i_x++, i++) {

				// determine coordinates of grid point
				double x = (_shift_vec[0][0] + i_x * incr_lattice_vecs[0][0] + i_y * incr_lattice_vecs[1][0] + i_z * incr_lattice_vecs[2][0]);
				double y = (_shift_vec[0][1] + i_x * incr_lattice_vecs[0][1] + i_y * incr_lattice_vecs[1][1] + i_z * incr_lattice_vecs[2][1]);
				double z = (_shift_vec[0][2] + i_x * incr_lattice_vecs[0][2] + i_y * incr_lattice_vecs[1][2] + i_z * incr_lattice_vecs[2][2]);

				// evaluate basis functions at this point
				arma::vec funcs = _basis->eval_func(x, y, z);

				// add products (\phi_a * V_emb * \phi_b * volume_element) to the matrix of integrals
				// use symmetry of integrals: int(a,b)=int(b,a)
				for (size_t a=0; a<NBas; a++) {
					ints(a,a) += funcs(a) * funcs(a) * volume_element * _embedding_potential[i];
					for (size_t b=0; b<a; b++) {
						ints(a,b) += funcs(a) * funcs(b) * volume_element * _embedding_potential[i];
						ints(b,a) = ints(a,b);	
					}
				}
			}
		}
	}
	
	// print result for debugging
	std::cout << std::setprecision (15) << "One-electron Integrals: \n" << ints << std::endl;

	embedding_integrals.print_elapsed_time("  Calculate embedding integrals");

	return ints;
}


int EmbeddingData::compute_nr_grid_pts(std::vector<int>& grid_dimension) {
	return grid_dimension[0]*grid_dimension[1]*grid_dimension[2];
}

double EmbeddingData::compute_cell_volume(std::vector<std::vector<double>>& lattice_vecs) {

	// volume of cell spanned by lattice vectors a,b,c is:
	// |a * (b x c)|

	// calculate cross product beween b and c
	std::vector<double> cross_prod;
	cross_prod.push_back(lattice_vecs[1][1] * lattice_vecs[2][2] - lattice_vecs[1][2] * lattice_vecs[2][1]);
	cross_prod.push_back(lattice_vecs[1][2] * lattice_vecs[2][0] - lattice_vecs[1][0] * lattice_vecs[2][2]);
	cross_prod.push_back(lattice_vecs[1][0] * lattice_vecs[2][1] - lattice_vecs[1][1] * lattice_vecs[2][0]);

	// Calculate volume 
	return std::abs(cross_prod[0] * lattice_vecs[0][0] + cross_prod[1] * lattice_vecs[0][1] + cross_prod[2] * lattice_vecs[0][2]);
}

std::vector<double> EmbeddingData::read_embedding_potential(int nr_grid_points) {

	std::string filename = _globals->get_global_string("embedding potential file");
	std::ifstream in(filename.c_str());

	// read embedding potential if file successfully opened
	if(in.good()) {
		// reserve enough memory
		std::vector<double> embedding_potential;
		embedding_potential.reserve(nr_grid_points);

		double number;
		while (in >> number) {
			embedding_potential.push_back(number);
		}
		// check whether the correct number of values was read
		if (embedding_potential.size() != nr_grid_points) {
			std::ostringstream oss;
			oss << "Problem reading embedding potential: Number of values read is not the same as number of grid points!\n";
			throw std::runtime_error(oss.str());
		}
		return embedding_potential;
	}
	// throw error if reading not successful
	else {
		std::ostringstream oss;
		oss << "Could not open embedding file \""<<filename<<"\"!\n";
		throw std::runtime_error(oss.str());
	}
}

std::vector<std::vector<double>> EmbeddingData::convert_to_bohr(std::vector<std::vector<double>>& vecs) {
	std::vector<std::vector<double>> vecs_conv = vecs;
	for (size_t i=0; i<vecs.size(); i++){
		for (size_t j =0; j<vecs[i].size(); j++){
			vecs_conv[i][j] = vecs[i][j] * _bohr_to_angstr;
		}	
	}
	return vecs_conv;
}                                                                	
