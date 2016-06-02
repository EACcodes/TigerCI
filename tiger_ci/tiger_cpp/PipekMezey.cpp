/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * File:   pipek_mezey.cpp
 * Author: Francis Ricci
 *
 * Created on August 25, 2014, 11:41 AM
 */

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <iostream>

#include "PipekMezey.h"

PipekMezey::PipekMezey(std::vector<size_t>& bf_per_atom, arma::mat& mo_coeff,
		arma::mat& overlap,	size_t norbitals, size_t nfrozen):
				_bf_per_atom(bf_per_atom),
				_mo_coeff(mo_coeff),
				_overlap(overlap),
				_norbitals(norbitals-nfrozen),
				_natoms(bf_per_atom.size()){}

// Performs Pipek_Mezey localization of orbitals
void PipekMezey::localize(){
	size_t MAX_CYCLES = 100;

	double old_pm_val = functional();
	for (size_t i = 0; i <  MAX_CYCLES; i++){
		for (size_t j = 0; j < _norbitals; j++){
			for (size_t k = j + 1; k < _norbitals; k++){
				double angle = rot_ang(j, k);

				arma::vec j_col = _mo_coeff.col(j);
				arma::vec k_col = _mo_coeff.col(k);

				_mo_coeff.col(j) = j_col*cos(angle)+sin(angle)*k_col;
				_mo_coeff.col(k) = k_col*cos(angle)-sin(angle)*j_col;
			}
		}

		double pm_val = functional();
		double pm_diff = pm_val - old_pm_val;

		if (std::abs(pm_diff) < 1.0e-9){
			std::cout << "Pipek-Mezey Localization converged after " << i+1 << " cycles" << std::endl;
			return;
		}

		old_pm_val = pm_val;
	}

	throw std::runtime_error("Pipek-Mezey localization didn't converge after 100 cycles");
}

// calculate rotation angle
double PipekMezey::rot_ang(size_t a, size_t b){
	double Ast, Bst;

	Ast = 0;
	Bst = 0;

	for (size_t i = 0; i < _natoms; i++){
		auto rqast = qast(a, b, i);
		auto rqas  = qas(a, i);
		auto rqat  = qas(b, i);

		Ast += rqast*rqast - 0.25*(rqas-rqat)*(rqas-rqat);
		Bst += rqast*(rqas-rqat);
	}

	double rot_ang = 0.25*acos(-Ast/sqrt(Ast*Ast+Bst*Bst));
	if (Bst <= 0)
		rot_ang *= -1;

	return rot_ang;
}

double PipekMezey::qast(size_t a, size_t b, size_t atom){

	//get bf starting position
	size_t start = accumulate(_bf_per_atom.begin(), _bf_per_atom.begin() + atom, 0);
	size_t end = start + _bf_per_atom.at(atom) - 1;
	size_t bas = _mo_coeff.n_rows;

	double val = 0;
	double* over_ptr = _overlap.colptr(start);
	double* mo_a_ptr = _mo_coeff.colptr(a);
	double* mo_a_old = mo_a_ptr;
	double* mo_a_sub_ptr = mo_a_ptr + start;
	double* mo_b_ptr = _mo_coeff.colptr(b);
	double* mo_b_old = mo_b_ptr;
	double* mo_b_sub_ptr = mo_b_ptr + start;

	for (size_t i = start; i <= end; i++){
		double tmp_a = 0;
		double tmp_b = 0;
		// multiply column of mos by a column of the overlap matrix
		for (size_t j = 0; j < bas; j++){
			tmp_a += (*over_ptr) * (*mo_a_ptr++);
			tmp_b += (*over_ptr++) * (*mo_b_ptr++);
		}
		// multiply these results by a subset of the original mo column
		val += tmp_a * (*mo_b_sub_ptr++) + tmp_b * (*mo_a_sub_ptr++);
		mo_a_ptr = mo_a_old;
		mo_b_ptr = mo_b_old;
	}

	return val/2;
}

double PipekMezey::qas(size_t a, size_t atom){

	//get bf starting position
	size_t start = accumulate(_bf_per_atom.begin(), _bf_per_atom.begin() + atom, 0);
	size_t end = start + _bf_per_atom.at(atom) - 1;
	size_t bas = _mo_coeff.n_rows;

	double val = 0;
	double* over_ptr = _overlap.colptr(start);
	double* mo_ptr = _mo_coeff.colptr(a);
	double* mo_old = mo_ptr;
	double* mo_sub_ptr = mo_ptr + start;

	for (size_t i = start; i <= end; i++){
		double tmp = 0;
		for (size_t j = 0; j < bas; j++){
			// multiply column of mos by a column of the overlap matrix
			tmp += (*over_ptr++) * (*mo_ptr++);
		}
		// multiply these results by a subset of the original mo column
		val += tmp * (*mo_sub_ptr++);
		mo_ptr = mo_old;
	}

	// This way is more concise, using armadillo's libraries, but takes more than twice as long, due to temporary creation. I'll leave it here in case someone in the future wants to work with it. (The same thing can be used in qast)
	//    double val = dot((_mo_coeff.col(a)).subvec(start,end), (_overlap.cols(start,end)).t() * _mo_coeff.col(a));

	return val;
}

double PipekMezey::functional(){
	double pm_val = 0.0;
	for (size_t i = 0; i < _norbitals; i++){
		for (size_t j = 0; j < _natoms; j++){
			double val = qas(i, j);
			pm_val += val*val;
		}
	}

	return pm_val;
}
