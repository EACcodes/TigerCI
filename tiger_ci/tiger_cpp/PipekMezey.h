/* 
 * File:   pipek_mezey.h
 * Author: fjricci
 *
 * Created on August 25, 2014, 11:41 AM
 */

#ifndef PIPEK_MEZEY_H
#define	PIPEK_MEZEY_H

#include <armadillo>
#include <vector>

class PipekMezey {
public:
	PipekMezey(std::vector<size_t>& bf_per_atom, arma::mat& mo_coeff, arma::mat& overlap,
			size_t norbitals, size_t nfrozen);

	// Performs Pipek_Mezey localization of orbitals
	void localize();

private:
	std::vector<size_t>& _bf_per_atom;
	arma::mat& _mo_coeff;
	arma::mat& _overlap;
	size_t _norbitals;
	size_t _natoms;

	//evaluate the current degree of localization
	double functional();

	//calculate orbital rotation angle for this iteration
	double rot_ang(size_t a, size_t b);

	//calculate mulliken charge across 2 orbitals
	double qast(size_t a, size_t b, size_t atom);

	//calculate gross mulliken charge on this orbital and atom
	double qas(size_t a, size_t atom);
};

#endif	/* PIPEK_MEZEY_H */

