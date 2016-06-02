/* 
 * File:   ErkaleData.h
 * Author: Francis Ricci
 *
 * Created on December 9, 2014, 2:11 PM
 */

#ifndef BASIS_DATA_H
#define	BASIS_DATA_H

#include <vector>
#include <armadillo>
#include <memory>

#include <erkale/basis.h>
#include <erkale/eriworker.h>

#include <tiger_cpp/TigerStructs.h>

class ErkaleData{
public:
	ErkaleData(std::shared_ptr<BasisSet> bas);

	double get_nuclear_repulsion();
	std::vector<tiger::atom_t> get_atoms();
	std::vector<tiger::coords_t> get_coordinates();
	std::vector<int> get_nFuncInShell();
	int get_nBas();
	int get_nShell();
	int get_nAtoms();
	std::vector<int> get_basis_map();
	std::vector<std::vector<size_t>> separate();
	std::vector<size_t> funcs();
    arma::vec eval_func(double x, double y, double z);
	void print_basis();

	// Shell information
	void print_shell(size_t ind);
	coords_t get_center(size_t ind);
	std::vector<contr_t> get_contr_normalized(size_t ind);
	std::vector<contr_t> get_contr(size_t ind);
	bool lm_in_use(size_t ind);
	arma::mat get_trans(size_t ind);
	int get_am(size_t ind);
	std::vector<shellf_t> get_cart(size_t ind);
	size_t get_Nbf(size_t ind);
	size_t get_first_ind(size_t ind);

private:
	std::shared_ptr<BasisSet> _basis;
};

#endif	/* BASIS_DATA_H */
