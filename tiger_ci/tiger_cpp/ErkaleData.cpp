/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * @file  BasisData.cpp
 * Author: Francis Ricci
 *
 * Created on December 9, 2014, 2:11 PM
 */

#include <stdexcept>

#include <tiger_cpp/ErkaleData.h>

using namespace std;

ErkaleData::ErkaleData(std::shared_ptr<BasisSet> bas){
    _basis = bas;
}

double ErkaleData::get_nuclear_repulsion(){
    return _basis->Enuc();
}

std::vector<tiger::atom_t> ErkaleData::get_atoms(){
	vector<tiger::coords_t> coords = get_coordinates();
	vector<tiger::atom_t> atoms;
	for (size_t i = 0; i < coords.size(); i++){
		tiger::atom_t atom;
		atom.coords.x = coords[i].x;
		atom.coords.y = coords[i].y;
		atom.coords.z = coords[i].z;
		atom.Z = _basis->get_Z(i);
		atoms.push_back(atom);
	}
	return atoms;
}

vector<tiger::coords_t> ErkaleData::get_coordinates(){
    arma::mat coord_mat = _basis->get_nuclear_coords();
    vector<tiger::coords_t> coords;
    for (size_t i = 0; i < coord_mat.n_rows; i++){
    	tiger::coords_t coord;
    	coord.x = coord_mat(i,0);
    	coord.y = coord_mat(i,1);
    	coord.z = coord_mat(i,2);
    	coords.push_back(coord);
    }
    return coords;
}

vector<int> ErkaleData::get_nFuncInShell(){
    int num_shells = get_nShell();
    vector<GaussianShell> shells = _basis->get_shells();
    
    vector<int> nFunc = vector<int>(num_shells);
    for (int i = 0; i < num_shells; i++){
        nFunc[i] = shells[i].get_Nbf();
    }
    return nFunc;
}

int ErkaleData::get_nBas(){
    return _basis->get_Nbf();
}

int ErkaleData::get_nShell(){
    return _basis->get_Nshells();
}

int ErkaleData::get_nAtoms(){
    return _basis->get_Nnuc();
}

vector<int> ErkaleData::get_basis_map(){
    int num_bas = get_nBas();
    vector<int> basis_map = vector<int>(num_bas);
    vector<GaussianShell> shells = _basis->get_shells();

    for (GaussianShell shell : shells){
        int first = shell.get_first_ind();
        int last = shell.get_last_ind();
        int center = shell.get_center_ind();
        
        for (int i = first; i <= last; i++){
            basis_map[i] = center + 1; //+1 because fortran array indices
        }
    }
    
    return basis_map;
}

// separate shells by angular momentum
vector<vector<size_t>> ErkaleData::separate(){

    vector<vector<size_t>> ang_bf;
    vector<size_t> temp_shells;

    size_t tot_nuc = get_nAtoms();

    for (size_t inuc = 0; inuc < tot_nuc; inuc++){
        vector<GaussianShell> shells  = _basis->get_funcs(inuc);
        size_t num_funcs = 0;
        int am = 0;
        for (GaussianShell shell : shells){
            if (shell.get_am() == am)
                num_funcs += shell.get_Nbf();
            else {
                am = shell.get_am();
                temp_shells.push_back(num_funcs);
                num_funcs = shell.get_Nbf();
            }
        }
        if (num_funcs != 0)
            temp_shells.push_back(num_funcs);
        ang_bf.push_back(temp_shells);
        temp_shells.clear();
    }

    return ang_bf;
}

//get number of basis functions per atom
vector<size_t> ErkaleData::funcs(){
    size_t natom = get_nAtoms();
    vector<size_t> bf_per_atom;
    for (size_t inuc = 0; inuc < natom; inuc++){
        vector<GaussianShell> shells = _basis->get_funcs(inuc);
        size_t num_funcs = 0;
        for (GaussianShell shell : shells)
            num_funcs += shell.get_Nbf();
        bf_per_atom.push_back(num_funcs);
    }
    return bf_per_atom;
}

//evaluate basis function at x,y,z
arma::vec ErkaleData::eval_func(double x, double y, double z){
    return _basis->eval_func(x, y, z);
}


// print basis set information
void ErkaleData::print_basis(){
	_basis->print();
}

// print shell with index ind
void ErkaleData::print_shell(size_t ind){
	GaussianShell shell = _basis->get_shell(ind);
	shell.print();
}

coords_t ErkaleData::get_center(size_t ind){
	return _basis->get_shell(ind).get_center();
}

std::vector<contr_t> ErkaleData::get_contr_normalized(size_t ind){
	return _basis->get_shell(ind).get_contr_normalized();
}

std::vector<contr_t> ErkaleData::get_contr(size_t ind){
	return _basis->get_shell(ind).get_contr();
}

bool ErkaleData::lm_in_use(size_t ind){
	return _basis->get_shell(ind).lm_in_use();
}

arma::mat ErkaleData::get_trans(size_t ind){
	return _basis->get_shell(ind).get_trans();
}

int ErkaleData::get_am(size_t ind){
	return _basis->get_shell(ind).get_am();
}

std::vector<shellf_t> ErkaleData::get_cart(size_t ind){
	return _basis->get_shell(ind).get_cart();
}

size_t ErkaleData::get_Nbf(size_t ind){
	return _basis->get_shell(ind).get_Nbf();
}

size_t ErkaleData::get_first_ind(size_t ind){
	return _basis->get_shell(ind).get_first_ind();
}


