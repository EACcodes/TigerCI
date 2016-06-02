/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
/*******************************************************************************
 *  Builds the local orthogonal virtual orbitals (LOVOs)
 * 
 *  Based on Fortran90 code by Jeremy Chwee, David Krisiloff,
 *          and Johannes M Dieterich
 * 
 *  Author: Francis Ricci
 ******************************************************************************/

#include <iostream>
#include <armadillo>
#include <stdexcept>
#include <sstream>

#include <tiger_cpp/local_ortho_virt.h>
#include <tiger_cpp/lovo_math.h>
#include <tiger_cpp/PipekMezey.h>
#include <tiger_cpp/OrbitalParsers.h>
#include <tiger_cpp/Timer.h>
#include <tiger_cpp/TigerStructs.h>

using namespace::std;
using namespace::arma;

/******************************************************************************/

// static function declarations
static mat min_bas_space(size_t min_bas, size_t num_bas, size_t num_internal,
        vector<size_t>& bf_per_atom, mat& overlap, mat& mos);

static vector<vector<size_t>> separate(BasisData *basis);

static vector<size_t> funcs(BasisData *basis);

static void localize_inactive_space(mat& orbitals, BasisData *basis,
        vector<size_t>& bf_per_atom, size_t num_inactive, size_t num_frozen);

static mat proto(size_t num_bas, size_t natom, mat& overlap,
        vector<vector<size_t>>& ang_bf, vector<vector<size_t>>& min_ang_bf);

static mat hard(size_t num_bas, size_t num_internal, size_t min_bas, 
        size_t natom, mat& mos, mat& overlap, mat& valence_virts, mat& protohard,
        vector<size_t>& bf_per_atom, vector<size_t>& min_bf_per_atom);

/******************************************************************************/

// Construct the local orthogonal virtual orbitals (LOVOs)
void ortho_virt(BasisData *basis){
    // THIS FUNCTION GENERATES ORTHONORMAL VIRTUAL ORBITALS
    // J. E. Subotnik, A. D. Dutoi, M. Head-Gordon, JCP, 123, 114108 (2005)
    //
    // Please note that if your AO space uses 5 d functions you need to make a
    // modification to the protohard virtuals
    
	if (basis->_globals->get_global_bool("skip_lovos")){
		cout << "WARNING: Skipping LOVO generation for a local calculation!" << endl;
		cout << "         If this is not desired, check your input file!" << endl;
		return;
	}

	Timer timer;

    // construct minimal basis
	BasisData *min_basis = new ErkaleBasis(*basis->_globals, tiger::minimal);

    cout << "Starting LOVO generation" << endl;
    
    // get some system information
    size_t natom = basis->get_nAtoms();
    size_t min_bas = min_basis->get_nBas();
    
    // get some global variables
    size_t num_internal = basis->_globals->get_global_int("num_internal");
    size_t num_bas = basis->_globals->get_global_int("number_bas");
    size_t num_inactive = basis->_globals->get_global_int("num_inactive");
    size_t num_frozen = basis->_globals->get_global_int("num_frozen");
    
    if (min_bas < num_internal){
        stringstream error;
        error << "LOVO ERROR" << '\n';
        error << "The size of the internal orbitals is greater than the";
        error << " size of the minimal basis space" << '\n';
        error << "Either decrease the number of references or increase";
        error << " the size of the minimal basis" << endl;
        throw runtime_error(error.str());
    }
    
    // preparation for minimal basis space calculation
    size_t tot_bas = min_bas + num_bas;
    mat overlap = zeros(tot_bas,tot_bas);
    overlap.submat(0,0,num_bas-1,num_bas-1) = basis->get_overlap();
    overlap.submat(num_bas,num_bas,tot_bas-1,tot_bas-1) = min_basis->get_overlap();
    overlap.submat(0,num_bas,num_bas-1,tot_bas-1) = basis->get_overlap(min_basis);
    overlap.submat(num_bas,0,tot_bas-1,num_bas-1) = basis->get_overlap(min_basis).t();
    
    mat mos = get_MOs();

    if (mos.n_cols != num_bas)
      throw runtime_error("Number of basis functions does not match MOs");
    
    vector<size_t> bf_per_atom = basis->funcs();
    vector<size_t> min_bf_per_atom = min_basis->funcs();
    vector<vector<size_t>> ang_bf = basis->separate();
    vector<vector<size_t>> min_ang_bf = min_basis->separate();
    
    mat valence_virts = min_bas_space(min_bas, num_bas, num_internal,
            bf_per_atom, overlap, mos);

    mat protohard = proto(num_bas, natom, overlap, ang_bf, min_ang_bf);
    
    mat hard_virts = hard(num_bas, num_internal, min_bas, natom, mos, overlap,
            valence_virts, protohard, bf_per_atom, min_bf_per_atom);
    
    if (basis->_globals->get_global_bool("localize_inactives"))
    	localize_inactive_space(mos, basis, bf_per_atom, num_inactive, num_frozen);
    else
    	cout << "WARNING: Inactive orbital localization has been turned off!" << endl;
    
    mos.cols(num_internal,min_bas-1) = valence_virts;
    mos.cols(min_bas,num_bas-1) = hard_virts;
    
    set_MOs(mos);

    delete min_basis;

    cout << endl;
    cout << "*******************************" << endl;
    cout << "* Orbital localization complete" << endl;
    cout << "*******************************" << endl;
    cout << endl;
    
    timer.print_elapsed_time("LOVO calculation    wall time:");
}

/******************************************************************************/

// Performs minimal basis space calculation
static mat min_bas_space(size_t min_bas, size_t num_bas, size_t num_internal,
        vector<size_t>& bf_per_atom, mat& overlap, mat& mos){
    
    size_t tot_bas = min_bas + num_bas;
    size_t red_bas = min_bas - num_internal;
    
    // Orthogonalize the current AO basis 
    mat ortho_vecs = eye(num_bas, num_bas);
    mat over_sub = overlap.submat(0,0,num_bas-1,num_bas-1);
    
    gram_schmidt(num_bas, ortho_vecs, over_sub);

    // Project minimal basis onto orthonormal AO basis, forming new overlap matrix
    over_sub = overlap.submat(0,num_bas,num_bas-1,tot_bas-1);
    mat proj_min_ao = project(min_bas, num_bas, over_sub, ortho_vecs);

    mat mos_sub = mos.cols(0,num_internal-1);
    over_sub = overlap.submat(0,0,num_bas-1,num_bas-1);
    mat proj_occ = proj_min_ao - (mos_sub * mos_sub.t() * over_sub * proj_min_ao);

    mat proj_over = proj_occ.t() * over_sub * proj_occ;
    proj_over = proj_over.submat(0,0,min_bas-1,min_bas-1);
    
    // Diagonalize overlap matrix, convert into valence virtual orbitals
    mat eigenvecs(min_bas,min_bas);
    vec eigenvals(min_bas);
    
    if (!eig_sym(eigenvals, eigenvecs, proj_over))
        throw runtime_error("Armadillo diagonalization failed.");
    
    mat virts = zeros(num_bas, red_bas);
    for (size_t i = 0; i < red_bas; i++)
        virts.col(i) = proj_occ*eigenvecs.col(i+num_internal);
    
    // normalization of valence virtuals
    normalize(virts, over_sub);
    for (size_t i = 0; i < red_bas; i++)
        virts.col(i) /= sqrt(dot((over_sub * virts.col(i)), virts.col(i)));

    PipekMezey pipek(bf_per_atom, virts, over_sub, red_bas, 0);
    // localize the valence virtuals, runtime_error if no convergence
    try {
      pipek.localize();
    }
    catch (runtime_error& err){
        cout << err.what() << endl;
        throw runtime_error("Problem localizing in min_bas_space!");
    }
    
    return virts;
}

/******************************************************************************/

static void localize_inactive_space(mat& orbitals, BasisData *basis,
    vector<size_t>& bf_per_atom, size_t num_inactive, size_t num_frozen){
    mat overlap = basis->get_overlap();
    mat orb_sub = orbitals.cols(num_frozen, num_inactive-1);
    
    PipekMezey pipek(bf_per_atom, orb_sub, overlap, num_inactive, num_frozen);
    try {
      pipek.localize();
    }
    catch (runtime_error& err){
        cout << err.what() << endl;
        throw runtime_error("Problem localizing in min_bas_space!");
    }
    
    orbitals.cols(num_frozen, num_inactive-1) = orb_sub;
}

/******************************************************************************/

// generate protohard virtual orbitals
static mat proto(size_t num_bas, size_t natom, mat& overlap,
        vector<vector<size_t>>& ang_bf, vector<vector<size_t>>& min_ang_bf){
    
    mat protohard = zeros(num_bas, num_bas);
    size_t iBas  = 0;
    size_t iMin  = 0;
    size_t iProt = 0;
      
    // iterate over all nuclei
    size_t inuc = 0;
    for (vector<size_t> nuc_bf : ang_bf){
        // looping over sets of shells having the same angular momentum component
        // different from the Fortran looping (ie. looping over px,py,pz separately)
        vector<size_t> min_nuc_bf = min_ang_bf.at(inuc);
        size_t tot_am = 0;
        for (size_t shell_funcs : nuc_bf){
            // orthonormalize AOs
            size_t end_bas = iBas+shell_funcs;
            mat ortho_vecs = eye(shell_funcs, shell_funcs);
            mat over_sub = overlap.submat(iBas,iBas,end_bas-1,end_bas-1);
            gram_schmidt(shell_funcs, ortho_vecs, over_sub);
 
            // project out minimal basis functions
            if (min_nuc_bf.size() > tot_am){
                size_t min_funcs = min_nuc_bf.at(tot_am);
                size_t end_min = iMin+num_bas+min_funcs;
                size_t red_funcs = shell_funcs - min_funcs;
                
                // project minimal AOs onto current AOs
                over_sub = overlap.submat(iBas,iMin+num_bas,end_bas-1,end_min-1);
                mat proj = project(min_funcs, shell_funcs, over_sub, ortho_vecs);
    
                // orthonormalize these AOs
                over_sub = overlap.submat(iBas, iBas, end_bas-1, end_bas-1);
                gram_schmidt(min_funcs, proj, over_sub);
            
                // project out the minimal basis components from current AO space
                mat density = proj * proj.t();
                proj = project_out(shell_funcs, shell_funcs, over_sub, density);

                // form overlap matrix
                mat proj_over = proj.t() * over_sub * proj;
                                
                // diagonalize the matrix
                mat eigenvecs(shell_funcs, shell_funcs);
                vec eigenvals(shell_funcs);    
                if (!eig_sym(eigenvals, eigenvecs, proj_over))
                    throw runtime_error("Armadillo diagonalization failed.");
                
                // form the new vectors
                mat vecs = zeros(shell_funcs, red_funcs);
                for (size_t i = 0; i < red_funcs; i++)
		  vecs.col(i) = proj * eigenvecs.col(min_funcs+i);
                
                // normalize the protohards
                normalize(vecs, over_sub);
                
                protohard.submat(iBas,iProt,end_bas-1,iProt+red_funcs-1) = vecs;
                
                iBas  += shell_funcs;
                iMin  += min_funcs;
                iProt += red_funcs;
            }            
            // no minimal basis functions to project out here
            else {
                size_t end_prot = iProt+shell_funcs;
		protohard.submat(iBas,iProt,end_bas-1,end_prot-1) = ortho_vecs;

                iBas  = end_bas;
                iProt = end_prot;
            }
            tot_am++;
        }
        inuc++;
    }
    
    return protohard;
}

/******************************************************************************/

// generate hard virtual orbitals
static mat hard(size_t num_bas, size_t num_internal, size_t min_bas, 
        size_t natom, mat& mos, mat& overlap, mat& valence_virts, mat& protohard,
        vector<size_t>& bf_per_atom, vector<size_t>& min_bf_per_atom){
    
    // The hard virtuals are constructed for each atom individually
    size_t iBas = 0;
    size_t hv_count = 0;
    mat hard_virts = zeros(num_bas, num_bas);
        
    /* Project all atomic orbitals for each atom onto the
       orthogonal complement of the minimal basis space. */
        
    // Make the projection operator
    mat proj_tmp = zeros(num_bas, min_bas);
    proj_tmp.cols(0,num_internal-1) = mos.cols(0,num_internal-1);
    proj_tmp.cols(num_internal,min_bas-1) = valence_virts.cols(0,min_bas-num_internal-1);
        
    mat proj_op = proj_tmp * proj_tmp.t();

    // iterate over all nuclei
    for (size_t inuc = 0; inuc < natom; inuc++){
        size_t nuc_funcs = bf_per_atom.at(inuc);
        size_t min_funcs = min_bf_per_atom.at(inuc);
        size_t hv_funcs  = nuc_funcs-min_funcs;

        // project out the minimal basis space from the current AO space
        mat over_sub = overlap.submat(0,iBas,num_bas-1,iBas+nuc_funcs-1);
	mat proj = -1 * proj_op * over_sub;

	for (size_t i = iBas; i < iBas+nuc_funcs; i++)
	  proj(i,i-iBas) = proj(i,i-iBas) + 1;

        over_sub = overlap.submat(0,0,num_bas-1,num_bas-1);
        mat proj_over = proj.t() * over_sub * proj;

        // diagonalize the matrix
        mat eigenvecs(nuc_funcs, nuc_funcs);
        vec eigenvals(nuc_funcs);    
        if (!eig_sym(eigenvals, eigenvecs, proj_over))
            throw runtime_error("Armadillo diagonalization failed.");
        
        // form and normalize initial hard virtuals
        for (size_t i = 0; i < hv_funcs; i++){
	  for (size_t j = 0; j < nuc_funcs; j++)
	    hard_virts.col(hv_count+i) += eigenvecs(j,i+min_funcs) * proj.col(j);
	  hard_virts.col(hv_count+i) /= sqrt(eigenvals(i+min_funcs));
        }
        
        iBas += nuc_funcs;
        hv_count += hv_funcs;
    }

    /* Now we want to form the transformation matrix, Z, which will orthogonalize
       and rotate the initial hard virtuals into maximum overlap with the
       protohard virtuals.
       Z = (R^-1 T)(T^T R^-1 T)^(-1/2)  */
    hv_count = 0;
    
    // Form the overlap matrix of the initial hard virtuals, R
    // Form the overlap matrix between the protohard and initial hard virtuals, T
    for (size_t inuc = 0; inuc < natom; inuc++){
        size_t nuc_funcs = bf_per_atom.at(inuc);
        size_t min_funcs = min_bf_per_atom.at(inuc);
        size_t hard_funcs = nuc_funcs - min_funcs;
        
        mat over_sub = overlap.submat(0,0,num_bas-1,num_bas-1);
        mat hard_sub = hard_virts.cols(hv_count,hv_count+hard_funcs-1);
        mat proto_sub = protohard.cols(hv_count,hv_count+hard_funcs-1);
        mat over_hard = over_sub * hard_sub;

        mat T = (over_hard).t() * proto_sub;
        mat R = (over_hard).t() * hard_sub;
        
        // find inverse of R
        mat R_inv = zeros(hard_funcs,hard_funcs);
        try{
            R_inv = inv(R, "std");
        }
        catch (runtime_error& err){
            cout << err.what() << endl;
            throw runtime_error("Armadillo failed to invert matrix R!");
        }
        
        // form T^T R^-1 T
        mat RT = R_inv * T;
        mat TRT = T.t() * RT;

        // take inverse sqrt
        vec eigenvals = zeros(hard_funcs);
        mat eigenvecs = zeros(hard_funcs,hard_funcs);
        diag_inv_sqrt(hard_funcs, eigenvals, eigenvecs, TRT);
        
        // form Z
        mat Z = RT * eigenvecs;
        
        // get the hard virtuals
        mat temp_hard = zeros(num_bas,hard_funcs);
        for (size_t i = 0; i < hard_funcs; i++){
            for (size_t j = 0; j < hard_funcs; j++)
                temp_hard.col(i) += Z(i,j) * hard_virts.col(hv_count+j);
        }
        
        hard_virts.cols(hv_count,hv_count+hard_funcs-1) = temp_hard.cols(0,hard_funcs-1);
	hv_count += hard_funcs;
    }        

    // normalize the new hard virtuals
    mat over_sub = overlap.submat(0,0,num_bas-1,num_bas-1);
    hard_virts.resize(num_bas,hv_count);
    normalize(hard_virts, over_sub);
        
    // form the overlap matrix
    mat over_hard = hard_virts.t() * over_sub * hard_virts;
        
    vec eigenvals = zeros(hv_count);
    mat eigenvecs = zeros(hv_count, hv_count);
    diag_inv_sqrt(hv_count, eigenvals, eigenvecs, over_hard);

    // Form new orthogonalized hard virtuals
    hard_virts *= eigenvecs;
    
    return hard_virts;
}
