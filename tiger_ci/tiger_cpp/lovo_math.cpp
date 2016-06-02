/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
/* 
 * File:   lovo_math.cpp
 * Author: Francis Ricci
 *
 * Created on August 4, 2014, 11:52 AM
 */

#include <math.h>
#include <numeric>
#include <stdexcept>

#include "lovo_math.h"

using namespace::arma;
using namespace::std;

// Performs modified Gram Schmidt orthogonalization on vecs.
// \param num_vec Number of vectors
// \param vecs    Vectors to be orthogonalized
// \param overlap Vector basis overlap
void gram_schmidt(size_t num_vec, mat& vecs, mat& overlap){
    
    for (size_t i = 0; i < num_vec; i++){

        // Normalize vector i, in orthogonalized basis. (Requires temporary mem)
        vecs.col(i) /= sqrt(dot((overlap * vecs.col(i)), vecs.col(i)));
        vec prod = overlap * vecs.col(i);
    
        // Remove projection of this vector from all remaining vectors
        for (size_t j = i + 1; j < num_vec; j++)
            vecs.col(j) -= dot(vecs.col(j), prod) * vecs.col(i);
    }
}

// Projection of matrix a onto matrix b
// \param len_a Size of matrix a basis
// \param len_b Size of matrix b basis
// \param a     Matrix a
// \param b     Matrix b
mat project(size_t len_a, size_t len_b, mat& a, mat& b){
    mat project = zeros(len_b, len_a);
    
    for (size_t i = 0; i < len_a; i++){
        for (size_t j = 0; j < len_b; j++){
            project.col(i) += dot(b.col(j), a.col(i)) * b.col(j);
        }
    }
    
    return project;
}

// Project out matrix b, with overlap given by matrix a
mat project_out(size_t len_a, size_t len_b, mat& a, mat& b){
    mat project = zeros(len_b, len_a);
    
    for (size_t i = 0; i < len_a; i++){
        for (size_t j = 0; j < len_b; j++){
            if (i == j)
                project(j,i) = 1 - dot(b.row(j), a.col(i));
            else
                project(j,i) = -1 * dot(b.row(j), a.col(i));
        }
    }
    
    return project;
}

// Normalize input matrix, with given overlap
void normalize(mat& matrix, mat& overlap){

  size_t ncols = matrix.n_cols;
    for (size_t i = 0; i < ncols; i++)
        matrix.col(i) /= sqrt(dot(matrix.col(i), overlap * matrix.col(i)));
}

// Diagonalize matrix by solving for the eigenvalues, forming inverse of matrix
void diag_inv_sqrt(size_t len, vec& eigenvals, mat& eigenvecs, mat& matrix){
    // Find inverse of input matrix
    mat mat_inv = zeros(len, len);
    try{
        mat_inv = inv(matrix, "std");
    }
    catch (runtime_error& err){
        cout << err.what() << endl;
        throw runtime_error("Armadillo failed to invert matrix!");
    }
    
    // diagonalize inverted matrix
    if (!eig_sym(eigenvals, eigenvecs, mat_inv))
        throw runtime_error("Armadillo diagonalization failed.");
    
    for (size_t i = 0; i < len; i++)
        eigenvals(i) = sqrt(eigenvals(i));
    
    // inverse square root matrix
    eigenvecs = eigenvecs * diagmat(eigenvals) * eigenvecs.t();
}
