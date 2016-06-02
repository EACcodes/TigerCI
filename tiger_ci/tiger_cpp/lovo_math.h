/* 
 * File:   lovo_math.h
 * Author: Francis Ricci
 *
 * Created on August 4, 2014, 11:52 AM
 */

#ifndef LOVO_MATH_H
#define	LOVO_MATH_H

#include <armadillo>
#include <vector>

using namespace::arma;
using namespace::std;

// Performs modified Gram Schmidt orthogonalization on vecs.
void gram_schmidt(size_t num_vec, mat& vecs, mat& overlap);

// Projection of matrix b, with overlap given by matrix a
mat project(size_t len_a, size_t len_b, mat& a, mat& b);

// Project out matrix b, with overlap given by matrix a
mat project_out(size_t len_a, size_t len_b, mat& a, mat& b);

// Diagonalize matrix by solving for the eigenvalues, forming inverse of matrix
void diag_inv_sqrt(size_t len, vec& eigenvals, mat& eigenvecs, mat& matrix);

// Normalize input matrix, with given overlap
void normalize(mat& matrix, mat& overlap);

#endif	/* LOVO_MATH_H */

