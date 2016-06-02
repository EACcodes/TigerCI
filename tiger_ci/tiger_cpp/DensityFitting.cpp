/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <assert.h>
#include <inttypes.h>
#ifdef TIGER_USE_OMP
#include <omp.h>
#endif

#include <tiger_cpp/Timer.h>
#include <tiger_cpp/DensityFitting.h>

using namespace std;

extern "C" void dposv_(char *uplo, int64_t *n, int64_t *nrhs, double* a, int64_t *lda, double* b, int64_t *ldb, int64_t *info);

DensityFitting::DensityFitting(BasisData *bas, const arma::mat& orbitals):
        basis(bas),
        nbf(basis->get_nBas()),
        naux(basis->get_nAux()),
        norb(orbitals.n_cols),
        MOs(orbitals){}

void DensityFitting::compute_fitting(){    
    Timer df_time; 
    cout << endl;
    cout << "***********************" << endl;
    cout << "*                      " << endl;
    cout << "*   Density Fitting    " << endl;
    cout << "*                      " << endl;
    cout << "***********************" << endl;
    cout << endl;
    cout << "  Number of auxiliary basis functions = " << naux << endl;

    // Calculate (A|B) and (A|B)^(1/2)
    arma::mat ab;
    basis->compute_two_center(ab);
    arma::mat ab_sqrt;
    compute_two_center_chol(ab, ab_sqrt);
    
    // Calculate the expansion coefficients
    calculate_aux_expansion_coeff(ab);
    
    // Calculate the contraction of (A|B)^(1/2) and the coefficients
    compute_symmetric_contraction(ab_sqrt);
    df_time.print_elapsed_time("  Density Fitting");
    
    // Transform to the MO basis 
    transform_to_MO_basis();
}

void DensityFitting::compute_two_center_chol(const arma::mat& ab, arma::mat& ab_sqrt) const{
    Timer t;
    ab_sqrt = arma::chol(ab);
    t.print_elapsed_time("  (A|B) cholesky time");
}



void DensityFitting::transform_to_MO_basis(){
    /*
     * Transform the auxiliary basis expansion coefficients to the MO basis
     * Using two steps, similar to how the cholesky vectors are transformed
     */
    cout << endl;
    cout << "****************************************" << endl;
    cout << "*                                       " << endl;
    cout << "*   AO->MO Density Fitting Transform    " << endl;
    cout << "*                                       " << endl;
    cout << "****************************************" << endl;
    cout << endl;
    
    Timer df_transform_time;
    size_t norb_pairs = norb * (norb + 1)/2;
    PairIndex map_munu(nbf, nbf, true);
    PairIndex map_ij(norb, norb, true);
   
    // Equation (1)
    // X^A_mui = \sum_nu c_inu * C^A_munu
    arma::Cube<double> tmp(naux,nbf, norb);
    
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tmp, map_munu)
#endif
    for (size_t i=0; i<norb; i++){
        for (size_t mu=0; mu<nbf; mu++){
            tmp.slice(i).col(mu).zeros();
            for (size_t nu=0; nu<nbf; nu++){
                size_t addr = map_munu.get_index(nu,mu);
                tmp.slice(i).col(mu) += MOs(nu,i) * coefficients.col(addr);
            }
        }
    }
   
    // Equation (2) 
    // Y^A_ji = \sum_mu c_jmu * X^A_mui
    coefficients.set_size(naux, norb_pairs);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tmp, map_ij)
#endif
    for (size_t i=0; i<norb; i++){
        for (size_t j=0; j<=i; j++){
            size_t mo_index = map_ij.get_index(i,j);
            coefficients.col(mo_index).zeros();
            for (size_t mu=0; mu<nbf; mu++){
                coefficients.col(mo_index) += MOs(mu,j) * tmp.slice(i).col(mu); 
            }
        }
    }
 
    df_transform_time.print_elapsed_time("  DF AO->MO transform");
    cout << endl;
}

void DensityFitting::compute_symmetric_contraction(const arma::mat& ab_sqrt){
    Timer t;   
    PairIndex map(nbf, nbf, true);

    arma::mat trans_ab = ab_sqrt.t();

   // coeff^A_munu = \sum_B tmp^B_munu * (A|B)^(1/2)
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(ab_sqrt,map)
#endif
    for (size_t mu=0; mu<nbf; ++mu){
        for (size_t nu=0; nu<=mu; ++nu){
	    size_t index = map.get_index(mu,nu);
	    coefficients.col(index) = ab_sqrt * coefficients.col(index);
        }
    }
    t.print_elapsed_time("  symmetric contraction time");
}


void DensityFitting::calculate_aux_expansion_coeff(arma::mat& ab){
    PairIndex addr(nbf, nbf, true);
    
    // compute the integrals (A|\mu \nu) and store them in coefficients (to save memory)
    Timer t1;
    basis->compute_three_center(coefficients, addr);
    t1.print_elapsed_time("  (A|mu nu) time");
    
    // We solve for the coefficients via the equation
    // \f$\sum_A (A|B) C^A_\mu\nu = (A|\mu \nu)\f$
    // This is just a linear equation of the form 
    // A * X = B
    // 
    // To solve this reliably we are going to cholesky
    // decompose A = L * L^T  (effectively an LU decomp for 
    // a symmetric matrix)
    // (L * L^T) * X = (L * U) * X = B
    //
    // Then we solve in two steps, let Y = U * X.
    // Solve for Y
    // L * Y = B
    // Then solve for X from Y
    // U * X = Y

    Timer t2;
    lapack_solve(ab);
    t2.print_elapsed_time("  expansion coeff time");
}

void DensityFitting::lapack_solve(arma::mat& ab){
  size_t nbfpairs = nbf * (nbf+1) / 2 ; 
  char uplo = 'L';
  int64_t naux_8 = (int64_t) naux;

  // If in parallel, split calculation up into blocks of equal size
#ifdef _OPENMP
  int blockID = 0;
  int numThreads = omp_get_max_threads();
  int colPerThread = nbfpairs / numThreads;

  vector<pair<int, int>> endpoints;

  for (int thread = 0; thread < numThreads; ++thread){
    int end;
    int start = blockID;
    blockID += colPerThread;
    if (thread == numThreads - 1)
      end = coefficients.n_cols - 1;
    else
      end = blockID - 1;
    pair<int,int> ends(start, end);
    endpoints.push_back(ends);    
  }

#pragma omp parallel default(none) shared(endpoints, nbfpairs, blockID, uplo, naux_8, ab)
  {
    pair<int,int> ends = endpoints.at(omp_get_thread_num());
    int start = ends.first;
    int end = ends.second;
    arma::mat my_ab(ab);
    int64_t nbfpairs_8 = (int64_t) end-start+1;
    int64_t info = 0;
    dposv_(&uplo, &naux_8, &nbfpairs_8, my_ab.memptr(), &naux_8, coefficients.colptr(start), &naux_8, &info);
    if (info < 0){
      ostringstream oss;
      oss << info;
      throw runtime_error("Argument " + oss.str() + " to dposv had an illegal value.");
    }
    if (info > 0){
      ostringstream oss;
      oss << info;
      throw runtime_error("Leading minor of order " + oss.str() + " of AB is not positive definite.");
    }
  }
#else
  int64_t nbfpairs_8 = (int64_t) nbfpairs;
  int64_t info = 0;
  dposv_(&uplo, &naux_8, &nbfpairs_8, ab.memptr(), &naux_8, coefficients.memptr(), &naux_8, &info);
  if (info < 0){
    ostringstream oss;
    oss << info;
    throw runtime_error("Argument " + oss.str() + " to dposv had an illegal value.");
  }
  if (info > 0){
    ostringstream oss;
    oss << info;
    throw runtime_error("Leading minor of order " + oss.str() + " of AB is not positive definite.");
  }
#endif

}


void DensityFitting::get_DF_block(size_t i, size_t j, std::vector<double>& values) const{   
    PairIndex map(norb, norb, true);
    size_t addr = map.get_index(i,j);
    assert(naux==values.size());
    for(size_t a=0; a<naux; a++){
        values[a] = coefficients(a, addr);
    }
}

size_t DensityFitting::get_naux() const{ 
    return naux;
}

size_t DensityFitting::get_norb() const{
    return norb;
}

void DensityFitting::debug_MO_3cent_contracted(arma::mat& C, PairIndex& map){   
    // Building the 2-e ints in the MO basis from the DF quantities
    // note that I'm not directly comparing to the analytic ints here
    // usually I print those out in an another job and compare by hand
    double y;
      
    for (size_t i =0; i <norb; i++){
        for (size_t j=0 ; j<=i ; j++){
            for (size_t k=0; k<norb; k++){
                for(size_t l=0; l<=k; l++){
                    
                    y = 0.0;
                    
                    size_t ind_ij = map.get_index(i,j);
                    size_t ind_kl = map.get_index(k,l);
                    for (size_t a=0; a<naux; a++){
                        y += C(a,ind_ij)*C(a,ind_kl);
                    }
                    
    cout << endl;
    cout << "( " << i << " " << j << " | ";
    cout << k << " " << l << " )" << " ";
    cout << " = " << y << endl;                           
                }
            }
        }
    }           
}



void DensityFitting::debug_AO_3cent_contracted(arma::mat& B, PairIndex& map){
   // Build the exact integrals and then check how well we approximate them

    double x,y;
	vector<int> nFunc = basis->get_nFuncInShell();
    for (size_t i =0; i <nbf; i++){
		int iShell = bas_to_shell(i, nFunc);
		int iIdx = bas_in_shell(i, nFunc);
        for (size_t j=0 ; j<=i ; j++){
			int jShell = bas_to_shell(j, nFunc);
			int jIdx = bas_in_shell(j, nFunc);
            for (size_t k=0; k<nbf; k++){
				int kShell = bas_to_shell(k, nFunc);
				int kIdx = bas_in_shell(k, nFunc);
                for(size_t l=0; l<=k; l++){
					int lShell = bas_to_shell(l, nFunc);
					int lIdx = bas_in_shell(l, nFunc);
                    
					arma::mat eris = basis->eval_ijkl(iShell, jShell, kShell, lShell);
					x = eris(iIdx * nFunc[iShell] + jIdx, kIdx * nFunc[kShell] + lIdx);
                    y = 0.0;
                    size_t ind_ij = map.get_index(i,j);
                    size_t ind_kl = map.get_index(k,l);
                    for (size_t a=0; a<naux; a++){
                        y += B(a,ind_ij)*B(a,ind_kl);
                    }
                    // output here
                    cout << endl;
                    cout << "( " << i << " " << j << " | ";
                    cout << k << " " << l << " )" << " ";
                    cout << "exact = " << x << " df = " << y << endl;
                }
            }
        }
    }
}

void DensityFitting::debug_AO_3cent(arma::mat& C, arma::mat& ab, PairIndex& map){
   // Build the exact integrals and then check how well we approximate them

    double x,y;
	vector<int> nFunc = basis->get_nFuncInShell();
    for (size_t i =0; i <nbf; i++){
		int iShell = bas_to_shell(i, nFunc);
		int iIdx = bas_in_shell(i, nFunc);
        for (size_t j=0 ; j<=i ; j++){
			int jShell = bas_to_shell(j, nFunc);
			int jIdx = bas_in_shell(j, nFunc);
            for (size_t k=0; k<nbf; k++){
				int kShell = bas_to_shell(k, nFunc);
				int kIdx = bas_in_shell(k, nFunc);
                for(size_t l=0; l<=k; l++){
					int lShell = bas_to_shell(l, nFunc);
					int lIdx = bas_in_shell(l, nFunc);
                    
					arma::mat eris = basis->eval_ijkl(iShell, jShell, kShell, lShell);
					x = eris(iIdx * nFunc[iShell] + jIdx, kIdx * nFunc[kShell] + lIdx);
					y = 0.0;
                    
                    size_t ind_ij = map.get_index(i,j);
                    size_t ind_kl = map.get_index(k,l);
                    for (size_t a=0; a<naux; a++){
                        double c_aij = C(a, ind_ij);
                        for (size_t b=0; b<naux; b++){
                            double c_bkl = C(b, ind_kl);
                            y += c_aij * ab(a,b) * c_bkl;
                        }
                    }                   
                    cout << endl;
                    cout << "( " << i << " " << j << " | ";
                    cout << k << " " << l << " )" << " ";
                    cout << "exact = " << x << " df = " << y << endl;                           
                }
            }
        }
    }           
}
