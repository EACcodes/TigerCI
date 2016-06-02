/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/*
 * ErkaleDF.cpp
 *
 *  Created on: Jan 24, 2015
 *      Author: fjricci
 */

#include <erkale/density_fitting.h>
#include <erkale/eriworker.h>

#include <tiger_cpp/ErkaleDF.h>
#include <tiger_cpp/Timer.h>

ErkaleDF::ErkaleDF(std::shared_ptr<BasisSet> bas, std::shared_ptr<BasisSet> auxbas):
	_basis(bas),
	_auxbasis(auxbas){}

void ErkaleDF::compute_three_center_integrals(arma::mat& integrals, PairIndex& addr){
    /* This code is almost entirely a copy/paste from Erkale
     * Normally I wouldn't do this, but Erkale doesn't provide
     * a good API for accessing this info, (A|\mu\nu), directly
     */
	int nbf = _basis->get_Nbf();
	int naux = _auxbasis->get_Nbf();

    size_t nbf_pairs = nbf * (nbf+1) / 2;
    integrals.set_size(naux, nbf_pairs);

    size_t maxam = std::max(_basis->get_max_am(), _auxbasis->get_max_am());
    size_t maxcontr = std::max(_basis->get_max_Ncontr(), _auxbasis->get_max_Ncontr());
    const std::vector<double>*erip;

    auto orbshells=_basis->get_shells();
    auto auxshells=_auxbasis->get_shells();
    auto dummy=dummyshell();
    auto orbpairs=_basis->get_unique_shellpairs();

#ifdef _OPENMP
#pragma omp parallel default(none) shared(addr, integrals, maxam, maxcontr, orbpairs, dummy, orbshells, auxshells) private(erip)
    {
    ERIWorker eri(maxam,maxcontr);
#pragma omp for
#else
    ERIWorker eri(maxam,maxcontr);
#endif
    for (size_t ia = 0; ia < auxshells.size(); ia++){
        for (size_t ip = 0; ip < orbpairs.size(); ip++) {
            // Shells in question are
            size_t imu = orbpairs[ip].is;
            size_t inu = orbpairs[ip].js;

            // Amount of functions

            size_t Na = auxshells[ia].get_Nbf();
            size_t Nmu = orbshells[imu].get_Nbf();
            size_t Nnu = orbshells[inu].get_Nbf();

            // Compute (a|mn)
            eri.compute(&auxshells[ia], &dummy, &orbshells[imu], &orbshells[inu]);
            erip = eri.getp();

            // Store integrals
            for (size_t af = 0; af < Na; af++) {
                size_t inda = auxshells[ia].get_first_ind() + af;

                for (size_t muf = 0; muf < Nmu; muf++) {
                    size_t indmu = orbshells[imu].get_first_ind() + muf;

                    for (size_t nuf = 0; nuf < Nnu; nuf++) {
                        size_t indnu = orbshells[inu].get_first_ind() + nuf;

                        integrals(inda, addr.get_index(indmu, indnu)) = (*erip)[(af * Nmu + muf) * Nnu + nuf];
                    }
                }
            }
        }
    }
#ifdef _OPENMP
    }
#endif
}

void ErkaleDF::compute_two_center_integrals(arma::mat& ab){
    // This code is more or less directly from Erkale
    Timer t;
    size_t maxauxam=_auxbasis->get_max_am();
    size_t maxauxcontr=_auxbasis->get_max_Ncontr();
    auto aux_pairs = _auxbasis->get_unique_shellpairs();
    auto aux_shells = _auxbasis->get_shells();
    auto dummy = dummyshell();

	int naux = _auxbasis->get_Nbf();
    ab.set_size(naux, naux);

    ERIWorker eri(maxauxam,maxauxcontr);
    for (auto pair: aux_pairs){
        size_t A = pair.is;
        size_t B = pair.js;

        // compute for shells A,B (A|B)
        eri.compute(&aux_shells[A],&dummy,&aux_shells[B],&dummy);
        auto erip = eri.getp();

        // store the integrals in ab
        size_t Ni = aux_shells[A].get_Nbf();
        size_t Nj = aux_shells[B].get_Nbf();
        for (size_t ii=0; ii<Ni; ii++){
           size_t a=aux_shells[A].get_first_ind()+ii;
           for(size_t jj=0;jj<Nj;jj++) {
               size_t b=aux_shells[B].get_first_ind()+jj;
               ab(a,b)=(*erip)[ii*Nj+jj];
	       ab(b,a)=(*erip)[ii*Nj+jj];
	}
      }
    }
    t.print_elapsed_time("  (A|B) time");
}

int ErkaleDF::get_nAux(){
	return _auxbasis->get_Nbf();
}
