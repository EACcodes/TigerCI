/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * @file   IntegralScreen.cpp
 * @author Francis Ricci
 *
 * Created on November 24, 2014, 11:29 AM
 */

#include <iostream>

#include "IntegralScreen.h"
#include "FunctionMap.h"

#include <math.h>

using namespace std;

IntegralScreen::IntegralScreen(double ao_thresh, double cd_thresh, BasisData *data,
                               RestrictedIndex *pair_map){
    _cd_thresh = cd_thresh;
    _ao_thresh = ao_thresh;
    _nScreen = 0;
    _basis_data = data;
    _pair_map = pair_map;
    _nFunc = _basis_data->get_nFuncInShell();
    _idx = new PairIndex(data->get_nBas(), data->get_nBas(), true);
    calc_diag_ints();
    setup_screening();
}

IntegralScreen::~IntegralScreen(){
  delete _idx;
}

/* @brief This subroutine sets up the data structure used for
 *        Cauchy-Schwarz (CS) screening
 * @param diag       The diagonal two electron integrals (ij|ij)
 * @param integrals  An object for getting data bout integrals
 */
void IntegralScreen::setup_screening(){
    //allocate and zero screening arrays
    int nBas = _basis_data->get_nBas();
    int nShell = _basis_data->get_nShell();
    _Q_int.zeros(nBas, nBas);
    _Q_shell.zeros(nShell, nShell);
    
    // First step, building the Q screening quantities at the individual
    // basis set level
    for (int iBas = 0; iBas < nBas; iBas++){
        for (int jBas = 0; jBas <= iBas; jBas++){
            int iDiag = diag_address(iBas, jBas);
            double tmp = sqrt(_diag[iDiag]);
            
            _Q_int(iBas, jBas) = tmp;
            _Q_int(jBas, iBas) = tmp;
        }
    }
    
    vector<int> shell_map(nBas);
    vector<int> nFunc = _basis_data->get_nFuncInShell();
    int base = 0;
    for (int i = 0; i < nShell; i++){
        for (int j = base; j < base + nFunc[i]; j++)
            shell_map[j] = i;
	base += nFunc[i];
    }

    // Now we build the Q screening quantities at the shell level
    for (int iBas = 0; iBas < nBas; iBas++){
        int iShell = shell_map[iBas];
        for (int jBas = 0; jBas < nBas; jBas++){
            int jShell = shell_map[jBas];

            double screen = _Q_int(iBas, jBas);
            if (screen > _Q_shell(iShell, jShell)){
                _Q_shell(iShell, jShell) = screen;
                _Q_shell(jShell, iShell) = screen;
            }
        }
    }
}

/* @brief Answers the question is this integral shell block screened
 * @param i      The I in the integral block (IJ|KL) to be screened
 * @param j      The J in the integral block (IJ|KL) to be screened
 * @param k      The K in the integral block (IJ|KL) to be screened
 * @param l      The L in the integral block (IJ|KL) to be screened
 * @param nFunc  The number of AO functions in each shell
 * @return bool   Is integral shell screened?
 */
bool IntegralScreen::is_screened(int i, int j, int k, int l){
    double bound = _Q_shell(i,j) * _Q_shell(k,l);
    if (bound > _ao_thresh) return false;
    //    _nScreen += _nFunc.at(i) * _nFunc.at(j) * _nFunc.at(k) * _nFunc.at(l);
    return true;
}

/* @brief Has this pair been prescreened out
 * @param i   First member of pair
 * @param j   Second member of pair
 * @return bool   Is this pair screened?
 */
bool IntegralScreen::is_screened(int i, int j){
  int idx = _idx->get_index(i,j);
  if (_prescreened.count(idx) == 0) return false;
  return true;
}
   
    
/* @brief Prints out stats on the number of integrals which weren't
 *        calculated due to CS screening */
void IntegralScreen::print_screen_stats(){
  //    cout << "Integrals screened by CS = " << _nScreen << endl;
}

int IntegralScreen::prescreen(){
    double max_diag = *max_element(_diag.begin(), _diag.end());
    double test = (_cd_thresh * _cd_thresh) / max_diag;
    int nBas = _basis_data->get_nBas();

    /* Now we loop over all the diagonal elements and figure out which (if any)
     * we can safely ignore
     * iDiag points to the diagonal elements 
     * ij    keeps track of how many diagonal elements we have kept so far
     * If an element passes screening we add it to the pairMap
     * otherwise we add it to prescreened
     */
    int ij = 0;
    for (int iBas = 0; iBas < nBas; iBas++){
        for (int jBas = 0; jBas <= iBas; jBas++){
            int iDiag = diag_address(iBas, jBas);
            if (_diag[iDiag] > test)
                _pair_map->set_pair(iBas, jBas, ij++);
            else{
	        int idx = _idx->get_index(iBas, jBas);
                _prescreened.insert(idx);
	    }
        }
    }
    
    cout << " Prescreening threshold      = " << test << endl;
    cout << " Basis pairs                 = " << nBas * (nBas + 1) / 2 << endl;
    cout << " Basis pairs after screening = " << ij << endl;
    cout << endl;
    
    return ij;
}

/* @brief Calculates the diagonal integrals (ij|ij)
 */
void IntegralScreen::calc_diag_ints(){
    vector<int> nFunc = _basis_data->get_nFuncInShell();
    int nShell = _basis_data->get_nShell();
    int nBas = _basis_data->get_nBas();
    _diag.resize(nBas * nBas);
    
    for (int iShell = 0; iShell < nShell; iShell++){
        for (int jShell = 0; jShell <= iShell; jShell++){
            
            //reorder shells for Erkale/libint if necessary
            int erk_iShell = iShell;
            int erk_jShell = jShell;
            
            int am_i = funcs_to_am(nFunc[erk_iShell]);
            int am_j = funcs_to_am(nFunc[erk_jShell]);
            
            if (am_i < am_j){
                int tmp = erk_iShell;
                erk_iShell = erk_jShell;
                erk_jShell = tmp;
            }
            
            vector<double> eris = _basis_data->eval_ijkl(erk_iShell, erk_jShell,
                                                   erk_iShell, erk_jShell);
            
            //We need to place the integrals into the diagonal vector
            //We do this according to Jeremy's original formulas
            for (int iBas = 0; iBas < nFunc[erk_iShell]; iBas++){
                int jLimit = nFunc[erk_jShell] - 1;
                if (erk_iShell == erk_jShell) jLimit = iBas;
                for (int jBas = 0; jBas <= jLimit; jBas++){
                    int i = accumulate(nFunc.begin(), nFunc.begin() + erk_iShell, 0) + iBas;
                    int j = accumulate(nFunc.begin(), nFunc.begin() + erk_jShell, 0) + jBas;
                    
                    int ij = diag_address(i, j);
                    int ijij = ((iBas * nFunc[erk_jShell] + jBas)
                               * nFunc[erk_iShell] + iBas)
                               * nFunc[erk_jShell] + jBas;

                    _diag[ij] = eris[ijij];
                }
            }            
        }
    }
}
