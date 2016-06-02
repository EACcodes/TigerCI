/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * File:   Integral.cpp
 * Author: Francis Ricci
 *
 * Created on June 20, 2014, 9:32 AM
 */

#include <stdexcept>

#include "ErkaleIntegral.h"

using namespace std;

ErkaleIntegral::ErkaleIntegral(std::shared_ptr<BasisSet> bas) :
			_worker(bas->get_max_am(), bas->get_max_Ncontr()){
	_basis = bas;
	_shells = _basis->get_shells();
	_iShell = NULL;
	_jShell = NULL;
	_kShell = NULL;
	_lShell = NULL;
}

vector<double> ErkaleIntegral::eval_ijkl(int i, int j, int k, int l){
	_iShell = &_shells[i];
	_jShell = &_shells[j];
	_kShell = &_shells[k];
	_lShell = &_shells[l];
	//calculate and store two electron integrals (ERIs)
	_worker.compute(_iShell, _jShell, _kShell, _lShell);

	return _worker.get();
}

arma::mat ErkaleIntegral::get_one_ints(){
	return _basis->kinetic() + _basis->nuclear();
}

arma::mat ErkaleIntegral::get_overlap(){
	return _basis->overlap();
}

arma::mat ErkaleIntegral::get_overlap(ErkaleIntegral *bas){
	return _basis->overlap(*bas->_basis);
}
