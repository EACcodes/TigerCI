/* 
 * File:   Integral.h
 * Author: Francis Ricci
 *
 * Created on June 20, 2014, 9:32 AM
 */

#ifndef INTEGRAL_H
#define	INTEGRAL_H

#include <vector>
#include <armadillo>
#include <memory>

#include "basis.h"
#include "eriworker.h"

class ErkaleIntegral{
public:
	ErkaleIntegral(std::shared_ptr<BasisSet> bas);

	std::vector<double> eval_ijkl(int i, int j, int k, int l);
	arma::mat get_one_ints();
	arma::mat get_overlap();
	arma::mat get_overlap(ErkaleIntegral *bas);

private:    
	ERIWorker _worker;
	std::vector<GaussianShell> _shells;
	std::shared_ptr<BasisSet> _basis;
	GaussianShell *_iShell;
	GaussianShell *_jShell;
	GaussianShell *_kShell;
	GaussianShell *_lShell;
};

#endif	/* INTEGRAL_H */
