/*
 *  @file ErkaleDF.h
 *
 *  Created on: Jan 24, 2015
 *      Author: Francis Ricci
 */

#ifndef TIGER_CPP_ERKALEDF_H_
#define TIGER_CPP_ERKALEDF_H_

#include <erkale/basis.h>

#include <tiger_cpp/PairIndex.h>

#include <memory>

/**
 * @class ErkaleDF
 * @brief Implementation of necessary density fitting data, using Erkale
 */

class ErkaleDF {
public:
	ErkaleDF(std::shared_ptr<BasisSet> bas, std::shared_ptr<BasisSet> auxbas);

    /**
     * @brief Calculates the three center integrals \f$(A|\mu \nu)\f$ between
     *        the auxiliary basis and the normal basis
     * @param integrals Matrix that stores the integrals on exit
     * @param addr An index for the addresses of the \f$(\mu \nu)\f$ pairs
     */
    void compute_three_center_integrals(arma::mat& integrals, PairIndex& addr);

    /**
     * @brief Calculates the two center integrals \f$(A|B)\f$ involving just
     *        functions from the auxiliary basis set
     * @param ab The final \f$(A|B)\f$ matrix
     */
    void compute_two_center_integrals(arma::mat& ab);

    //return number of auxiliary basis functions
    int get_nAux();

private:
    std::shared_ptr<BasisSet> _basis;
    std::shared_ptr<BasisSet> _auxbasis;
};

#endif /* TIGER_CPP_ERKALEDF_H_ */
