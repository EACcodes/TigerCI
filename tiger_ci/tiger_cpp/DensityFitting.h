/**
 * @file DensityFitting.h
 * @brief Contains the TigerCI density fitting implementation
 * 
 */

#ifndef DF_H
#define	DF_H

#include <memory>
#include <cstdint>
#include <armadillo>

#include <IO_buffer/fortran_api.hpp>

#include <tiger_cpp/PairIndex.h>
#include <tiger_cpp/ConvertKeys.h>
#include <tiger_cpp/BasisData.h>
#include <tiger_cpp/FunctionMap.h>

/**
 * @class DensityFitting
 * @brief Performs a density fitting calculation
 *
 * This class calculates the density fitted (DF) two-electron integrals. There 
 * are two major differences compared to the DF in erkale. 
 * 1. We solve the expansion coefficient equations using a cholesky 
 *     decomposition instead of using an explicit inverse. This improves 
 *     numerical stability
 * 
 * 2. We write the final expansion coefficients using a symmetric
 *     contraction with \f$(A|B)^{1/2}\f$. This gives us a three index 
 *     quantity functionally identical to a Cholesky decomposition of the 
 *     two-electron integrals.
 * 
 * For details on this approach see: 
 * <a href=" http://dx.doi.org/10.1063/1.4820484">Epifanovsky, E., Zuev, D., Feng, X., Khistyaev, K., Shao, Y., & Krylov, A. I; J. Chem Phys 139 (2013)</a>
 */

class DensityFitting{
public:
    
    /**
     * @brief Constructs a density fitting calculation 
     * @param bas  The basis set
     * @param auxbas The auxiliary basis set
     * @param orbitals    The molecular orbitals
     */
    DensityFitting(BasisData *bas, const arma::mat& orbitals);
    
    /**
     * @brief Cleans up a density fitting calculation
     */
    virtual ~DensityFitting() {};
    
    /**
     * @brief Computes the density fitting 
     * 
     * Computes first the fitting in the AO basis then transforms into the MO 
     * basis
     */
    void compute_fitting();

    /**
     * @brief Copies the density fitted integrals (MO basis) for a particular
     *        (ij) pair 
     * @param i The i in (ij)
     * @param j The j in (ij)
     * @param values The memory where the density fitted integrals are stored. 
     *               We assume (and double check) that values already has the 
     *               correct size.
     */
    void get_DF_block(size_t i, size_t j, std::vector<double>& values) const;
    
    /**
     * 
     * @return the number of auxiliary functions
     */
    size_t get_naux() const;
    
    /**
     * 
     * @return the number of orbitals
     */
    size_t get_norb() const;


        
private:

    /// Normal basis set
    BasisData *basis;
    /// Number of functions in the normal basis set
    const size_t nbf;
    /// Number of functions in the auxiliary basis set
    const size_t naux;
    /// Number of molecular orbitals
    const size_t norb;
    
    /// Molecular orbitals
    const arma::mat MOs;
    
    /// auxiliary basis expansion coefficients (the final DF results))
    arma::mat coefficients;
    
    /**
     * @brief Calculates the cholesky decomposition of \f$(A|B)\f$
     * @param ab The matrix \f$(A|B)\f$
     * @param ab_sqrt One exit contains the cholesky decomposition, \f$(A|B)^{1/2}\f$
     */
    void compute_two_center_chol(const arma::mat& ab, arma::mat& ab_sqrt) const;

    /**
     * @brief Calculates the expansion coefficients \f$C^A{_\mu\nu}\f$ describing 
     *        how each \f$(\mu\nu)\f$ pair density is expanded in the auxiliary 
     *        basis set
     * @param ab is the matrix \f$(A|B)\f$
     */
    void calculate_aux_expansion_coeff(arma::mat& ab);

    void lapack_solve(arma::mat& ab);  

    /**
     * @brief Calculates the symmetric contraction of \f$C^A_{\mu\nu}\f$ with 
     *        \f$(A|B)^{1/2}\f$ yielding a result functionally equivalent to 
     *        a Cholesky decomposition of the two-electron integrals
     * @param ab_sqrt \f$(A|B)^{1/2}\f$ 
     */
    void compute_symmetric_contraction(const arma::mat& ab_sqrt);
    
    /**
     * @brief Transforms the density fitted integrals to the MO basis set
     */
    void transform_to_MO_basis();
    
    // The following are just debugging code and should never be run
    // in production code.
    void debug_MO_3cent_contracted(arma::mat& C, PairIndex& map);
    void debug_AO_3cent_contracted(arma::mat& B, PairIndex& map);
    void debug_AO_3cent(arma::mat& C, arma::mat& ab, PairIndex& map);
    void double_check(arma::mat& C);
};

#endif	/* DF_H */
