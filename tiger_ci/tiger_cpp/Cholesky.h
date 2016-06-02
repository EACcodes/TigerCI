/* 
 * @file   Cholesky.h
 * @author Francis Ricci
 *
 * Created on November 21, 2014, 1:47 PM
 * 
 *  For references on this work see:
 * (1) 1 H. Koch, A. Sánchez de Merás, and T.B. Pedersen, J. Chem. Phys. 118, 9481 (2003)
 * (2) 1 D.B. Krisiloff, J.M. Dieterich, F. Libisch, and E.A. Carter, in Prog. Math. Comput. Model. Nat. Soc. Sci., edited by R. Melnick (20XX)
 *      
 * Yes, the year on ref 2 is still unclear ....
 */

#ifndef CHOLESKY_H
#define	CHOLESKY_H

#include <utility>
#include <armadillo>
#include <map> //note that we can't use unordered maps for pairs unless we want
#include <set> //to define our own hash functions. This is almost as fast though
#include <vector>
#include <assert.h>
#include <exception>

#include "IntegralScreen.h"
#include "PairIndex.h"
#include "BasisData.h"
#include "ConvertKeys.h"
#include "CholeskyUtils.h"

using namespace std;

/*
 * @class CholeskyOutofMem
 * @brief Our custom exception for running out of memory during the CD
 */
class CholeskyOutOfMem : public bad_alloc{
    virtual const char* what() const throw() {
        return "Cholesky decomposition ran out of memory";
    }
};

/*
 * @class TranformationOutofMem
 * @brief Our custom exception for running out of memory during the Transformation to MO
 */
class TransformationOutOfMem : public bad_alloc{
    virtual const char* what() const throw() {
        return "AO->MO Transformation ran out of memory";
    }
};


/*
 *@class Cholesky
 *@brief This is the abstract base class that all of our Cholesky implementations
 *       inherit from.
 */
class Cholesky{
public:  
    /*
     * @brief Performs a Cholesky decomposition of the two-electron integrals
     * @param bas  The AO basis set
     * @param MOs  The molecular orbitals
     * @param global_vars  The global settings for this calculation
     * 
     * This routine performs the cholesky decomposition by constructing the 
     * two-electron integral matrix (subject to screening) and running the 
     * decomposition via LAPACK
     */
    Cholesky(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars);
    
    /*
     *@brief We have a non-default destructor because we currently have some mem on the heap.
     */
    virtual ~Cholesky();
    
    /*
     *@brief Runs the two-electron integral decomposition
     * 
     * This function is virtual and must be overriden in different Cholesky implementations
     * (i.e. how the decomposition is done distiniguishes different Cholesky codes)
     */
    virtual void decompose() = 0;
    
    /*
     *@return Returns the number of cholesky vectors
     */
    size_t get_ncho() const { return _ncho;};
    
    /*
     * @return Returns true if we have enough memory to run the CD in core. 
     *   
     * The exact logic varies on implementation
     */
    virtual bool enough_mem() = 0;
       
    /*
     * @brief Returns a block of MO CD quantities for a fixed i and j, \f$L^(:)_ij\f$
     * @param i The molecular orbital i
     * @param j The molecular orbital j 
     * @param values On exit contains the CD quantities. We assume 
     *               (and double check) that values already has the correct size.
     */
    void get_CD_block(const size_t i, const size_t j, vector<double>& values) const;
    
    /*
     * @brief Transforms the Cholesky decomposed integrals into the MO basis
     * 
     * This cannot be run before decompose().
     */
    void transfrom_to_MO_basis();
    
    /*
     * @brief Compute the two-electron integrals direct from libint and compares against
     *        integrals calculated from the cholesky vectors
     * 
     * Note that we must have already decomposed, but not transformed.
     */
    void check_against_analytic();
    
    /*
     * @brief Writes out the cholesky decomposition to a file in ascii
     * @param filename The name of the file we will write to
     * 
     * For debugging only!
     */
    void debug_print_in_ascii(std::string filename);
    
    /*
     * @brief Prints out a description of the method
     */
    virtual void print_header(){};
    
	/*
 	 * @brief Get number of active pairs
 	 */ 
	int get_n_activePairs();

	/*
 	 * @brief Return Cholesky vector i
 	 */
	void get_cho_vector(const size_t ind, std::vector<double>& values);

	/*
     * @brief Return active pair with index ij
     */
	std::vector<size_t> get_active_pair(size_t ij);  

    /*
     * @brief Returns the location of the pair ij after we have decomposed and pivoted the integrals
     * 
     */
    int get_pair_after_pivot(const int ij);

protected:
    BasisData *_basis_data;
    IntegralScreen *_screen;
    RestrictedIndex *_pair_map;
	const ConvertKeys *_globals;
    double _cd_thresh;
    int64_t _activePairs;
    int64_t _mem; //lets not overload the 32-bit ints :D
    string _scratch_directory;
	size_t _nthreads;
    arma::mat V;
    bool _prescreen_done;
    size_t _ncho;
    size_t _nbf;
    size_t _norb;
    const arma::mat& MOs;
    
    vector<int> P;    // I tend to get confused about the pivots. The easiest way to think about 
                      // P is that P[i] tells you the pair currently located at position i
                      // if you want to know where i is, call find_pair_after_pivot
    
    /*
     * @brief Prescreens the two-electron integrals using the formula provided 
     *        in the original Pedersen paper
     */    
    void prescreen();    

    /*
     * @brief Returns the location of the pair ij after we have decomposed and pivoted the integrals
     * 
     *  This operation is O(N). It's a bit tricky because the pivot array P tells us what is currently 
     *  located in position i (i.e. if P[i] == j then j is currently in position i). To find where ij 
     *  is we need to search all of P.
     */
    int find_pair_after_pivot(const int ij);
};


/*
 * @class Cho_WholeMat
 * @brief A cholesky implementation where we calculate (after screening) the whole 
 *        two-electron integral matrix and then decompose via LAPACK.
 * 
 * The pro's of this method are (1) we get a tighter decomposition and (2) we
 * aren't responsible for mainting the logic of the decomposition (LAPACK does 
 * it for us!). The con's are (1) memory usage is rather large (we need all of V)
 * and (2) we calculate more integrals than we need to (we need to fill all of V).
 */
class Cho_WholeMat: public Cholesky{
public:
    
    /*
     * @brief Performs a Cholesky decomposition of the two-electron integrals
     * @param bas  The AO basis set
     * @param MOs  The molecular orbitals
     * @param global_vars  The global settings for this calculation
     * 
     * This routine performs the cholesky decomposition by constructing the 
     * two-electron integral matrix (subject to screening) and running the 
     * decomposition via LAPACK
     */
    Cho_WholeMat(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars);
    
    /*
     *@brief Runs the two-electron integral decomposition
     */
    void decompose();    
    
    /*
     * @return Returns true if we have enough memory to run the CD in core. 
     *   
     * The exact logic varies on implementation
     */
    virtual bool enough_mem();
         
protected:
    
    /* @brief Calculates the lower triangular part of the two-electron integral matrix
     * @return arma::mat  The two-electron integral matrix 
     */
    arma::mat calculate_eris();
    
    /* @brief Wrapper for the LAPACK pivoted cholesky decomposition routine DPSTF2
     * @param A   The matrix to be decomposed
     * @return arma::mat  Permutation matrix used for the pivoted decomposition
     */
    arma::mat dpstf2_wrapper(arma::mat A);   
};


/* @ChoModel
 * @brief A model cholesky decomposition code. Used for testing/debugging. 
 */
class ChoModel : public Cho_WholeMat{
public: 
    ChoModel(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars):
                 Cho_WholeMat(bas, orbitals, global_vars){};
    void decompose();
    int find_pivot(std::vector<double>& diag, int p);
};

/*
 * @class Cho_Pedersen
 * @brief Cholesky decomposition implemented as in the Pedersen et al. paper
 * 
 * Calculates a CD of the two-electron integral matrix  by a modified pivoting
 * scheme which should prevent calculating extraneous two-electron integrals
 */
class Cho_Pedersen : public Cholesky{
public:
    /*
     * @brief Performs a Cholesky decomposition of the two-electron integrals
     * @param bas  The AO basis set
     * @param MOs  The molecular orbitals
     * @param global_vars  The global settings for this calculation
     * 
     * This routine performs the cholesky decomposition by constructing the 
     * two-electron integral matrix (subject to screening) and running the 
     * decomposition via LAPACK
     */
    Cho_Pedersen(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars);
                 
    /*
     *@brief Runs the two-electron integral decomposition
     */
    void decompose();      
    
    /*
     * @return Returns true if we have enough memory to run the CD in core. 
     *   
     * The exact logic varies on implementation
     */
    bool enough_mem();
    
    /*
     * @brief Calculates the diagonal integrals \f$(ij|ij)\f$
     * @param integrals Returns a vector of diagonal integrals
     */
    void diagonal_integrals(vector<double>& integrals);
    
    
private:
    
    /*
     * @brief Finds the next pivot (the next largest diagonal element)
     * @param integrals A list of all the diagonal integrals
     * @param start     The place to start searching for the next pivot (integral[start:])
     */
    int find_pivot(vector<double>& integrals, int start);
    
  
    /*
     * @brief Finds the next pivot within the already-calculated columns 
     * @param integrals The current diagonal integrals
     * @param shell_cols The calculated columns
     * @param pivots     The current pivots
     * 
     */
    int find_pivot_in_block(vector<double>& integrals, const ChoShellCol &shell_cols, const PivotHandler& pivots);
    
    /*
     * @brief Allocates more memory during the decomposition
     * 
     *  If we run out of memory, this will throw a bad_alloc. (It might be 
     *  cool to use that to hook into an out-of-core algorithm later on)
     */
    void get_more_mem();
    
    /*
     * @brief checks for convergence of the CD
     * @param diag The remaining diagonal integrals
     * @param p  The current pivot
     */
    bool is_converged(vector<double>& diag, int p);
    
    /* @brief Decomposes a single column of two-electron integrals
     * @param i The column to decompose
     * @param column The (already pivoted) column i
     * @param diag   The diagonal integrals
     */
    void decompose_column(const int i, const std::vector<double>& column, const std::vector<double>& diag);
    
    /* @brief Updates the diagonal integrals after decomposing column i
     * @param diag The diagonal integrals
     * @param i    The decomposed column
     */
    void update_diag(vector<double>& diag, const int i);    
};


#endif	/* CHOLESKY_H */

