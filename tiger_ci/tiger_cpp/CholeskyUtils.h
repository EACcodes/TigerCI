// A bunch of utility classes for the cholesky decomposition

#ifndef CHOLESKYUTILS_H
#define	CHOLESKYUTILS_H

#include "BasisData.h"
#include "PairIndex.h"
#include "IntegralScreen.h"
#include <vector>
#include<armadillo>

/* @class PivotHandler
 * @brief Handles storing/applying pivoting during an incomplete
 *        Cholesky decomposition
 */
class PivotHandler{
public:
    
    /* @param N The size of the pivot array (number of rows in the matrix 
     *          being decomposed)
     */
    PivotHandler(int N);
    
    /* @brief Returns the current pivot array P
     * 
     * If there was a column/row in the original matrix, in the pivoted matrix
     * the column/row is now located at P[i]
     */
    std::vector<int> get_pivots() const;
    
    /* @brief Applies a new pivoting, switching i and j 
     * @param i First row/column to be pivoted
     * @param j Second row/column to be pivoted
     */
    void pivot(int i, int j);
    
    /* @brief Updates the matrix A(mxn) with the last pivot applied via the function "pivot"
     * @param A The matrix to be pivoted
     *
     * Note that we use this during the cholesky decomposition after finding a 
     * new pivot. For that use case, all this function needs to do is swap the
     * rows of A.
     */
    void update_matrix(arma::mat& A) const;
    
    /* @brief Updates the vector V with the last pivot applied via the function "pivot"
     * @param v The vector to be pivoted
     */   
    void update_vector(std::vector<double>& v) const;
    
    /* @brief Applies the pivot array P to the in vector and places the result in out
     * @param in The input vector
     * @param out The output vector
     */
    void pivot_vector(const std::vector<double>& in, std::vector<double>& out);
    
    int get_curr_val(const int i) const;
    
    
    int get_position_of(const int i) const;
        
private:
    int last_i;
    int last_j;
    std::vector<int> P;
};


/* @class ChoShellCol 
 * @brief Handles all the logic for calculating columns \f$(**|IJ)\f$ corresponding
 *        to the shell pair \f$IJ\f$ and for pivoting within that shell pair block
 */
class ChoShellCol {
public:
    ChoShellCol(BasisData *bas, const int activePairs,  RestrictedIndex& map,
                 IntegralScreen& iscreen);
    void calculate_columns_including(const int ij, const PivotHandler& pivots);
    const std::vector<double>& get_column(const int ij, const PivotHandler& pivots);
    std::vector<int> get_ij_in_buffer() const;
    
private:
    std::vector<int> ij_in_buffer;
    std::vector<int> buffer_index;
    std::vector<int> basis_2_shell;
    std::vector<std::vector<double>> buffer;
    
    int n;
    BasisData *_basis_data;
    int max_size;
    int curr_ishell, curr_jshell;
    int curr_buff_size;
    RestrictedIndex& pairmap;
    IntegralScreen& screen;
    
    void clean_data_structs();
    void setup_for_columns_including(const int ibas, const int jbas);
    void print_buffer_status();
};

#endif	/* CHOLESKYUTILS_H */

