/* 
 * @file   IntegralScreen.h
 * @author Francis Ricci
 *
 * Created on November 24, 2014, 11:28 AM
 */

#ifndef INTEGRALSCREEN_H
#define	INTEGRALSCREEN_H

#include <armadillo>
#include <vector>
#include <unordered_set>

#include "BasisData.h"
#include "PairIndex.h"

using namespace std;

class IntegralScreen{
public:
  IntegralScreen(double ao_thresh, double cd_thresh, BasisData *data,
		 RestrictedIndex *pair_map);
    ~IntegralScreen();
    
    /* @brief Answers the question is this integral shell block screened */
    bool is_screened(int iShell, int jShell, int kShell, int lShell);
    
    /* @brief Has this pair been prescreened out */
    bool is_screened(int i, int j);
    
    /* @brief Prints out stats on the number of integrals which weren't
     *        calculated due to CS screening */
    void print_screen_stats();
    
    int prescreen();
    
private:
    arma::mat _Q_int;
    arma::mat _Q_shell;
    double _ao_thresh;
    double _cd_thresh;
    int _nScreen;
    ErkaleIntegral *_ints;    
    BasisData *_basis_data;
    RestrictedIndex *_pair_map;
    PairIndex *_idx;
    vector<double> _diag;
    vector<int> _nFunc;
    unordered_set<int> _prescreened; //the (i,j) pairs removed during prescreen
    
    void setup_screening();
   
    void calc_diag_ints();
};

#endif	/* INTEGRALSCREEN_H */

