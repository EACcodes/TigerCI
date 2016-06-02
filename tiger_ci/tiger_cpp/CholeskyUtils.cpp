/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <tiger_cpp/CholeskyUtils.h>
#include <tiger_cpp/FunctionMap.h>

#include <set>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

using namespace std;

PivotHandler::PivotHandler(int N) {
    last_i = 0;
    last_j = 0;

    for (int i = 0; i < N; ++i)
        P.push_back(i);
}

void PivotHandler::pivot(int i, int j) {
    int tmp;
    tmp = P[i];
    P[i] = P[j];
    P[j] = tmp;
    last_i = i;
    last_j = j;
}

vector<int> PivotHandler::get_pivots() const {
    vector<int> new_P = P;
    return new_P;
}

void PivotHandler::pivot_vector(const vector<double>& in, vector<double>& out) {
    for (size_t i = 0; i < in.size(); ++i) {
        out.at(i) = in.at(P.at(i));
    }
}

int PivotHandler::get_curr_val(const int i) const {
    return P.at(i);
}

int PivotHandler::get_position_of(const int i) const {
    auto it = find(P.begin(), P.end(), i); // find the location of i in P
    int p = distance(P.begin(), it);
    return p;
}

void PivotHandler::update_matrix(arma::mat& A) const {
    // this may not be the fastest way to implement this ...
    A.swap_rows(last_i, last_j);
}

void PivotHandler::update_vector(vector<double>& v) const {
    double tmp = v.at(last_i);
    v.at(last_i) = v.at(last_j);
    v.at(last_j) = tmp;
}

ChoShellCol::ChoShellCol(BasisData *bas, const int activePairs,
        RestrictedIndex& map, IntegralScreen& iscreen) :
pairmap(map), screen(iscreen), _basis_data(bas) {
    n = activePairs;
    int nshells = bas->get_nShell();
    vector<int> func_in_shell = bas->get_nFuncInShell();

    for (int shell = 0; shell < nshells; ++shell)
        for (int bfunc = 0; bfunc < func_in_shell.at(shell); ++bfunc)
            basis_2_shell.push_back(shell);

    for (int i = 0; i < n; ++i) {
        buffer_index.push_back(-1);
    }

    int max_func_in_shell = *max_element(func_in_shell.begin(), func_in_shell.end());
    max_size = max_func_in_shell * max_func_in_shell;

    for (int i = 0; i < max_size; ++i) {
        buffer.push_back(vector<double>());
        for (int j = 0; j < n; ++j) {
            buffer.at(i).push_back(0.0);
        }
    }

    curr_buff_size = 0;
}

void ChoShellCol::clean_data_structs() {
    for (int i = 0; i < n; ++i)
        buffer_index.at(i) = -1;

    ij_in_buffer.clear();
    curr_buff_size = 0;
}

void ChoShellCol::setup_for_columns_including(const int ibas, const int jbas) {
    curr_ishell = basis_2_shell.at(ibas);
    curr_jshell = basis_2_shell.at(jbas);

    clean_data_structs();

    set<int> pairs;
    vector<int> func_in_shell = _basis_data->get_nFuncInShell();
    int i_start = accumulate(func_in_shell.begin(), func_in_shell.begin() + curr_ishell, 0);
    int j_start = accumulate(func_in_shell.begin(), func_in_shell.begin() + curr_jshell, 0);
    for (int i = 0; i < func_in_shell.at(curr_ishell); ++i) {
        for (int j = 0; j < func_in_shell.at(curr_jshell); ++j) {
            int ibas = i_start + i;
            int jbas = j_start + j;
            if (screen.is_screened(ibas, jbas))
                continue;
            int ij = pairmap.get_index(ibas, jbas);
            pairs.insert(ij);
        }
    }

    curr_buff_size = 0;
    for (auto ij : pairs) {
        buffer_index.at(ij) = curr_buff_size;
        ij_in_buffer.push_back(ij);
        curr_buff_size++;
    }
}

void ChoShellCol::print_buffer_status() {
    cout << endl;
    cout << "     Buffer Status" << endl;
    cout << "     -------------" << endl;
    cout << "     max buffer size = " << max_size << endl;
    cout << "     current buffer size = " << curr_buff_size << endl;
    cout << "     current ij in the buffer: ";
    for (auto x : ij_in_buffer)
        cout << " " << x;
    cout << endl;
    cout << "     ij to buffer index mapping ..." << endl;
    for (size_t i = 0; i < ij_in_buffer.size(); ++i) {
        int ij = ij_in_buffer.at(i);
        cout << "     ij:" << ij_in_buffer.at(i) << " -> ibuff:" << buffer_index.at(ij) << endl;
    }
}

void ChoShellCol::calculate_columns_including(const int p, const PivotHandler& pivots) {

    // get the current ij pair at the position p
    int ij = pivots.get_curr_val(p);
    pair<int, int> x = pairmap.get_pair(ij);
    const int ibas = x.first;
    const int jbas = x.second;
    vector<int> func_in_shell = _basis_data->get_nFuncInShell();
    int nshell = _basis_data->get_nShell();
    setup_for_columns_including(ibas, jbas);

#ifdef _OPENMP
#pragma omp parallel default(none) shared(nshell, func_in_shell)
	{
#endif
		vector<double> eris;
		BasisData *data = _basis_data->duplicate();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
		for (int kshell = 0; kshell < nshell; kshell++) {
            for (int lshell = 0; lshell <= kshell; lshell++) {

                if (screen.is_screened(curr_ishell, curr_jshell, kshell, lshell))
                    continue;

                //reorder shells for proper erkale am ordering
                int erk_iShell = curr_ishell;
                int erk_jShell = curr_jshell;
                int erk_kShell = kshell;
                int erk_lShell = lshell;

                int am_i = funcs_to_am(func_in_shell.at(erk_iShell));
                int am_j = funcs_to_am(func_in_shell.at(erk_jShell));
                int am_k = funcs_to_am(func_in_shell.at(erk_kShell));
                int am_l = funcs_to_am(func_in_shell.at(erk_lShell));

                if (am_i < am_j)
                    swap(erk_iShell, erk_jShell);

                if (am_k < am_l)
                    swap(erk_kShell, erk_lShell);

                // evaluate integrals
                eris = data->eval_ijkl(erk_iShell, erk_jShell, erk_kShell, erk_lShell);

                // the following are the number of basis functions 
                // before this shell, you can consider them conceptually
                // to be "pointers" to the current position in the
                // basis set
                vector<int>::iterator start = func_in_shell.begin();
                int ip = accumulate(start, start + erk_iShell, 0);
                int jp = accumulate(start, start + erk_jShell, 0);
                int kp = accumulate(start, start + erk_kShell, 0);
                int lp = accumulate(start, start + erk_lShell, 0);

                //loop over all basis functions in the shells we evaluated
                for (int iFunc = 0; iFunc < func_in_shell.at(erk_iShell); iFunc++) {
                    int p_i = iFunc * func_in_shell.at(erk_jShell);
                    for (int jFunc = 0; jFunc < func_in_shell.at(erk_jShell); jFunc++) {
                        int p_ij = (p_i + jFunc) * func_in_shell.at(erk_kShell);
                        for (int kFunc = 0; kFunc < func_in_shell.at(erk_kShell); kFunc++) {
                            int p_ijk = (p_ij + kFunc) * func_in_shell.at(erk_lShell);
                            for (int lFunc = 0; lFunc < func_in_shell.at(erk_lShell); lFunc++) {

                                //Move the pointer to the next integral
                                int ijkl = p_ijk + lFunc;

                                //these are the absolute basis function numbers
                                //for this integral
                                int i = ip + iFunc;
                                int j = jp + jFunc;
                                int k = kp + kFunc;
                                int l = lp + lFunc;

                                //avoid (ij) or (kl) if prescreened
                                if (screen.is_screened(i, j) || screen.is_screened(k, l))
                                    continue;

                                //Get the positions in the V matrix for
                                //(ij) and (kl), which is non-trivial
                                //due to the prescreening
                                int ij_idx = pairmap.get_index(i, j);
                                int kl_idx = pairmap.get_index(k, l);
                                int ij_buff = buffer_index.at(ij_idx);

                                (buffer.at(ij_buff)).at(kl_idx) = eris.at(ijkl);
                            }
                        }
                    }
                }
            }
        }
		delete data;

#ifdef _OPENMP
    } // close parallel region
#endif
    //  print_buffer_status();  
}

const vector<double>& ChoShellCol::get_column(const int p, const PivotHandler& pivots) {
    int ij = pivots.get_curr_val(p);
    int buffer_addr = buffer_index.at(ij);
    if (ij < 0) {
        cout << "Was asked for column " << p << " which, given the current pivoting, is actually " << ij << endl;
        cout << "which isn't actually in the current batch of columns. Fatal Error." << endl;
        cout << endl;
        exit(1);
    }
    return buffer.at(buffer_addr);
}

vector<int> ChoShellCol::get_ij_in_buffer() const {
    vector<int> vals = ij_in_buffer;
    return vals;
}
