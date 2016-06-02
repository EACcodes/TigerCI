/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * @file   Cholesky.cpp
 * Author: Francis Ricci
 *
 * Created on November 21, 2014, 1:47 PM
 */

#include <tiger_cpp/Cholesky.h>
#include <tiger_cpp/Timer.h>
#include <tiger_cpp/FunctionMap.h>
#include <tiger_cpp/c_io.h>
#include <tiger_cpp/Lapack_Defs.h>
#include <tiger_cpp/FunctionMap.h>
#include <tiger_cpp/ConvertKeys.h>

#include <new>
#include <numeric>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;
using namespace arma;

Cholesky::Cholesky(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars) :
				_nbf(bas->get_nBas()), _norb(orbitals.n_cols), MOs(orbitals) {
	_basis_data = bas;
	_cd_thresh = global_vars.get_global_double("cd_thresh");
	_mem = ((int64_t) global_vars.get_global_int("max_mem") * 1024 * 1024) / 8;
	_scratch_directory = global_vars.get_global_string("scratch_directory");
	_nthreads = global_vars.get_global_int("numThreads");
	_pair_map = new RestrictedIndex(bas->get_nBas(), bas->get_nBas(), true);
	double ao_thresh = global_vars.get_global_double("ao_integral_threshold");
	_screen = new IntegralScreen(ao_thresh, _cd_thresh, _basis_data, _pair_map);
	_globals = &global_vars;
}

Cholesky::~Cholesky() {
	delete _screen;
	delete _pair_map;
}

void Cholesky::prescreen() {
	_activePairs = _screen->prescreen();
}

int Cholesky::find_pair_after_pivot(const int ij) {
	auto it = find(P.begin(), P.end(), ij + 1); // find the current vector holding this pair P[it] = ipair
	auto p = distance(P.begin(), it);
	return p;
}

int Cholesky::get_pair_after_pivot(const int ij) {
	return find_pair_after_pivot(ij);
}



void Cholesky::debug_print_in_ascii(std::string filename) {
	ofstream myfile;
	myfile.open(filename);
	for (int i = 0; i < V.n_rows; ++i) {
		for (int j = 0; j < V.n_cols; ++j) {
			myfile << V(i, j) << " ";
		}
		myfile << endl;
	}
	myfile.close();
}

void Cholesky::transfrom_to_MO_basis() {

	// This might be a really, really bad way to do things. Just an FYI
	Timer time;
	cout << endl;
	cout << "*********************************" << endl;
	cout << "*                                       " << endl;
	cout << "*   AO->MO Cholesky Transform    " << endl;
	cout << "*                                       " << endl;
	cout << "*********************************" << endl;
	cout << endl;

	size_t norb = MOs.n_cols;
	size_t norb_pairs = norb * (norb + 1) / 2;
	// estimate memory and if it exceeds Cholesky mem switch to write L^A_munu
	// and switch to Fortran out-of-core
	//
	// L + V (V has two different sizes, take the max) + tmp
	size_t mem_usage = max((size_t) (_activePairs * _activePairs),
			_ncho * norb_pairs) + _ncho * _nbf * norb;
	double mem_usage_mb = mem_usage * 8 / 1024. / 1024.;
	cout << " Roughly max memory usage (MB) = " << mem_usage_mb << endl;

	if (mem_usage > _mem) {
		V.resize(_activePairs, _ncho);
		throw TransformationOutOfMem();	
	}

	// Step #1
	// Reorder our data in way which makes future computation really efficient
	// This is ugly because
	//   (1) We transpose the entire cholesky decomposition
	Timer t1;
	//size_t nbf_pairs = _nbf * (_nbf + 1) / 2;
	PairIndex ind(_nbf, _nbf, true);
	V.resize(_activePairs, _ncho);
	V = V.t();
	t1.print_elapsed_time(" step 1");

	Timer t2;
//	size_t norb = MOs.n_cols;
//	size_t norb_pairs = norb * (norb + 1) / 2;
	PairIndex map_munu(_nbf, _nbf, true);
	PairIndex map_ij(norb, norb, true);

	// Equation (1)
	// X^A_mui = \sum_nu c_inu * L^A_munu
	Cube<double> tmp(_ncho, _nbf, norb);
	// Note for this loop that we are pulling out columns from the pivoted matrix V
	// be careful to avoid columns (equivalently ij pairs) which were screened out
	// load balancing may be an issue due to the prescreening
	//auto start = P.begin();
	//auto end = P.end();
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tmp, map_munu, norb) schedule(dynamic)
#endif
	for (size_t i = 0; i < norb; i++) {
		for (size_t mu = 0; mu < _nbf; mu++) {
			tmp.slice(i).col(mu).zeros();
			for (size_t nu = 0; nu < _nbf; nu++) {
				if (_screen->is_screened(mu, nu)) continue;
				auto ipair = _pair_map->get_index(nu, mu); // this pair
				int pos = find_pair_after_pivot(ipair);
				tmp.slice(i).col(mu) += MOs(nu, i) * V.col(pos);
			}
		}
	}
	t2.print_elapsed_time(" step 2");

	// Equation (2)
	// Y^A_ji = \sum_mu c_jmu * X^A_mui
	Timer t3;
	V.set_size(_ncho, norb_pairs);
#ifdef _OPENMP
#pragma omp parallel for default(none) shared(tmp, map_ij, norb)
#endif
	for (size_t i = 0; i < norb; i++) {
		for (size_t j = 0; j <= i; j++) {
			size_t mo_index = map_ij.get_index(i, j);
			V.col(mo_index).zeros();
			for (size_t mu = 0; mu < _nbf; mu++) {
				V.col(mo_index) += MOs(mu, j) * tmp.slice(i).col(mu);
			}
		}
	}
	_norb = norb;
	t3.print_elapsed_time(" step 3");
	time.print_elapsed_time(" AO->MO Cholesky Transform B");
}

void Cholesky::get_CD_block(const size_t i, const size_t j, std::vector<double>& values) const {
	PairIndex mo_map(_norb, _norb, true);
	size_t addr = mo_map.get_index(i, j);
	assert(_ncho == values.size());
	for (size_t a = 0; a < _ncho; a++) {
		values[a] = V(a, addr);
	}
}

// Build the exact integrals and then check how well we approximate them
void Cholesky::check_against_analytic() {
	cout << "final pivots are " << endl;
	for (auto x : P)
		cout << x << endl;
	cout << endl;

	double x, y;
	int nbf = _basis_data->get_nBas();
	vector<int> nFunc = _basis_data->get_nFuncInShell();
	for (int i = 0; i < nbf; i++) {
		int iShell = bas_to_shell(i, nFunc);
		int iIdx = bas_in_shell(i, nFunc);
		for (int j = 0; j <= i; j++) {
			int jShell = bas_to_shell(j, nFunc);
			int jIdx = bas_in_shell(j, nFunc);
			for (int k = 0; k < nbf; k++) {
				int kShell = bas_to_shell(k, nFunc);
				int kIdx = bas_in_shell(k, nFunc);
				for (int l = 0; l <= k; l++) {
					int lShell = bas_to_shell(l, nFunc);
					int lIdx = bas_in_shell(l, nFunc);

					mat eris = _basis_data->eval_ijkl(iShell, jShell, kShell, lShell);
					x = eris(iIdx * nFunc[iShell] + jIdx, kIdx * nFunc[kShell] + lIdx);
					y = 0.0;
					if ((!_screen->is_screened(i, j)) && (!_screen->is_screened(k, l))) {
						int ind_ij = _pair_map->get_index(i, j);
						int ind_kl = _pair_map->get_index(k, l);

						int ij = find_pair_after_pivot(ind_ij);
						int kl = find_pair_after_pivot(ind_kl);

						for (int a = 0; a < _ncho; a++) {
							y += V(ij, a) * V(kl, a);
						}
					} else {
						cout << "   next integral was prescreened " << endl;
					}
					// output here
					cout << endl;
					cout << "( " << i << " " << j << " | ";
					cout << k << " " << l << " )" << " ";
					cout << "exact = " << x << " CD = " << y << endl;
				}
			}
		}
	}
}


int Cholesky::get_n_activePairs() {
	return _activePairs;
}

void Cholesky::get_cho_vector(const size_t ind, std::vector<double>& values) {
	assert(_activePairs == values.size());
	for (size_t a=0; a<_activePairs; a++){
		values[a]=V(a,ind);
	}
}

std::vector<size_t> Cholesky::get_active_pair(size_t ij){
	std::vector<size_t> pair;
	pair.push_back(_pair_map->get_pair(ij).first);
	pair.push_back(_pair_map->get_pair(ij).second);
	return pair;
}

bool Cho_WholeMat::enough_mem() {
	// I need to able to store the whole _activePairs by _activePairs
	// matrix for this to work
	if (_activePairs * _activePairs < _mem) {
		return true;
	}
	return false;
}

Cho_WholeMat::Cho_WholeMat(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars) :
				Cholesky(bas, orbitals, global_vars) {
	cout << endl;
	cout << "*********************************" << endl;
	cout << "*                                " << endl;
	cout << "*   Cholesky decomposition       " << endl;
	cout << "*   Fully pivoted                " << endl;
	cout << "*                                " << endl;
	cout << "*********************************" << endl;
	cout << endl;
	prescreen();
}

void Cho_WholeMat::decompose() {
	/* This cholesky routine is designed to be as simple as possible
	 * We are going to build the two-electron integral matrix (or those
	 * pieces which survive screening) and then let LAPACK deal with
	 * doing the heavy math.
	 */

	Timer decomp_timer;
	if (!enough_mem())
	throw CholeskyOutOfMem();

	V = calculate_eris();

	//call lapack to do the cholesky factorization, requires weird int size setup
	int_lapack n = (int_lapack) _activePairs;
	int_lapack piv[n];
	int_lapack rank = 0;
	int_lapack ierr = 0;

	char uplo = 'L';
	double scratch[n * 2];
	dpstf2_(&uplo, &n, V.memptr(), &n, piv, &rank, &_cd_thresh, scratch, &ierr);

	cout << " Cholesky decomposition threshold = " << _cd_thresh << endl;
	if (ierr < 0){
		ostringstream oss;
		oss << "DPSTRF reports error # " << ierr;
		throw runtime_error(oss.str());
	}
	int P2[n];
	for (int i = 0; i < n; i++) {
		P.push_back((int) piv[i]);
		P2[n] = (int) piv[i];
	}

	_ncho = (size_t) rank;
	cout << " Number of cholesky vectors = " << _ncho << endl;
	decomp_timer.print_elapsed_time(" total two-electron integral cholesky decomposition");
	cout << endl;
}

mat Cho_WholeMat::calculate_eris() {
	Timer eri_time;

	vector<int> nFunc = _basis_data->get_nFuncInShell();

	//zeroing is required if any integrals are screened
	mat V = zeros<mat>(_activePairs, _activePairs);

	/* Integral Evaluations
	 *************************/

	int nShell = _basis_data->get_nShell();

#ifdef _OPENMP
#pragma omp parallel default(none) shared(nShell, nFunc, V)
	{
#endif
		vector<double> eris;
		BasisData *data = _basis_data->duplicate();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
		for (int iShell = 0; iShell < nShell; iShell++) {
			for (int jShell = 0; jShell <= iShell; jShell++) {
				for (int kShell = 0; kShell < nShell; kShell++) {
					for (int lShell = 0; lShell <= kShell; lShell++) {

						int ijShell = iShell * (iShell + 1) / 2 + jShell;
						int klShell = kShell * (kShell + 1) / 2 + lShell;

						//By permutational symmetry I only need ijShell <= klShell
						if (ijShell > klShell || _screen->is_screened(iShell, jShell, kShell, lShell))
							continue;

						//reorder shells for proper erkale am ordering
						int erk_iShell = iShell;
						int erk_jShell = jShell;
						int erk_kShell = kShell;
						int erk_lShell = lShell;

						int am_i = funcs_to_am(nFunc[erk_iShell]);
						int am_j = funcs_to_am(nFunc[erk_jShell]);
						int am_k = funcs_to_am(nFunc[erk_kShell]);
						int am_l = funcs_to_am(nFunc[erk_lShell]);

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
						vector<int>::iterator start = nFunc.begin();
						int ip = accumulate(start, start + erk_iShell, 0);
						int jp = accumulate(start, start + erk_jShell, 0);
						int kp = accumulate(start, start + erk_kShell, 0);
						int lp = accumulate(start, start + erk_lShell, 0);

						//loop over all basis functions in the shells we evaluated
						for (int iFunc = 0; iFunc < nFunc[erk_iShell]; iFunc++) {
							int p_i = iFunc * nFunc[erk_jShell];
							for (int jFunc = 0; jFunc < nFunc[erk_jShell]; jFunc++) {
								int p_ij = (p_i + jFunc) * nFunc[erk_kShell];
								for (int kFunc = 0; kFunc < nFunc[erk_kShell]; kFunc++) {
									int p_ijk = (p_ij + kFunc) * nFunc[erk_lShell];
									for (int lFunc = 0; lFunc < nFunc[erk_lShell]; lFunc++) {

										//Move the pointer to the next integral
										int ijkl = p_ijk + lFunc;

										//these are the absolute basis function numbers
										//for this integral
										int i = ip + iFunc;
										int j = jp + jFunc;
										int k = kp + kFunc;
										int l = lp + lFunc;

										//avoid (ij) or (kl) if prescreened
										if (_screen->is_screened(i, j) || _screen->is_screened(k, l))
											continue;

										//Get the positions in the V matrix for
										//(ij) and (kl), which is non-trivial
										//due to the prescreening
										int ij_idx = _pair_map->get_index(i, j);
										int kl_idx = _pair_map->get_index(k, l);

										if (ij_idx < kl_idx)
											swap(ij_idx, kl_idx);
										V(ij_idx, kl_idx) = eris[ijkl];
									}
								}
							}
						}
					}
				}
			}
		}
		delete data;
#ifdef _OPENMP
	}
#endif
	//_screen->print_screen_stats();
	eri_time.print_elapsed_time(" (ij|kl) evaluations");

	return V;
}

Cho_Pedersen::Cho_Pedersen(BasisData *bas, const arma::mat& orbitals, const ConvertKeys& global_vars) :
				Cholesky(bas, orbitals, global_vars) {

	cout << endl;
	cout << "**************************************" << endl;
	cout << "*                                     " << endl;
	cout << "*   Cholesky decomposition            " << endl;
	cout << "*   Partial pivoting with shell blocks" << endl;
	cout << "*                                     " << endl;
	cout << "**************************************" << endl;
	cout << endl;
	prescreen();
}

bool Cho_Pedersen::enough_mem() {
	// I'm going to need somewhere between
	// _activePairs * _nbf
	// and
	// _activePairs * (_nbf)^2
	// memory. The variability comes from the fact that I don't
	// know how many CD vectors I'll have when I'm done

	if (_mem < _activePairs * _nbf) {
		cout << " Not enough memory for the Pedersen CD" << endl;
		return false;
	}
	if (_mem < _activePairs * _nbf * _nbf) {
		cout << " I might have enough memory for the Pedersen CD" << endl;
		return true;
	}
	return true;
}

void Cho_Pedersen::decompose() {
	Timer t;

	if (!enough_mem())
		throw CholeskyOutOfMem();

	// Initial memory allocation for V
	// We start with _nactivePairs by _nbf, which is a minimal lower bound.
	// Later on we will increase memory by doubling the number of columns
	// each time we run out
	//
	// For other comp sci problems the above strategy usually strikes a
	// decent balance between number of times you need to allocate more mem
	// and actual mem usage. We'll see if it works here
	//V.resize(_activePairs, _nbf);
	V.resize(_activePairs, _nbf);
	vector<double> diag;
	diagonal_integrals(diag);

	PivotHandler pivots(_activePairs);
	ChoShellCol shell_cols(_basis_data, _activePairs, *_pair_map, *_screen);
	double x_max = *(max_element(diag.begin(), diag.end()));

	vector<double> tmp;
	for (int i = 0; i < _activePairs; ++i)
		tmp.push_back(i);

	cout << " Cholesky decomposition threshold = " << _cd_thresh << endl;

	_ncho = -1;
	while (_ncho + 1 < _activePairs) {
		_ncho++;
		int p = find_pivot(diag, _ncho);

		if (is_converged(diag, p)) {
			_ncho--;
			break;
		}

		if (_ncho == V.n_cols)
			get_more_mem();

		// We need to decompose another column, calculate the batch of columns  containing p
		shell_cols.calculate_columns_including(p, pivots);

		while (true) {
			const vector<double>& col = shell_cols.get_column(p, pivots); // do this before pivoting to avoid pivoting twice
			pivots.pivot(p, _ncho);
			pivots.update_matrix(V);
			pivots.update_vector(diag);
			pivots.pivot_vector(col, tmp);

			decompose_column(_ncho, tmp, diag);
			update_diag(diag, _ncho);

			p = find_pivot_in_block(diag, shell_cols, pivots);

			if (p < 0)
				break; // we don't have any more columns left in shell_cols
			if (is_converged(diag, p))
				break; // we converged
			x_max = *(max_element(diag.begin() + _ncho + 1, diag.end()));
			if (diag.at(p) < (x_max / 1000))
				break; // Pedersen et al. condition to stop decomposing more columns in this shell pair

			// we are decomposing another column in the shell_cols
			_ncho++;
			if (_ncho == V.n_cols)
				get_more_mem();
		}
	}

	P = pivots.get_pivots();

	// this is a bit weird, but our other CD uses LAPACK which starts with 1 instead of 0
	// just add 1 to everything in P for compatability later on
	for (size_t i = 0; i < P.size(); ++i)
		P[i] += 1;

	cout << " Number of cholesky vectors = " << _ncho << endl;
	t.print_elapsed_time(" total two-electron integral cholesky decomposition");
}

void Cho_Pedersen::diagonal_integrals(vector<double>& integrals) {
	Timer t;
	vector<int> func_in_shell = _basis_data->get_nFuncInShell();
	int nshell = _basis_data->get_nShell();

	//zeroing is required if any integrals are screened
	for (int i = 0; i < _activePairs; ++i)
		integrals.push_back(0.0);

#ifdef _OPENMP
#pragma omp parallel default(none) shared(nshell, func_in_shell, integrals)
	{
#endif
		vector<double> eris;
		BasisData *data = _basis_data->duplicate();
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
		for (int iShell = 0; iShell < nshell; iShell++) {
			for (int jShell = 0; jShell <= iShell; jShell++) {
				if (_screen->is_screened(iShell, jShell, iShell, jShell))
					continue;

				int am_i = funcs_to_am(func_in_shell.at(iShell));
				int am_j = funcs_to_am(func_in_shell.at(jShell));

				// Need these if we swap. DON'T DELETE
				int ishell2 = iShell;
				int jshell2 = jShell;

				if (am_i < am_j)
					swap(ishell2, jshell2);

				// evaluate integrals
				eris = data->eval_ijkl(ishell2, jshell2, ishell2, jshell2);

				// the following are the number of basis functions
				// before this shell, you can consider them conceptually
				// to be "pointers" to the current position in the
				// basis set
				vector<int>::iterator start = func_in_shell.begin();
				int ip = accumulate(start, start + ishell2, 0);
				int jp = accumulate(start, start + jshell2, 0);

				//loop over all basis functions in the shells we evaluated

				for (int iFunc = 0; iFunc < func_in_shell.at(ishell2); iFunc++) {
					int p_i = iFunc * func_in_shell.at(jshell2);
					for (int jFunc = 0; jFunc < func_in_shell.at(jshell2); jFunc++) {
						int p_ij = (p_i + jFunc) * func_in_shell.at(ishell2);
						int ijkl = (p_ij + iFunc) * func_in_shell.at(jshell2) + jFunc;

						//these are the absolute basis function numbers
						//for this integral
						int i = ip + iFunc;
						int j = jp + jFunc;

						//avoid (ij) if prescreened
						if (_screen->is_screened(i, j))
							continue;

						//Get the positions in the V matrix for
						//(ij) and (kl), which is non-trivial
						//due to the prescreening
						int ij_idx = _pair_map->get_index(i, j);
						integrals.at(ij_idx) = eris.at(ijkl);
					}
				}
			}
		}
		delete data;
#ifdef _OPENMP
	} //close parallel region
#endif
	t.print_elapsed_time(" (ij|ij) evaluations");
}

int Cho_Pedersen::find_pivot(vector<double>& integrals, int start) {
	double largest = integrals.at(start);
	int pos = start;
	for (int i = start + 1; i < _activePairs; ++i) {
		if (integrals.at(i) > largest) {
			largest = integrals.at(i);
			pos = i;
		}
	}
	return pos;
}

int Cho_Pedersen::find_pivot_in_block(vector<double>& integrals, const ChoShellCol &shell_cols, const PivotHandler& pivots) {
	auto shell_pair_list = shell_cols.get_ij_in_buffer();
	double max_val = -10000;
	int pos = -1; // if we never find a new pivot in this block we return -1
	for (int ij : shell_pair_list) {
		//     cout << "*** DEBUG processing pair " << ij << endl;
		int ij_addr = pivots.get_position_of(ij);
		//     cout << "  which is located at " << ij_addr << endl;
		bool decomposed = (ij_addr < _ncho + 1);
		//     cout << "  which was already decomposed ? " << decomposed << endl;
		if (!decomposed) {
			double val = integrals.at(ij_addr);
			//         cout << "  which has a value of " << val << endl;
			if (max_val < val) {
				max_val = val;
				pos = ij_addr;
				//             cout << "  new max val! "<< endl;
			}
		}
	}
	return pos;
}

bool Cho_Pedersen::is_converged(vector<double>& diag, int p) {
	if (diag.at(p) < _cd_thresh) {
		return true;
	}
	return false;
}

void Cho_Pedersen::get_more_mem() {
	//  Our current matrix is MxN and we want to resize to Mx(2N))
	int m = _activePairs;
	int n = V.n_cols;
	int new_mem = m * 2 * n;
	int required = new_mem + m*n;

	// if I can't get everything I want, how many extra columns can I get?
	int extra = _mem - 2 *m*n;
	extra /= m;

	/*
    cout << " MEM DEBUG" << endl;
    cout << " available mem = " << _mem <<endl;
    cout << " used          = " << m*n << endl;
    cout << " required for expansion " << required << endl;
    cout << " extra " << extra << endl;
	 */
	if (required < _mem) {
		mat tmp;
		tmp.set_size(m, n);
		tmp = V;
		V.zeros(m, 2 * n);
		V(span::all, span(0, n - 1)) = tmp(span::all, span::all);
	} else {
		if (extra == 0) {
			//  We can't fit even one more column into memory
			cout << endl;
			cout << " Ran out of memory " << endl;
			cout << " Current matrix size: m = " << m << " n = " << n << endl;
			throw CholeskyOutOfMem();
		} else {
			// we can't double our space, but we can allocate a bit more
			mat tmp;
			cout << " Warning! The CD is using the last of its available memory " << endl;
			cout << " Current matrix size: m = " << m << " n = " << V.n_cols + extra << endl;
			tmp.set_size(m, n);
			tmp = V;
			V.zeros(m, V.n_cols + extra);
			V(span::all, span(0, n - 1)) = tmp(span::all, span::all);
		}
	}
}

void Cho_Pedersen::decompose_column(const int i, const vector<double>& column, const vector<double>& diag) {
	for (int j = i + 1; j < _activePairs; ++j)
		V(j, i) = column.at(j);

	V(i, i) = sqrt(diag.at(i));

	// Careful, FORTRAN math routines start with 1 not 0!
	char trans = 'N';
	int_lapack m = (int_lapack) _activePairs - i;
	int_lapack n = (int_lapack) i;
	double alpha = -1.0;
	int_lapack LDA = (int_lapack) _activePairs;
	int_lapack incx = (int_lapack) _activePairs;
	double beta = 1.0;
	int_lapack incy = 1.0;
	double *A = V.colptr(0) + i + 1;
	double *x = V.colptr(0) + i;
	double *y = V.colptr(i) + i + 1;

	dgemv_(&trans, &m, &n, &alpha, A, &LDA, x, &incx, &beta, y, &incy);


	for (int j = i + 1; j < _activePairs; ++j) {
		V(j, i) /= V(i, i);
	}

}

void Cho_Pedersen::update_diag(vector<double>& diag, const int cho) {
	for (int i = cho + 1; i < _activePairs; ++i)
		diag.at(i) -= V(i, cho) * V(i, cho);
}

void ChoModel::decompose() {
	// A model CD
	int n = _activePairs;
	V = calculate_eris();

	// include the upper triangular (makes pivoting easier for testing)
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < i; ++j) {
			V(j, i) = V(i, j);
		}
	}

	vector<double> diag;
	for (int i = 0; i < n; ++i) {
		P.push_back(i);
		diag.push_back(V(i, i));
	}

	cout << "starting model CD" << endl;
	cout << "CD threshold = " << _cd_thresh << endl;
	cout << "size of V(nxn), n = " << n << endl;
	cout << endl;

	// Pivoted Cholesky Decomposition
	for (int j = 0; j < n; ++j) {

		cout << "Building cholesky vector " << j << endl;

		// find the next pivot
		int p = find_pivot(diag, j);

		// stop if we are done
		if (diag.at(p) < _cd_thresh) {
			cout << "Hit CD threshold" << endl;
			break;
		}

		// Pivot everything
		int tmp = P.at(j);
		P.at(j) = P.at(p);
		P.at(p) = tmp;

		double d = diag.at(j);
		diag.at(j) = diag.at(p);
		diag.at(p) = d;

		V.swap_cols(j, p);
		V.swap_rows(j, p);

		cout << "Decomposing column " << j << " pivot is " << p << endl;
		cout << "after pivot" << endl;
		for (int i = 0; i < n; ++i)
			cout << "  " << V(i, j) << endl;
		//

		//        // This is a tricky point, we need to swap both columns p,j and rows p,j
		//        // but we are only storing the lower triangular part of our matrix.
		//
		//        {
		//            // CALL DSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
		//            auto tmp = V(j, span(0,j-1));
		//            V(j, span(0,j-1)) = V(p, span(0,j-1));
		//            V(p, span(0,j-1)) = tmp;
		//
		//            // CALL DSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
		//            auto tmp2 = V(span(p+1, n-1), j);
		//            V(span(p+1, n-1), j) = V(span(p+1, n-1), p);
		//            V(span(p+1, n-1), p) = tmp;
		//
		//            // CALL DSWAP( PVT-J-1, A( J+1, J ), 1, A( PVT, J+1 ), LDA )
		//            auto tmp3 = V(span(j+1, p), j);
		//            V(span(j+1, p), j) = V(p, span(j+1, p)).t();
		//            V(p, span(j+1, p)) = tmp3.t();
		//        }

		V(j, j) = sqrt(diag.at(j));
		for (int i = j + 1; i < n; ++i) {
			for (int k = 0; k < j; k++) {
				V(i, j) = V(i, j) - V(i, k) * V(j, k);
			}
			V(i, j) = V(i, j) / V(j, j);
		}

		// just for checking later on, zero out the upper diagonal
		for (int i = 0; i < j; ++i) {
			V(i, j) = 0.0;
		}

		_ncho = j;

		// update diagonals
		for (int i = j + 1; i < n; ++i) {
			diag.at(i) -= V(i, j) * V(i, j);
		}


		cout << "cholesky vector " << _ncho << endl;
		cout << V(span::all, j) << endl;

	}

	// this is a bit weird, but our other CD uses LAPACK which starts with 1 instead of 0
	// just add 1 to everything in P for compatability later on
	for (size_t i = 0; i < P.size(); ++i)
		P[i] += 1;

	cout << " final pivots " << endl;
	for (auto x : P)
		cout << x << endl;


	cout << endl;
	cout << "Number of cholesky vectors = " << _ncho << endl;
}

int ChoModel::find_pivot(std::vector<double>& diag, int p) {
	double largest = diag.at(p);
	int pos = p;
	for (int i = p + 1; i < _activePairs; ++i) {
		if (diag.at(i) > largest) {
			largest = diag.at(i);
			pos = i;
		}
	}
	return pos;
}
