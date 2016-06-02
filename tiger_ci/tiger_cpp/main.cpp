/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
/*
 * File: main.cpp
 * Authors: Francis Ricci, David Krisiloff, Johannes M Dieterich
 *
 */

#include <iostream>
#include <string>
#include <stdexcept>

#include <tiger_cpp/InputReader.h>
#include <tiger_cpp/Keyword.h>
#include <tiger_cpp/ConvertKeys.h>
#include <tiger_cpp/OrbitalParsers.h>
#include <tiger_cpp/Timer.h>
#include <tiger_cpp/MemUtils.h>

#include <tiger_cpp/local_ortho_virt.h>
#include <tiger_cpp/DensityFitting.h>
#include <tiger_cpp/Cholesky.h>
#include <tiger_cpp/ThreeIndexWriter.h> 
#include <tiger_cpp/CholeskyVectorWriter.h>
#include <tiger_cpp/PairIndexWriter.h>

#include <tiger_cpp/BasisData.h>
#include <tiger_cpp/TigerStructs.h>

#ifdef _OPENMP
#include <omp.h>
#endif 

using namespace std;

// prevent C++ name-mangling
extern "C"
{
// Fortran subroutine for collection of global variables from ConvertKeys
void grab_globals();
}

void density_fitting(BasisData *basis, arma::mat& MOs);

void cholesky(BasisData *bas, arma::mat& MOs, ConvertKeys global_vars);

int main(int argc, char** argv) {
    
        if(argc == 1){
		cout << "ERROR: TigerCI must be called with an input file as first argument." << endl;
		return 1;
	}
        
	Timer total_time;
        
        cout << endl;
        cout << "********************************************************************************" << endl;
        cout << "********************************************************************************" << endl;
        cout << "**********                                                            **********" << endl;
	cout << "                                     TigerCI                                    " << endl;
        cout << "           A Local Multireference Symmetric Group Graphical Approach            " << endl;
        cout << "                 Singles and Doubles Configuration Interaction                  " << endl;
        cout << "                                Computer Program                                " << endl;
        cout << "                     based on the work by Duch and Karwoski                     " << endl;
        cout << endl;
        cout << "                                     AUTHORS                                    " << endl;
        cout << endl;
        cout << "        Derek Walter, Arun Venkathanathan, Jeremy Chwee, David Krisiloff,       " << endl;
        cout << "       Johannes Dieterich, Andrew Szilva, Francis Ricci, Florian Libisch,       " << endl;
        cout << "                                Caroline Krauter                                " << endl;
        cout << "**********                                                            **********" << endl;
	cout << "********************************************************************************" << endl;
        cout << "********************************************************************************" << endl;
	cout << endl;
        cout << "********************************************************************************" << endl;
        cout << "**********                                                            **********" << endl;
        cout << "                              RECOMMENDED CITATIONS                             " << endl;
        cout << endl;
        cout << "Basic TigerCI usage (any calculation type):" << endl;
        cout << " * W. Duch, J. Karwowski, Comput. Phys. Rep. , 2, 93 (1985)" << endl;
        cout << " * D. Walter and E. A. Carter, Chem. Phys. Lett., 346, 177 (2001)" << endl;
        cout << " * D. Walter, A. Szilva, K. Niedfeldt, E. A. Carter," << endl;
        cout << "         J. Chem. Phys., 117, 1982 (2002)" << endl;
        cout << " * D. Walter, A. Venkatnathan, E. A. Carter, J. Chem. Phys., 118, 8127 (2003)" << endl;
        cout << " * D. B. Krisiloff, J. M. Dieterich, F. Libisch, E. A. Carter" << endl;
        cout << "   Numerical Challenges in a Cholesky-Decomposed Local Correlation Quantum Chemistry Framework" << endl;
        cout << "   in: Mathematical and Computational Modeling, pp. 59-91," << endl;
        cout << "   Editor: R. Melnick, John Wiley & Sons, Inc., (2015)" << endl;
        cout << endl;
        cout << "Additional citations based on calculation type:" << endl;
        cout << endl;
        cout << "For parallel TigerCI calculations:" << endl;
        cout << " * J. M. Dieterich, D. B. Krisiloff, A. Gaenko, F. Libisch, T. Windus," << endl;
        cout << "         M. S. Gordon, E. A. Carter, Comput. Phys. Commun., 185, 3175 (2014)" << endl;
        cout << endl;
        cout << "For local MRSDCI calculations:" << endl;
        cout << " * T. S. Chwee, A. B. Szilva, R. Lindh, E. A. Carter," << endl;
        cout << "         J. Chem. Phys., 128, 224106 (2008)" << endl;
        cout << endl;
        cout << "For (local) MRACPF(2) calculations:" << endl;
        cout << " * A. Venkatnathan, A. B. Szilva, D. Walter, R. J. Gdanitz, E. A. Carter," << endl;
        cout << "         J. Chem. Phys., 120, 1693 (2004)" << endl;
        cout << " * D. B. Krisiloff and E. A. Carter, Phys. Chem. Chem. Phys., 14, 7710 (2012)" << endl;;
        cout << " * D. B. Krisiloff, V. B. Oyeyemi, F. Libisch, and E. A. Carter," << endl;
        cout << "         J. Chem. Phys., 140, 024102 (2014)" << endl;
        cout << endl;
        cout << "For local excited states calculations:" << endl;
        cout << " * T. S. Chwee and E. A. Carter, J. Chem. Theory Comput., 7, 103 (2011)" << endl;
        cout << endl;
        cout << "For Cholesky-decomposed calculations" << endl;
        cout << " * T. S. Chwee and E. A. Carter, J. Chem. Phys., 132, 074104 (2010)" << endl;
        cout << endl;
        cout << "For acCD calculations:" << endl;
        cout << " * T. S. Chwee and E. A. Carter, Molecular Physics, 108, 2519 (2010)" << endl;
        cout << endl;
        cout << "For density-fitting calculations:" << endl;
        cout << " * D. B. Krisiloff, C. M. Krauter, F. Ricci, E. A. Carter, submitted (2015)" << endl;
        cout << endl;
        cout << "For integral-direct calculations:" << endl;
        cout << " * J. M. Dieterich and E. A. Carter, Comp. Theor. Chem., 1051, 47 (2015)" << endl;
        cout << endl;
        cout << "**********                                                            **********" << endl;
        cout << "********************************************************************************" << endl;
        cout << endl;

	ConvertKeys global_vars;
	try {
		// read input file
		InputReader setup;
		setup.read_tiguar_input(argv[1]);
		setup.output_keywords();

		// convert keywords into global variables
		global_vars.keys_to_globals(setup);
	}
	catch (runtime_error &e){
		cerr << "FATAL ERROR! " << e.what() << endl;
		return tiger::INPUT_FAILURE;
	}

#ifdef _OPENMP
	// tell OMP how many threads to use
	omp_set_num_threads(global_vars.get_global_int("numThreads"));
	cout << endl << "Running with OpenMP using " << global_vars.get_global_int("numThreads");
	cout << " threads" << endl;
#endif

	BasisData *basis;
	try{
		basis = new ErkaleBasis(global_vars, tiger::full);
	}
	catch (runtime_error &e){
		delete basis;
		cerr << "FATAL ERROR! " << e.what() << endl;
		return tiger::BASIS_FAILURE;
	}

	bool spherical = global_vars.get_global_bool("spherical");
	MoldenOrbitalParser parser = MoldenOrbitalParser(spherical);
	try{
		parser.parse(global_vars.get_my_string("orb_file"), basis);
	}
	catch (runtime_error &e){
		delete basis;
		cerr << "FATAL ERROR! " << e.what() << endl;
		return tiger::MO_PARSE_FAILURE;
	}

	try{
		bool nonlocal = (global_vars.get_global_int("nonlocal_flag") == 1);
		if (!nonlocal){
			ortho_virt(basis);
			parser.writeMoldenFile(global_vars.get_my_string("orb_file"));
		}
		parser.check_norms(basis);
	}
	catch (runtime_error &e){
		delete basis;
		cerr << "FATAL ERROR! " << e.what() << endl;
		return tiger::ORBITAL_FAILURE;
	}

	try{
		arma::mat MOs = get_MOs();
		// Density fitting or Cholesky decomposition
		if (global_vars.get_global_bool("DENSITY FITTING"))
			density_fitting(basis, MOs);
		else
			cholesky(basis, MOs, global_vars);
	}
	catch (runtime_error &e){
		delete basis;
		cerr << "FATAL ERROR! " << e.what() << endl;
		if (global_vars.get_global_bool("DENSITY FITTING"))
			return tiger::DF_FAILURE;
		else
			return tiger::CHOLESKY_FAILURE;
	}


	// Command line argument to debug (or time) only the C++ code
	if (argc > 2){
		string debug = string(argv[2]);
		if (debug == "DEBUG"){
			delete basis;
			total_time.print_elapsed_time("C++ EXECUTION TIME ");
			cout << endl;
			cout << "**********************************" << endl;
			cout << "Stopping early for debugging..." << endl;
			cout << "**********************************" << endl;
			cout << endl;
			return 0;
		}
	}

	// export global variables to fortran
	cout << endl;
	cout << "**********************" << endl;
	cout << "Entering Fortran..." << endl;
	cout << "**********************" << endl;
	cout << endl;
        
        cout << "Current memory usage is " << process_mem_usage() /1024. << " in MB" << endl;

	try {
		grab_globals();
	}
	catch (runtime_error &e){
		delete basis;
		cerr << "FATAL ERROR! " << e.what() << endl;
		return tiger::TIGER_FAILURE;
	}

	delete basis;
	// Just before exiting print total execution time
	total_time.print_elapsed_time("TOTAL EXECUTION TIME ");
	return tiger::SUCCESS;
}

void density_fitting(BasisData *basis, arma::mat& MOs){

	// Compute the density fitting
	DensityFitting dfit(basis, MOs);
	dfit.compute_fitting();

	// Write the density fitting to the IO Buffer for later use
	size_t nthreads = basis->_globals->get_global_int("numThreads");
	size_t naux = basis->get_nAux();
	size_t norb = MOs.n_cols;
	std::string workdir = basis->_globals->get_global_string("scratch_directory");
	ThreeIndexWriter io(305, naux, norb, nthreads, workdir);

	std::vector<double> array;
	for (int i=0; i<naux ; ++i)
		array.push_back(0);

	for (size_t i=0; i<norb; i++){
		for (size_t j=0; j<=i; j++){
			dfit.get_DF_block(i,j,array);
			io.write_array(i,j,array);
		}
	}
	basis->_globals->add_global_bool("CPP_TRANSFORMED_INTS", true);
	io.communicate_size_data(*basis->_globals);
}

void cholesky(BasisData *bas, arma::mat& MOs, ConvertKeys global_vars) {
	Cholesky* cho_impl;
	if (global_vars.get_global_bool("PEDERSEN CD"))
		cho_impl = new Cho_Pedersen(bas, MOs, global_vars);
	else
		cho_impl = new Cho_WholeMat(bas, MOs, global_vars);

	try {
		cho_impl->decompose();
		cho_impl->transfrom_to_MO_basis();
	}
	catch (CholeskyOutOfMem&) {
		// we ran out of memory somewhere
		// make sure the Fortran knows that it needs to do the CD
		cout << endl << "Not enough memory for the CD ... delaying to the Fortran out-of-core CD" << endl;
		global_vars.add_global_bool("CPP_DECOMPOSED_INTS", false);
		global_vars.add_global_bool("CPP_TRANSFORMED_INTS", false);
		return;
	}

	catch (TransformationOutOfMem&) {
		// we ran out of memory in the AO->MO transformation
		// make sure the Fortran knows that it needs to do the transformation but not the CD
		cout << endl << "Not enough memory for the AO-MO transformation ... delaying to the Fortran transformation" << endl;

		size_t ncho = cho_impl->get_ncho();
		size_t activePairs = cho_impl->get_n_activePairs();
		size_t nthreads = global_vars.get_global_int("numThreads");
		std::string workdir = global_vars.get_global_string("scratch_directory");
		std::string filename = "cvec.dat";

		// Write decomposed Cholesky Vectors to IO Buffer
		CholeskyVectorWriter io(302, ncho, activePairs, nthreads, workdir,filename);
		std::vector<double> array;
		for (size_t i = 0; i < activePairs; ++i)
			array.push_back(0);

		for (size_t i = 0; i < ncho; i++) {
			cho_impl->get_cho_vector(i, array);
			io.write_array(i, array);
		}
		io.communicate_size_data(global_vars);

		// Write Indices of active pairs to IO Buffer
		std::string filename_index = "cvec_index.dat";
		PairIndexWriter io_index(304,nthreads, workdir,filename_index);

		std::vector<size_t> ij(2,0);
		std::vector<std::vector<size_t>> vec_active_pairs(activePairs,ij); // initialize vector of activePairs (0,0) pairs

		for (size_t i=0; i<activePairs; i++) {
			ij = cho_impl->get_active_pair(i);        // get pair with index i
			size_t pos = cho_impl->get_pair_after_pivot(i); // where this pair is now located
			vec_active_pairs[pos]=ij;  // insert pair with index i at new position
		}

		for (size_t i = 0; i<activePairs;i++){
			ij = vec_active_pairs[i];
			io_index.write_pair(i,ij);
		}

		io_index.communicate_size_data(global_vars);

		delete cho_impl;

		global_vars.add_global_bool("CPP_TRANSFORMED_INTS", false);
		return;
	}  

	// at this point we know that the CD and the AO->MO transform succeeded
	// store the results and inform the Fortran code that it has decomposed integrals!
	size_t ncho = cho_impl->get_ncho();
	size_t nthreads = global_vars.get_global_int("numThreads");
	size_t norb = MOs.n_cols;
	std::string workdir = global_vars.get_global_string("scratch_directory");

	ThreeIndexWriter io(305, ncho, norb, nthreads, workdir);
	std::vector<double> array;
	for (int i = 0; i < ncho; ++i)
		array.push_back(0);

	for (size_t i = 0; i < norb; i++) {
		for (size_t j = 0; j <= i; j++) {
			cho_impl->get_CD_block(i, j, array);
			io.write_array(i, j, array);
		}
	}
	io.communicate_size_data(global_vars);
	delete cho_impl;
	global_vars.add_global_bool("CPP_TRANSFORMED_INTS", true);
	
}
