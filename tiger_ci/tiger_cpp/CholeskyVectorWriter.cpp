/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "CholeskyVectorWriter.h"

CholeskyVectorWriter::CholeskyVectorWriter(const size_t file_number, const size_t nCho, const size_t activePairs,
                                                 const size_t nthread, const std::string& workdir, const std::string& filename):
                                  unit_number(file_number), nCho(nCho), activePairs(activePairs) {

    pool_id = -1000;
    FINT maxmem, policy, threads;
    maxmem = 10000; // just for now, should be an input value later
    threads = (FINT) nthread;
    policy = 0;
 
    std::string fname = workdir + filename;
    for_double_buf_construct(maxmem, activePairs, policy, threads, pool_id);
    for_double_buf_openfile(pool_id, unit_number, &fname[0], fname.length());
}

size_t CholeskyVectorWriter::get_pool_number(){ 
    return pool_id;
}

size_t CholeskyVectorWriter::get_file_number(){
    return unit_number;
}

void CholeskyVectorWriter::write_array(const size_t ind, const std::vector<double>& a){
    FINT pos = (FINT) ind;
    pos++;      // We are numbering like FORTRAN (starting with 1)
    for_double_buf_writeblock(unit_number, pos, &(a[0]), (FINT)1);
}

void CholeskyVectorWriter::communicate_size_data(ConvertKeys global){ 
	global.add_global_int("Density Fitting Unit Number", unit_number);
	global.add_global_int("Density Fitting Pool ID", pool_id);
	global.add_global_int("Density Fitting naux", nCho);
	global.add_global_int("Nr of active Pairs",activePairs);
	global.add_global_bool("CPP_DECOMPOSED_INTS", true);
}

