/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include "PairIndexWriter.h"

PairIndexWriter::PairIndexWriter(const size_t file_number, const size_t nthread, const std::string& workdir, const std::string& filename):
                                  unit_number(file_number){

    pool_id = -1000;
    FINT maxmem, policy, threads;
    maxmem = 10000; // just for now, should be an input value later
    threads = (FINT) nthread;
    policy = 0;
 
    std::string fname = workdir + filename;
    for_int_buf_construct(maxmem, (FINT) 2, policy, threads, pool_id);
    for_int_buf_openfile(pool_id, unit_number, &fname[0], fname.length());
}

size_t PairIndexWriter::get_pool_number(){ 
    return pool_id;
}

size_t PairIndexWriter::get_file_number(){
    return unit_number;
}

void PairIndexWriter::write_pair(const size_t ind, const std::vector<size_t>& pair){
    FINT pos = (FINT) ind;
    pos++;      // We are numbering like FORTRAN (starting with 1)
	std::vector<FINT> array;
	array.push_back((FINT) (pair[0]+1));
	array.push_back((FINT) (pair[1]+1));

    for_int_buf_writeblock(unit_number, pos, &(array[0]), (FINT)1);
}

void PairIndexWriter::communicate_size_data(ConvertKeys global){
    global.add_global_int("Pair Index Unit Number", unit_number);
    global.add_global_int("Pair Index Pool ID", pool_id);
}


