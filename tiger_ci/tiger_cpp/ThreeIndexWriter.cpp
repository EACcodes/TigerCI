/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include "ThreeIndexWriter.h"

ThreeIndexWriter::ThreeIndexWriter(const size_t file_number, const size_t nAux, const size_t nOrb,
                                                 const size_t nthread, const std::string& workdir):
                                  unit_number(file_number), naux(nAux), norb(nOrb), index(nOrb, nOrb, true){

    pool_id = -1000;
    FINT maxmem, policy, threads;
    maxmem = 10000; // just for now, should be an input value later
    threads = (FINT) nthread;
    policy = 0;
 
    std::string fname = workdir + "mo_int.dat";
    for_double_buf_construct(maxmem, naux, policy, threads, pool_id);
    for_double_buf_openfile(pool_id, unit_number, &fname[0], fname.length());
}

size_t ThreeIndexWriter::get_pool_number(){ 
    return pool_id;
}

size_t ThreeIndexWriter::get_file_number(){
    return unit_number;
}

void ThreeIndexWriter::write_array(const size_t i, const size_t j, const std::vector<double>& a){
    FINT pos = (FINT) index.get_index(i,j);
    pos++;      // We are numbering like FORTRAN (starting with 1)
    for_double_buf_writeblock(unit_number, pos, &(a[0]), (FINT)1);
}

void ThreeIndexWriter::communicate_size_data(ConvertKeys global){
    global.add_global_int("Density Fitting Unit Number", unit_number);
    global.add_global_int("Density Fitting Pool ID", pool_id);
    global.add_global_int("Density Fitting naux", naux);
    global.add_global_bool("CPP_DECOMPOSED_INTS", true);
}




//void DensityFittingIOHandler::write_three_center_by_aux(DensityFitting& dfit){
//    PairIndex map(norb, norb, true);
//    
//    std::vector<double> array;
//    for (int i=0; i<naux ; ++i)
//        array.push_back(0);
//    
//    for (size_t i=0; i<norb; i++){
//        for (size_t j=0; j<=i; j++){
//            dfit.get_DF_block(i,j,array);
//            FINT pos = (FINT) map.get_index(i,j);
//            pos ++; // We are numbering like FORTRAN (starting with 1)
//            for_double_buf_writeblock(unit_number, pos, &(array[0]), (FINT)1);
//        }
//    }
//
//}
