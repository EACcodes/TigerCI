/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */
#include "c_sort.h"

void c_sort_real_array(double *vector, const int64_t& n){
    std::sort(vector, vector+n);
}

void c_sort_int_array(int64_t *vector, const int64_t& n){
    std::sort(vector, vector+n);
}

void c_sort_real_array_with_index(double *vector, int64_t *index, const int64_t& n){
    // This is a bit unusual. First we sort the index, based on the elements in vector
    // (noting that the index comes from Fortran and we want 0 based indexing in C++)
    std::sort(index, index+n, [&vector](int a, int b)->bool {
        return vector[a-1] < vector[b-1];});
        
    // Now that the index is sorted just sort the original vector
    std::sort(vector, vector+n);   
}

void c_sort_int_array_with_index(int64_t *vector, int64_t *index, const int64_t& n){
    // This is a bit unusual. First we sort the index, based on the elements in vector
    // (noting that the index comes from Fortran and we want 0 based indexing in C++)
    std::sort(index, index+n, [&vector](int a, int b)->bool {
        return vector[a-1] < vector[b-1];});
        
    // Now that the index is sorted just sort the original vector
    std::sort(vector, vector+n);   
}

