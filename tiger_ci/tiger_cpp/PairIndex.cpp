/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

/* 
 * @file   PairIndex.cpp
 * @author David Krisiloff and Francis Ricci
 * 
 * Created on August 4, 2014, 2:35 PM
 */

#include "PairIndex.h"

PairIndex::PairIndex(size_t max_dim_1, size_t max_dim_2, bool symmetric): sym(symmetric){
    index.set_size(max_dim_1, max_dim_2);

    if (!sym){
        // For the non-symmetric case we'll need to enumerate the pairs
        // and store the addresses
        index.resize(max_dim_1, max_dim_2);
        int count = 0;
        for (size_t i = 0; i< max_dim_2; ++i){
            for (size_t j=0; j< max_dim_1; ++j){
                index(j,i) = count;
                count++;
            }
        }
    }
}

size_t PairIndex::get_index(size_t i, size_t j){
    if (sym){
        if (i<j) std::swap(i,j);
        /*
         * We use the formula out of erkale. Our pairs are stored 
         * as 
         * (0,0)
         * (0,1)
         * (1,1)
         * (1,2)
         * (2,2)
         * (0,3)
         * 
         * We calculate the start of the ith block first then add j 
         */
        size_t block_start = i*(i+1)/2 ;
        return block_start + j;
    }
    return index(i,j);
}

RestrictedIndex::RestrictedIndex(size_t max_dim_1, size_t max_dim_2, bool symmetric):
  PairIndex(max_dim_1, max_dim_2, symmetric){
  _pair2idx.resize(max_dim_1 * max_dim_2);
  _idx2pair.resize(max_dim_1 * max_dim_2);
}

int RestrictedIndex::get_index(int i, int j){
  //  if (i > j) swap(i,j);
  return _pair2idx[PairIndex::get_index(i,j)];
}

std::pair<int, int> RestrictedIndex::get_pair(int idx){
  return _idx2pair[idx];
}

void RestrictedIndex::set_pair(int i, int j, int ij){
    _pair2idx[PairIndex::get_index(i,j)] = ij;
    std::pair<int, int> ij_pair(i,j);
    _idx2pair[ij] = ij_pair;
}
