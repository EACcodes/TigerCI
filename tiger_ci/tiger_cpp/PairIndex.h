/**
 * @file PairIndex.h
 * @brief Provides an index for pairs (i.e. (ij) orbital pairs)
 */

#ifndef PAIRINDEX_H
#define	PAIRINDEX_H

#include <armadillo>
#include <unordered_map>

/**
 * @class PairIndex
 * @brief Provides an index for pairs of object stored in an array
 */
class PairIndex {
public:
    /* 
     * Gives a 1D index for a pair (i,j). 
     * We assume the first dimension (max_dim_1)
     * changes the fastest
     * 
     * If pairs are symmetric (i,j) == (j,i) we
     * only store the unique pairs
     */
    
    /**
     * @brief Constructs a new pair index (ij)
     * @param max_dim_1 Maximum number elements for the first index (i)
     * @param max_dim_2 Maximum number of elements for the second index (j)
     * @param sym Is the pair symmetric (true if (ij) == (ji)) This results
     *            in just storing the unique pairs ( if (ij) == (ji) then both
     *            pairs have the same address)
     */
    PairIndex(size_t max_dim_1, size_t max_dim_2, bool sym);
    
    /**
     * @brief gives the index of pair (ij)
     * @param i i in (ij)
     * @param j j in (ij)
     * @return the index/address of (ij) 
     */
    size_t get_index(size_t i, size_t j);
private:
    const bool sym;
    arma::Mat<size_t> index;
};


class RestrictedIndex : public PairIndex{
public:
    RestrictedIndex(size_t max_dim_1, size_t max_dim_2, bool sym);
    int get_index(int i, int j);
    std::pair<int, int> get_pair(int idx);
    void set_pair(int i, int j, int ij);
private:
    std::vector<int> _pair2idx;
    std::vector<std::pair<int,int>> _idx2pair;
};

#endif	/* PAIRINDEX_H */

