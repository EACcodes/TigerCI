// c_sort.h
//
// A quick interface to the STL to lets us sort vectors without having to write
// the sorting routines ourselves
//
// All sorts sort into ascending order
#ifndef C_SORT_DEF
#define C_SORT_DEF

#include <inttypes.h>
#include <algorithm>
using namespace std;

// the following routines are being directly accessed 
// from fortran and we need to avoid mangling the names
// **************************************************
extern "C" 
{
    // sort an array of doubles of size n 
    void c_sort_real_array(double *vector, const int64_t& n);
    
    // sort an array of integers of size n 
    void c_sort_int_array(int64_t *vector, const int64_t& n);
    
    // sorts a size n vector of doubles and it's associated index vector
    void c_sort_real_array_with_index(double *vector, int64_t *index, const int64_t& n);
    
    // sorts a size n vector of doubles and it's associated index vector
    void c_sort_int_array_with_index(int64_t *vector, int64_t *index, const int64_t& n);

  // sorts a size n vector doubles and then applies that permutation to an arbitrary vector
  void c_sort_real_array_with_arb_permute(double *vector, int64_t *avector, const int64_t& n);
}

#endif


