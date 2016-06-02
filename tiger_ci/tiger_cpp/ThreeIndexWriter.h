#ifndef THREEINDEXWRITER_H
#define	THREEINDEXWRITER_H

#include <memory>
#include "ConvertKeys.h"
#include "fortran_api.hpp"
#include "PairIndex.h"
#define FINT long long 

/**
 * @class ThreeIndexWriter
 * @brief Handles the writing three index quantities (DF/CD) to the I/O buffer
 *        for use in the Fortran parts of TigerCI.
 * 
 * This is just a simple series of wrappers over the I/O buffer functionality.
 * The I/O buffer is extremely useful, not just for the mem/disk buffering, but 
 * also for providing a simple way of keeping this data available for later 
 * Fortran code.
 * 
 * @todo Probably want to make the amount of memory used in the buffer related to "CD MEMORY"
 */
class ThreeIndexWriter {
public:
    
    /*
     * @brief The ThreeIndexWriter handles writing DF or CD results to the I/O buffer
     * @param file_number  File number (fortran style) for the I/O buffer. We always use 305
     * @param nAux         Number of auxiliary basis functions or cholesky vectors
     * @param nOrb         Number of molecular orbitals
     * @param nthread      Maximum number of threads that will read/write to this I/O buffer (in C++ or Fortran!)
     * @param workdir      Location of the work directory
     */
    ThreeIndexWriter(const size_t file_number, const size_t nAux, const size_t nOrb, const size_t nthread, const std::string& workdir);
    
    /*
     * @brief Given a three index quantity \f$A^P_ij\f$ writes out the vector \f$A^(:)_ij\f$ to the I/O buffer (i.e.
     *        we write out the vector for a fixed (ij) pair)
     * @param i  The index i of (ij)
     * @param j  The index j of (ij)
     * @param a  The vector containing data to be written
     */
    void write_array(const size_t i, const size_t j, const std::vector<double>& a);
    
    /*
     * @return Returns the file number for the file containing 3 index quantities
     */
    size_t get_file_number();
    
    /*
     *@return Returns the I/O buffer pool number for the buffer containing 3 index quantities
     */
    size_t get_pool_number();
    
    /*
     *@brief Correctly sets a few global variables pertaining to the 3 index quantities. 
     * 
     * This functions must be run before moving into the Fortran portion of TigerCI!
     */
    void communicate_size_data(ConvertKeys global);
    
private:
    const size_t unit_number, naux, norb;
    PairIndex index;
    FINT pool_id;
};


#endif	/* THREEINDEXWRITER_H */

