#ifndef CHOLESKYVECTORWRITER_H
#define	CHOLESKYVECTORWRITER_H

#include <memory>
#include "ConvertKeys.h"
#include "fortran_api.hpp"
#define FINT long long 

/**
 * @class CholeskyVectorWriter
 * @brief Handles the writing Cholesky vectors (usually before AO-MO tranbsformation) 
 *        to the I/O buffer for use in the Fortran parts of TigerCI.
 * 
 * This is just a simple series of wrappers over the I/O buffer functionality.
 * The I/O buffer is extremely useful, not just for the mem/disk buffering, but 
 * also for providing a simple way of keeping this data available for later 
 * Fortran code.
 * 
 * @todo Probably want to make the amount of memory used in the buffer related to "CD MEMORY"
 */
class CholeskyVectorWriter {
public:
    
    /*
     * @brief CholeskyVectorWriter handles writing of CD results before AO-MO transformation to the I/O buffer
     * @param file_number  File number (fortran style) for the I/O buffer. We always use 302
     * @param nCho         Number of cholesky vectors
     * @param activePairs  Number of active Pairs
     * @param nthread      Maximum number of threads that will read/write to this I/O buffer (in C++ or Fortran!)
     * @param workdir      Location of the work directory
     * @param filename     Name of output file
     */
    CholeskyVectorWriter(const size_t file_number, const size_t nCho, const size_t activePairs, const size_t nthread, const std::string& workdir, const std::string& filename);
    
    /*
     * @brief Given a matrix quantity \f$V_ij\f$ writes out the vector \f$V_j\f$ to the I/O buffer (i.e.
     *        we write out the vector for a fixed cholesky vector (i))
     * @param ind  The index ind 
     * @param a  The vector containing data to be written
     */
    void write_array(const size_t ind, const std::vector<double>& a);
    
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
    const size_t unit_number, nCho, activePairs;
    FINT pool_id;
};


#endif	/* CHOLESKYVECTORWRITER_H */

