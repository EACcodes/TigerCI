#ifndef PAIRINDEXWRITER_H
#define	PAIRINDEXWRITER_H

#include <memory>
#include "ConvertKeys.h"
#include "fortran_api.hpp"
#include "PairIndex.h"
#define FINT long long 

/**
 * @class PairIndexWriter
 * @brief Handles the writing Pair Indeces to the I/O buffer
 *        for use in the Fortran parts of TigerCI.
 * 
 * This is just a simple series of wrappers over the I/O buffer functionality.
 * The I/O buffer is extremely useful, not just for the mem/disk buffering, but 
 * also for providing a simple way of keeping this data available for later 
 * Fortran code.
 * 
 * @todo Probably want to make the amount of memory used in the buffer related to "CD MEMORY"
 */
class PairIndexWriter {
public:
    
    /*
     * @brief The PairIndexWriter handles writing CD results to the I/O buffer
     * @param file_number  File number (fortran style) for the I/O buffer. We always use 305
     * @param nthread      Maximum number of threads that will read/write to this I/O buffer (in C++ or Fortran!)
     * @param workdir      Location of the work directory
     * @param filename     Name of output file
     */
    PairIndexWriter(const size_t file_number, const size_t nthread, const std::string& workdir, const std::string& filename);
    
    /*
     * @brief Given a pair index quantity \f$ij\f$ writes out the pair indices i and j to the I/O buffer 
     * @param ind  The index (ij)
     * @param a  The vector containing data to be written
     */
    void write_pair(const size_t ind, const std::vector<size_t>& pair);
    
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
    const size_t unit_number;
    FINT pool_id;
};


#endif	/* PAIRINDEXWRITER_H */

