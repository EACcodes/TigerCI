// c_io.h
//
// A collection of routines to read/write big vectors to disk. Not in fortran because
// I/O in fortran is a pain (gfortran still uses a 32 bit interface for some I/O and until
// it implements that part of f2003 this is soooooo much easier)
//
#ifndef C_IO_DEF
#define C_IO_DEF

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <inttypes.h>
#include <unordered_map>
using namespace std;

// the following routines are being directly accessed 
// from fortran and we need to avoid mangling the names
// **************************************************
extern "C" 
{
  // opens a file for subsequent I/O
  void c_setup_vector_io(const char* filename);

  // closes all currently open files (if you open them again all the current info will be deleted)
  void c_close_files();

  // writes a vector to a file of vectors at a given location (i.e. write this vector at the position of the 4th vector)
  // numbering is base 0
  void c_write_vector(const char* filename, int64_t& position, double* vector, int64_t& length);

  // read a vector to a file of vectors at a given location (i.e. read a vector at the position of the 4th of vector)
  // numbering is base 0 
  void c_read_vector(const char* filename, int64_t& position, double* vector, int64_t& length);

  // writes a vector of ints to a file of vectors at a given location (i.e. write this vector at the position of the 4th vector)
  // numbering is base 0
  void c_write_ints(const char* filename, int64_t& position, int64_t* vector, int64_t& length);

  // read a vector of ints to a file of vectors at a given location (i.e. read a vector at the position of the 4th of vector)
  // numbering is base 0 
  void c_read_ints(const char* filename, int64_t& position, int64_t* vector, int64_t& length);

}

// remaining declarations ARE NOT called from fortran 
// **************************************************

// finds the file pointer corresponding to a file
FILE* get_file_pointer(const char* filename);

// checks for errors after a read/write operation
void error_processing(FILE* fp, const char* filename, const char* rw, size_t count, size_t length, size_t position);

// strips whitespace from the beginning/end of the filename
std::string file_name_clean(const char* filename);

#endif
