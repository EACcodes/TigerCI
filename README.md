# TigerCI

This is the stand-alone version of TigerCI, a local multi-reference singles doubles configuration
interaction code.

# Compilation
We require the following for compilation:
* cmake (version 2.5 or later)
* a Fortran compiler (tested: Intel Fortran and gfortran)
* a C++ compiler (tested: Intel CPP, g++, and clang++)
* a BLAS/LAPACK library (tested: Intel MKL and OpenBLAS)
* Armadillo
* GSL
* libint v1 API

The current code checkout contains parts of Erkale in the "erkale" subdirectory, which are under GPL2 or later.

To compile, please create a build directory and cmake in it. To adjust compile flags and/or provide paths other than the standard ones, use ccmake and/or -D macros to cmake.

# Use of TigerCI
TBD, please see documentation provided in documentation folder.
