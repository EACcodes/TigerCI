cmake_minimum_required(VERSION 3.0)
project(tiger_ci C CXX Fortran)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

option(OMP "if OMP is on then the program is compiled with OpenMP support" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as no other type was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

#**************************
# Compiler options
#**************************
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp") # make sure we are using the c preprocessor
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8 -r8 -warn all -nofor-main -gen-interfaces -assume no2underscores")
  set (CMAKE_AR "xiar")   # needed for libraries with -ipo
  if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
     set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ipo -xHost -assume buffered_io -opt-prefetch -auto -unroll4")
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ipo -xHost -opt-prefetch -unroll4")
   else()
     # intel will disable certain optimizations with checks on, so for Release
     # plus debug info we don't add the "-check all" but we do include it for
     # the full debug mode
     if (${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -traceback -fpe:0  ")
     else()
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -check all -traceback -fpe:0  ")
     endif()
   endif()
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8 -fdefault-real-8 -fdefault-double-8 -ffree-line-length-none -Wall -pedantic")
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fcheck=all")
  elseif(${CMAKE_BUILD_TYPE} STREQUAL "Release" OR ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror") # warnings are errors
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-maybe-uninitialized") # this is a rather pointless warning from gfortran wrt to allocatable arrays
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-function") # there are a few unused functions (mainly debug stuff which should be #ifdef'd out)
    if (NOT OMP)
      set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wno-unused-dummy-argument") # in a few spots, locks are handed over as dummies
    endif()
  endif()
endif()

if (${CMAKE_CXX_COMPILER_ID}  STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID}  STREQUAL "Clang")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
endif()

if (OMP)
  # should be using FindOpenMP
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp" )
  else() #applies to both g++ and clang++
    set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp" )
  endif()
endif()

# Make sure we user O0 for debug
if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

#**************************
# Preprocessor definitions
#**************************

add_definitions("-D_REENTRANT -D FORBUF_FAST") # I/O buffer
add_definitions("-DTIGER_USE_BLAS2") # probably shouldn't even be an option anymore
if (OMP)
 add_definitions("-DTIGER_USE_OMP")
endif()


#*****************
# Math libraries
#*****************

# BLAS/LAPACK
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mkl=sequential")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=sequential")
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DTIGER_LAPACK_INT4")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTIGER_LAPACK_INT4")
endif()


#*******************
# Other libraries
#*******************

# libint: can set LIBINT_ROOT_DIR as a helper
find_package(libint REQUIRED)
include_directories(${LIBINT_INCLUDE_DIRS})

# GSL: can set GSL_ROOT_DIR as a helper
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})

# Armadillo
find_package(Armadillo REQUIRED)
include_directories(BEFORE ${ARMADILLO_INCLUDE_DIRS})


IF(${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")
    # OpenMP, currently implies GNU compiler collection on FBSD (looking forward to clang+OMP and "flang"
    include_directories(/usr/local/lib/gcc5/gcc/x86_64-portbld-freebsd11.0/5.3.0/include)
    # FreeBSD specific code: compiling with gfortran5 requires knowledge of its rpath
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=/usr/local/lib/gcc5")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "FreeBSD")


#********************
# Building TigerCI
#********************

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/erkale)
include_directories(${CMAKE_SOURCE_DIR}/IO_buffer)
include_directories(${CMAKE_BINARY_DIR}/erkale/)  # needed for our custom, at compile time, basis library location

# Define all the sub-parts and TigerCI
add_subdirectory(erkale)
add_subdirectory(IO_buffer)
add_subdirectory(tiger_fortran)
add_subdirectory(tiger_cpp)
add_executable(tiger_ci.exe ${ERKALE_SOURCE} ${TIGER_FORTRAN_SOURCE} ${TIGER_CPP_SOURCE} ${IO_BUFFER_SOURCE})

# Link TigerCI to the external librarires erkale wants
target_link_libraries(tiger_ci.exe ${ARMADILLO_LIBRARY} ${LIBINT_LIBRARIES} ${GSL_LIBRARIES})

# Earlier CMAKE version don't have direct handling for C++11
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#set_property(TARGET tiger_ci.exe PROPERTY CXX_STANDARD_REQUIRED CXX_STANDARD 11)

# if we are using the GNU compiler we need to explicit ask to link to LAPACK/BLAS
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  target_link_libraries(tiger_ci.exe ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# note that our lives are a lot easier if intel links
# with ifort and not icpc and then manually add the 
# C++ standard lib
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  SET_TARGET_PROPERTIES(tiger_ci.exe PROPERTIES LINKER_LANGUAGE Fortran)
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -cxxlib")
endif()


