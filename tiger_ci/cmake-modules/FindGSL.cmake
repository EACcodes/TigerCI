# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set (_check_list)

find_library (GSL_LIBRARY NAMES gsl HINTS ${GSL_ROOT_DIR} PATH_SUFFIXES lib)
mark_as_advanced (GSL_LIBRARY)
list (APPEND GSL_LIBRARIES ${GSL_LIBRARY})
list (APPEND _check_list GSL_LIBRARY)

find_library (GSLCBLAS_LIBRARY NAMES gslcblas HINTS ${GSL_ROOT_DIR} PATH_SUFFIXES lib)
mark_as_advanced (GSLCBLAS_LIBRARY)
list (APPEND GSL_LIBRARIES ${GSLCBLAS_LIBRARY})
list (APPEND _check_list GSLCBLAS_LIBRARY)

find_path (GSL_INCLUDE_DIR gsl/gsl_rng.h 
  HINTS ${GSL_ROOT_DIR} PATH_SUFFIXES include)
mark_as_advanced (GSL_INCLUDE_DIR)
list(APPEND _check_list GSL_INCLUDE_DIR)

# Handle the QUIETLY and REQUIRED arguments and set GSL_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GSL DEFAULT_MSG ${_check_list})
