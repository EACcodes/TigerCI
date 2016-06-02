# Keep a list of variable names that we need to pass on to
# find_package_handle_standard_args().
set (_check_list)

find_library (LIBDERIV_LIBRARY NAMES deriv HINTS ${LIBINT_ROOT_DIR} PATH_SUFFIXES lib)
mark_as_advanced (LIBDERIV_LIBRARY)
list (APPEND LIBINT_LIBRARIES ${LIBDERIV_LIBRARY})
list (APPEND _check_list LIBDERIV_LIBRARY)

find_library (LIBINT_LIBRARY NAMES int HINTS ${LIBINT_ROOT_DIR} PATH_SUFFIXES lib)
mark_as_advanced (LIBINT_LIBRARY)
list (APPEND LIBINT_LIBRARIES ${LIBINT_LIBRARY})
list (APPEND _check_list LIBINT_LIBRARY)

# make sure to properly look for libint/libint.h as it is used in the code
find_path (LIBINT_INCLUDE_DIRS libint/libint.h 
  HINTS ${LIBINT_ROOT_DIR} PATH_SUFFIXES include)
mark_as_advanced (LIBINT_INCLUDE_DIRS)
list(APPEND _check_list LIBINT_INCLUDE_DIRS)

# Handle the QUIETLY and REQUIRED arguments and set LIBINT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (LIBINT DEFAULT_MSG ${_check_list})
