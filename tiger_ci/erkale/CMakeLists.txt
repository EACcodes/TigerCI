set(ERKALE_SOURCE
erkale/basis.cpp
erkale/boys.cpp
erkale/chebyshev.cpp 
erkale/elements.cpp 
erkale/erifit.cpp
erkale/eritable.cpp
erkale/eriworker.cpp  
erkale/eri_digest.cpp
erkale/lebedev.cpp  
erkale/lobatto.cpp  
erkale/obara-saika.cpp
erkale/solidharmonics.cpp  
erkale/basislibrary.cpp  
erkale/density_fitting.cpp
erkale/eriscreen.cpp  
erkale/eriworker_gentransform.cpp  
erkale/integrals.cpp  
erkale/linalg.cpp
erkale/mathf.cpp 
erkale/settings.cpp 
erkale/stringutil.cpp 
erkale/xyzutils.cpp
erkale/zmatrix.cpp
PARENT_SCOPE)

# custom basis set library location at compile time
set(ERKALE_BASIS_LOCATION "${CMAKE_SOURCE_DIR}/erkale/basis/")
configure_file(custom_basis_library.h.in ${CMAKE_BINARY_DIR}/erkale/custom_basis_library.h)
