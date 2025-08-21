# add library
add_library(glass_utils STATIC 
            glass_lisa.c glass_lisa.h 
            glass_data.h glass_data.c 
            glass_math.c glass_math.h 
            glass_gmm.c glass_gmm.h 
            glass_wavelet.c glass_wavelet.h 
            glass_galaxy.c glass_galaxy.h
            glass_model.c glass_model.h
            glass_sampler.c glass_sampler.h
            glass_catalog.c glass_catalog.h
            glass_constants.h glass_utils.c glass_utils.h)

include_directories ("${PROJECT_SOURCE_DIR}/external/")
include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS} ${BLAS_INCLUDE_DIRS})

install(DIRECTORY "./" DESTINATION include/glass FILES_MATCHING PATTERN "*.h")

