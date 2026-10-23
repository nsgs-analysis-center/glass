# add library

set(MBH_SRC
    glass_mbh_io.h glass_mbh_io.c
    glass_mbh_data.h glass_mbh_data.c
    glass_mbh_prior.h glass_mbh_prior.c
    glass_mbh_proposal.h glass_mbh_proposal.c
    glass_mbh_sampler.h glass_mbh_sampler.c
    glass_mbh_model.h glass_mbh_model.c
    glass_mbh_waveform.h glass_mbh_waveform.c
    glass_mbh_IMRPhenom.h
    glass_mbh_IMRPhenom.c
    glass_mbh.h
)

set(LAL_SRC
    IMRPhenomD_internals.c
    IMRPhenomD_internals.h
    IMRPhenomD.c
    IMRPhenomD.h
    IMRPhenomT.c
    IMRPhenomT.h
)

add_library(glass_mbh STATIC
    ${MBH_SRC}
    ${LAL_SRC}
)

include_directories ("${PROJECT_SOURCE_DIR}/external/")
include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories ("${PROJECT_SOURCE_DIR}/mbh/src/")
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
include_directories(SYSTEM ${GSL_INCLUDE_DIRS})

install(DIRECTORY "./" DESTINATION include/glass FILES_MATCHING PATTERN "*.h")

