# add library
add_library(glass_external STATIC 
            kissfft/kiss_fft.c kissfft/kiss_fft.h
            kissfft/kiss_fftr.c kissfft/kiss_fftr.h
            kissfft/_kiss_fft_guts.h
            kissfft/kiss_fft_log.h
            astrometry/astrometry.c astrometry/astrometry.h) 

# The whole project is built with -fopenmp. kissfft would otherwise fork an
# OpenMP parallel region on every FFT (called millions of times in the F-stat
# proposal and per-step waveforms); for GLASS's small narrow-band transforms
# that fork/join cost dwarfs the FFT and causes nested parallelism inside the
# already-parallel chain/F-stat loops. Keep kissfft strictly serial.
set_source_files_properties(kissfft/kiss_fft.c kissfft/kiss_fftr.c
            PROPERTIES COMPILE_DEFINITIONS GLASS_KISS_FFT_NO_OPENMP)

include_directories ("${PROJECT_SOURCE_DIR}/external/kissfft/" "${PROJECT_SOURCE_DIR}/external/astrometry/")

install(DIRECTORY "astrometry/" DESTINATION include/astrometry FILES_MATCHING PATTERN "*.h")
install(DIRECTORY "kissfft/" DESTINATION include/kissfft FILES_MATCHING PATTERN "kiss_fft.h")
install(DIRECTORY "kissfft/" DESTINATION include/kissfft FILES_MATCHING PATTERN "kiss_fftr.h")
