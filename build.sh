export OMPI_CC=/opt/homebrew/opt/llvm/bin/clang
CFLAGS="-I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib" \
        cmake ../ \
          -DCMAKE_INSTALL_PREFIX="$HOME/.local" \
          -DCMAKE_BUILD_TYPE=Debug \
          -DCMAKE_PREFIX_PATH="/opt/homebrew/opt/openblas;/opt/homebrew/opt/hdf5;/opt/homebrew/opt/open-mpi" \
          -DBLA_VENDOR=OpenBLAS \
          -DOpenMP_C_FLAGS="-fopenmp -L/opt/homebrew/opt/libomp/lib" \
          -DOpenMP_C_LIB_NAMES=omp \
          -DOpenMP_omp_LIBRARY=/opt/homebrew/opt/libomp/lib/libomp.dylib
