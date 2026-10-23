#!/bin/bash
set -o errexit
set -o nounset
set -o pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 INSTALL_PREFIX"
  exit 1
fi

INSTALL_PREFIX="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"
export CFLAGS="-I/opt/local/include/libomp"
export CXXFLAGS="-I/opt/local/include/libomp"
export LDFLAGS="-L/opt/local/lib/libomp"
export OPENBLAS_NUM_THREADS=1

rm -rf build
mkdir -p build
pushd build/
cmake .. \
	-DCMAKE_C_COMPILER=clang \
  	-DCMAKE_CXX_COMPILER=clang++ \
  	-DOpenMP_C_FLAGS="-Xpreprocessor -fopenmp" \
  	-DOpenMP_C_LIB_NAMES="omp" \
  	-DOpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp" \
  	-DOpenMP_CXX_LIB_NAMES="omp" \
  	-DOpenMP_omp_LIBRARY=/opt/local/lib/libomp/libomp.dylib \
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=true

cmake --build . -- VERBOSE=1
cmake --build . --target install
popd

echo ""
echo "*****************************************************************************"
echo "  DONE: GLASS built and installed to: "
echo "      ${INSTALL_PREFIX}"
echo "  To use: "
echo "      export PATH=${INSTALL_PREFIX}/bin:\${PATH}"
echo "*****************************************************************************"

