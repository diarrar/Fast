#!/bin/sh
# Recupere et installe Kokkos pour un backend donne.
#   build_kokkos.sh <prefix_install> openmp
#   build_kokkos.sh <prefix_install> cuda  <arch>   (ex: AMPERE80, HOPPER90)
#   build_kokkos.sh <prefix_install> hip   <arch>   (ex: AMD_GFX90A, AMD_GFX942)
# Variables : KOKKOS_SRC (defaut: ./kokkos-src), KOKKOS_VERSION (defaut 4.5.01)
set -e
PREFIX=${1:?prefix d installation}
BACKEND=${2:-openmp}
ARCH=${3:-}
SRC=${KOKKOS_SRC:-$(pwd)/kokkos-src}
VERSION=${KOKKOS_VERSION:-4.5.01}

if [ ! -d "$SRC" ]; then
  git clone --depth 1 --branch "$VERSION" https://github.com/kokkos/kokkos "$SRC"
fi
BUILD="$SRC/build-$BACKEND"
mkdir -p "$BUILD"
OPTS="-DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_CXX_STANDARD=17 -DKokkos_ENABLE_SERIAL=ON"
case "$BACKEND" in
  openmp) OPTS="$OPTS -DKokkos_ENABLE_OPENMP=ON" ;;
  cuda)   OPTS="$OPTS -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_CUDA=ON -DKokkos_ENABLE_CUDA_LAMBDA=ON -DKokkos_ARCH_${ARCH:?arch GPU requise}=ON -DCMAKE_CXX_COMPILER=$SRC/bin/nvcc_wrapper" ;;
  hip)    OPTS="$OPTS -DKokkos_ENABLE_OPENMP=ON -DKokkos_ENABLE_HIP=ON -DKokkos_ARCH_${ARCH:?arch GPU requise}=ON -DCMAKE_CXX_COMPILER=hipcc" ;;
  *) echo "backend inconnu: $BACKEND"; exit 1 ;;
esac
cmake -S "$SRC" -B "$BUILD" $OPTS
cmake --build "$BUILD" -j 8 --target install
echo "Kokkos ($BACKEND) installe dans $PREFIX"
