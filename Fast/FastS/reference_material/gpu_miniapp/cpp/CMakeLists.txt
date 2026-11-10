# Mini-app FastS GPU, phase 0 : version C++ / Kokkos
#   cmake -S . -B build -DKokkos_ROOT=<prefix Kokkos> && cmake --build build
# Le backend (OpenMP, CUDA, HIP, SYCL) est celui de l installation Kokkos.
cmake_minimum_required(VERSION 3.16)
project(fasts_gpu_miniapp CXX)
find_package(Kokkos REQUIRED)
add_executable(miniapp_kokkos main.cpp)
target_link_libraries(miniapp_kokkos Kokkos::kokkos)
target_compile_features(miniapp_kokkos PRIVATE cxx_std_17)
