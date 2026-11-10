#!/bin/sh
# Compile et execute les deux versions de la mini-app, puis compare les
# sommes de controle. A lancer depuis gpu_miniapp/.
#   KOKKOS_ROOT=<prefix> FORTRAN_COMPILER=gnu|nvhpc|amd|cray GPU_ARCH=... scripts/run_all.sh [ni nj nk nrep]
set -e
HERE=$(cd "$(dirname "$0")/.." && pwd)
NI=${1:-64}; NJ=${2:-64}; NK=${3:-64}; NREP=${4:-10}
FC=${FORTRAN_COMPILER:-gnu}
OUT=${OUT_DIR:-$HERE/results}
mkdir -p "$OUT"

echo "== Fortran ($FC)"
make -C "$HERE/fortran" COMPILER="$FC" GPU_ARCH="${GPU_ARCH:-cc90}" > "$OUT/build_fortran.log" 2>&1
"$HERE/fortran/miniapp_fortran" "$NI" "$NJ" "$NK" "$NREP" | tee "$OUT/fortran.out" | grep -E "INFO|TIME|DIFF"

echo "== Kokkos (KOKKOS_ROOT=${KOKKOS_ROOT:?prefix Kokkos requis})"
cmake -S "$HERE/cpp" -B "$OUT/build_kokkos" -DKokkos_ROOT="$KOKKOS_ROOT" -DCMAKE_BUILD_TYPE=Release > "$OUT/build_kokkos.log" 2>&1
cmake --build "$OUT/build_kokkos" >> "$OUT/build_kokkos.log" 2>&1
"$OUT/build_kokkos/miniapp_kokkos" "$NI" "$NJ" "$NK" "$NREP" | tee "$OUT/kokkos.out" | grep -E "INFO|TIME|DIFF"

echo "== Comparaison Fortran (reference CPU) / Kokkos"
python3 "$HERE/scripts/compare_checks.py" "$OUT/fortran.out" "$OUT/kokkos.out" \
  --map flux_cpu=flux_gather --map lusgs_cpu_dw=lusgs_hp_dw --map lusgs_cpu_rop1=lusgs_hp_rop1 --tol 1e-12
