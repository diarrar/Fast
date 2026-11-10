// Mini-app FastS GPU : description d une zone structuree "compacte"
// (meme disposition memoire que FastS : tableaux (ndimdx, nvar), l contigu)
#pragma once
#include <Kokkos_Core.hpp>
#include <vector>

namespace mini {

struct Params {
  int ni, nj, nk;        // cellules interieures
  int nx, ny, nz;        // dimensions avec 2 cellules fictives par cote
  int ndimdx;            // nx*ny*nz
  int inci, incj, inck;  // increments memoire
  double gamma, cv, psiroe, dt;
};

// indice lineaire 0-based d une cellule (i,j,k) 1-based (interieur 1..ni)
KOKKOS_INLINE_FUNCTION int inddm(const Params& p, int i, int j, int k) {
  return (i + 1) + (j + 1) * p.nx + (k + 1) * p.nx * p.ny;
}

// LayoutLeft : premier indice (l) contigu, comme rop(ndimdx, neq) en Fortran
using View2 = Kokkos::View<double**, Kokkos::LayoutLeft>;
using View1 = Kokkos::View<double*>;
using IView = Kokkos::View<int*>;

struct Zone {
  Params p;
  View2 rop, rop1, drodm, drodm_ref, drodm_out, coe, ti, tj, tk, fli, flj, flk;
  View1 vol;
  // hyperplans i+j+k = cte pour le LU-SGS
  int nplane = 0;
  std::vector<int> plane_ptr;   // host : debut de chaque plan dans plane_cell
  IView plane_cell;             // device : indices l des cellules, plan par plan
};

inline Zone make_zone(int ni, int nj, int nk) {
  Zone z;
  Params& p = z.p;
  p.ni = ni; p.nj = nj; p.nk = nk;
  p.nx = ni + 4; p.ny = nj + 4; p.nz = nk + 4;
  p.ndimdx = p.nx * p.ny * p.nz;
  p.inci = 1; p.incj = p.nx; p.inck = p.nx * p.ny;
  p.gamma = 1.4; p.cv = 2.5; p.psiroe = 0.1; p.dt = 2.0 / ni;
  const int n = p.ndimdx;
  z.rop       = View2("rop", n, 5);
  z.rop1      = View2("rop1", n, 5);
  z.drodm     = View2("drodm", n, 5);
  z.drodm_ref = View2("drodm_ref", n, 5);
  z.drodm_out = View2("drodm_out", n, 5);
  z.coe       = View2("coe", n, 5);
  z.ti        = View2("ti", n, 3);
  z.tj        = View2("tj", n, 3);
  z.tk        = View2("tk", n, 3);
  z.fli       = View2("fli", n, 5);
  z.flj       = View2("flj", n, 5);
  z.flk       = View2("flk", n, 5);
  z.vol       = View1("vol", n);
  return z;
}

}  // namespace mini
