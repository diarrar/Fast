// Donnees synthetiques, identiques a fortran/init.for (meme ordre des
// operations pour obtenir des entrees bit a bit identiques).
#pragma once
#include <cmath>
#include "mini_zone.hpp"

namespace mini {

inline void init_fields(Zone& z) {
  const Params& p = z.p;
  auto rop = Kokkos::create_mirror_view(z.rop);
  auto coe = Kokkos::create_mirror_view(z.coe);
  auto ti  = Kokkos::create_mirror_view(z.ti);
  auto tj  = Kokkos::create_mirror_view(z.tj);
  auto tk  = Kokkos::create_mirror_view(z.tk);
  auto vol = Kokkos::create_mirror_view(z.vol);

  const double dx = 1.0 / p.ni, dy = 1.0 / p.nj, dz = 1.0 / p.nk;
  const double twopi = 8.0 * std::atan(1.0);
  const double gam = p.gamma;
  const double rgp = p.cv * (gam - 1.0);
  const double dt = p.dt;

  for (int k = -1; k <= p.nk + 2; ++k)
    for (int j = -1; j <= p.nj + 2; ++j)
      for (int i = -1; i <= p.ni + 2; ++i) {
        const int l = inddm(p, i, j, k);
        const double xc = (i - 0.5) * dx, yc = (j - 0.5) * dy, zc = (k - 0.5) * dz;
        const double xf = (i - 1) * dx, yf = (j - 1) * dy, zf = (k - 1) * dz;
        rop(l, 0) = 1.0 + 0.1 * std::sin(twopi * xc) * std::cos(twopi * yc) * std::cos(twopi * zc);
        rop(l, 1) = 0.5 + 0.1 * std::cos(twopi * xc) * std::sin(twopi * yc);
        rop(l, 2) = 0.1 * std::sin(twopi * zc) * std::cos(twopi * xc);
        rop(l, 3) = 0.05 * std::cos(twopi * yc) * std::sin(twopi * zc);
        rop(l, 4) = 1.0 + 0.05 * std::cos(twopi * xc) * std::cos(twopi * yc) * std::sin(twopi * zc);

        ti(l, 0) = dy * dz * (1.0 + 0.1 * std::sin(twopi * yc) * std::cos(twopi * zc));
        ti(l, 1) = dy * dz * 0.05 * std::cos(twopi * xf);
        ti(l, 2) = dy * dz * 0.02 * std::sin(twopi * zc);
        tj(l, 0) = dx * dz * 0.05 * std::sin(twopi * yf);
        tj(l, 1) = dx * dz * (1.0 + 0.1 * std::cos(twopi * xc) * std::sin(twopi * zc));
        tj(l, 2) = dx * dz * 0.02 * std::cos(twopi * xc);
        tk(l, 0) = dx * dy * 0.02 * std::cos(twopi * zf);
        tk(l, 1) = dx * dy * 0.05 * std::sin(twopi * xc);
        tk(l, 2) = dx * dy * (1.0 + 0.1 * std::sin(twopi * xc) * std::cos(twopi * yc));
        vol(l) = dx * dy * dz * (1.0 + 0.1 * std::sin(twopi * xc) * std::sin(twopi * yc) * std::cos(twopi * zc));
        coe(l, 0) = 0.0; coe(l, 1) = 0.0; coe(l, 2) = 0.0; coe(l, 3) = 0.0; coe(l, 4) = 1.0;
      }

  for (int k = 1; k <= p.nk; ++k)
    for (int j = 1; j <= p.nj; ++j)
      for (int i = 1; i <= p.ni; ++i) {
        const int l = inddm(p, i, j, k);
        const double u = rop(l, 1), v = rop(l, 2), w = rop(l, 3), t = rop(l, 4);
        const double c = std::sqrt(gam * rgp * t);
        coe(l, 0) = dt / vol(l);
        double sx, sy, sz, qn;
        sx = 0.5 * (ti(l, 0) + ti(l + p.inci, 0));
        sy = 0.5 * (ti(l, 1) + ti(l + p.inci, 1));
        sz = 0.5 * (ti(l, 2) + ti(l + p.inci, 2));
        qn = u * sx + v * sy + w * sz;
        coe(l, 1) = 2.0 * (std::fabs(qn) + c * std::sqrt(sx * sx + sy * sy + sz * sz));
        sx = 0.5 * (tj(l, 0) + tj(l + p.incj, 0));
        sy = 0.5 * (tj(l, 1) + tj(l + p.incj, 1));
        sz = 0.5 * (tj(l, 2) + tj(l + p.incj, 2));
        qn = u * sx + v * sy + w * sz;
        coe(l, 2) = 2.0 * (std::fabs(qn) + c * std::sqrt(sx * sx + sy * sy + sz * sz));
        sx = 0.5 * (tk(l, 0) + tk(l + p.inck, 0));
        sy = 0.5 * (tk(l, 1) + tk(l + p.inck, 1));
        sz = 0.5 * (tk(l, 2) + tk(l + p.inck, 2));
        qn = u * sx + v * sy + w * sz;
        coe(l, 3) = 2.0 * (std::fabs(qn) + c * std::sqrt(sx * sx + sy * sy + sz * sz));
        coe(l, 4) = 1.0 + coe(l, 0) * (coe(l, 1) + coe(l, 2) + coe(l, 3)) * 0.5;
      }

  Kokkos::deep_copy(z.rop, rop);
  Kokkos::deep_copy(z.coe, coe);
  Kokkos::deep_copy(z.ti, ti);
  Kokkos::deep_copy(z.tj, tj);
  Kokkos::deep_copy(z.tk, tk);
  Kokkos::deep_copy(z.vol, vol);
  Kokkos::deep_copy(z.drodm, 0.0);
  Kokkos::deep_copy(z.drodm_out, 0.0);
  Kokkos::deep_copy(z.rop1, 0.0);
}

// Listes de cellules par hyperplan i+j+k = cte (ordre k,j,i dans le plan)
inline void build_planes(Zone& z) {
  const Params& p = z.p;
  z.nplane = p.ni + p.nj + p.nk - 2;
  std::vector<int> count(z.nplane + 1, 0);
  for (int k = 1; k <= p.nk; ++k)
    for (int j = 1; j <= p.nj; ++j)
      for (int i = 1; i <= p.ni; ++i) count[i + j + k - 2]++;
  z.plane_ptr.assign(z.nplane + 1, 0);
  for (int q = 1; q <= z.nplane; ++q) z.plane_ptr[q] = z.plane_ptr[q - 1] + count[q];
  std::vector<int> fill(z.plane_ptr.begin(), z.plane_ptr.end());
  z.plane_cell = IView("plane_cell", p.ni * p.nj * p.nk);
  auto cells = Kokkos::create_mirror_view(z.plane_cell);
  for (int k = 1; k <= p.nk; ++k)
    for (int j = 1; j <= p.nj; ++j)
      for (int i = 1; i <= p.ni; ++i) {
        const int q = i + j + k - 3;  // plan 0-based
        cells(fill[q]++) = inddm(p, i, j, k);
      }
  Kokkos::deep_copy(z.plane_cell, cells);
}

}  // namespace mini
