// Flux de Roe o3 Euler 3dfull : transcription de
// FastS/Compute/ROE/3dfull/fluFaceEuler_o3_3dfull_*.for et includes.
#pragma once
#include "mini_zone.hpp"

namespace mini {

// Flux sur la face (lf - 1/2) dans la direction d increment inc.
template <class V2>
KOKKOS_INLINE_FUNCTION void roe_o3_face(const Params& p, const V2& rop, const V2& tn,
                                        const int lf, const int inc,
                                        double& flu1, double& flu2, double& flu3,
                                        double& flu4, double& flu5) {
  const double gam = p.gamma;
  const double rgp = p.cv * (gam - 1.);
  const double gam1 = gam / (gam - 1.);
  const double psiroe = p.psiroe;
  double c6 = 1. / 6.;
  const double c4 = 5. * c6;
  const double c5 = 2. * c6;
  c6 = -1. * c6;

  const double tcx = tn(lf, 0);
  const double tcy = tn(lf, 1);
  const double tcz = tn(lf, 2);
  const double si = Kokkos::sqrt(tcx * tcx + tcy * tcy + tcz * tcz);

  const int nm = lf - inc;
  const int nm2 = lf - 2 * inc;
  const int np = lf + inc;

  // pentes : qm etat droit, qp etat gauche
  double qm1 = c4 * rop(lf, 0) + c5 * rop(nm, 0) + c6 * rop(np, 0);
  double qp1 = c4 * rop(nm, 0) + c6 * rop(nm2, 0) + c5 * rop(lf, 0);
  const double qm2 = c4 * rop(lf, 1) + c5 * rop(nm, 1) + c6 * rop(np, 1);
  const double qp2 = c4 * rop(nm, 1) + c6 * rop(nm2, 1) + c5 * rop(lf, 1);
  const double qm3 = c4 * rop(lf, 2) + c5 * rop(nm, 2) + c6 * rop(np, 2);
  const double qp3 = c4 * rop(nm, 2) + c6 * rop(nm2, 2) + c5 * rop(lf, 2);
  const double qm4 = c4 * rop(lf, 3) + c5 * rop(nm, 3) + c6 * rop(np, 3);
  const double qp4 = c4 * rop(nm, 3) + c6 * rop(nm2, 3) + c5 * rop(lf, 3);
  const double qm5 = c4 * rop(lf, 4) + c5 * rop(nm, 4) + c6 * rop(np, 4);
  const double qp5 = c4 * rop(nm, 4) + c6 * rop(nm2, 4) + c5 * rop(lf, 4);

  // etats gauche / droit et moyenne de Roe
  const double r1 = qp1;
  const double p1 = r1 * qp5 * rgp;
  const double h1 = gam1 * qp5 * rgp + .5 * (qp2 * qp2 + qp3 * qp3 + qp4 * qp4);
  const double r2 = qm1;
  const double p2 = r2 * qm5 * rgp;
  const double h2 = gam1 * qm5 * rgp + .5 * (qm2 * qm2 + qm3 * qm3 + qm4 * qm4);

  qp1 = Kokkos::sqrt(r1);
  qm1 = Kokkos::sqrt(r2);
  const double r = qp1 * qm1;
  const double r_1 = 1. / (qp1 + qm1);
  const double u = (qp1 * qp2 + qm1 * qm2) * r_1;
  const double v = (qp1 * qp3 + qm1 * qm3) * r_1;
  const double w = (qp1 * qp4 + qm1 * qm4) * r_1;
  const double h = (qp1 * h1 + qm1 * h2) * r_1;
  const double q = .5 * (u * u + v * v + w * w);
  const double c = Kokkos::sqrt((gam - 1.) * Kokkos::fabs(h - q));

  const double qn1 = qp2 * tcx + qp3 * tcy + qp4 * tcz;
  const double qn2 = qm2 * tcx + qm3 * tcy + qm4 * tcz;

  const double du = (qm2 - qp2) * si;
  const double dv = (qm3 - qp3) * si;
  const double dw = (qm4 - qp4) * si;
  const double dp = p2 - p1;
  const double dqn = qn2 - qn1;

  const double s_1 = 1.0 / Kokkos::fmax(si, 1.E-30);
  const double nx = tcx * s_1;
  const double ny = tcy * s_1;
  const double nz = tcz * s_1;

  const double qn = u * nx + v * ny + w * nz;

  double f1 = qn;
  const double f4 = psiroe * (Kokkos::fabs(f1) + c);

  double f2 = 0.5 + Kokkos::copysign(0.5, Kokkos::fabs(f1) - f4);
  double f3 = (f2 - 1.0) * (f1 * f1 / f4 + f4) * 0.5 - f2 * Kokkos::fabs(f1);

  f1 = (r2 - r1 - dp / (c * c)) * si;

  flu1 = f3 * (f1);
  flu2 = f3 * (f1 * u + r * (du - nx * dqn));
  flu3 = f3 * (f1 * v + r * (dv - ny * dqn));
  flu4 = f3 * (f1 * w + r * (dw - nz * dqn));
  flu5 = f3 * (f1 * q + r * (u * du + v * dv + w * dw - qn * dqn));

  f1 = qn + c;
  f2 = 0.5 + Kokkos::copysign(0.5, Kokkos::fabs(f1) - f4);
  f3 = (f2 - 1.) * (f1 * f1 / f4 + f4) * 0.5 - f2 * Kokkos::fabs(f1);

  f1 = (dp + r * c * dqn * s_1) / (2. * c * c) * si * f3;

  flu1 = flu1 + f1;
  flu2 = flu2 + f1 * (u + nx * c);
  flu3 = flu3 + f1 * (v + ny * c);
  flu4 = flu4 + f1 * (w + nz * c);
  flu5 = flu5 + f1 * (h + qn * c);

  f1 = qn - c;
  f2 = 0.5 + Kokkos::copysign(0.5, Kokkos::fabs(f1) - f4);
  f3 = (f2 - 1.) * (f1 * f1 / f4 + f4) * 0.5 - f2 * Kokkos::fabs(f1);

  f1 = (dp - r * c * dqn * s_1) / (2. * c * c) * f3 * si;

  flu1 = flu1 + f1;
  flu2 = flu2 + f1 * (u - nx * c);
  flu3 = flu3 + f1 * (v - ny * c);
  flu4 = flu4 + f1 * (w - nz * c);
  flu5 = flu5 + f1 * (h - qn * c);

  const double rou1 = r1 * qn1;
  const double rou2 = r2 * qn2;
  const double p1p2 = p1 + p2;

  flu1 = (flu1 + rou1 + rou2) * 0.5;
  flu2 = (flu2 + rou1 * qp2 + rou2 * qm2 + tcx * p1p2) * 0.5;
  flu3 = (flu3 + rou1 * qp3 + rou2 * qm3 + tcy * p1p2) * 0.5;
  flu4 = (flu4 + rou1 * qp4 + rou2 * qm4 + tcz * p1p2) * 0.5;
  flu5 = (flu5 + rou1 * h1 + rou2 * h2) * 0.5;
}

using Policy3 = Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;

// Formulation "gather" : chaque cellule calcule ses 6 flux de face.
inline void flux_gather(const Zone& z) {
  const Params p = z.p;
  const View2 rop = z.rop, ti = z.ti, tj = z.tj, tk = z.tk, drodm = z.drodm;
  Kokkos::parallel_for(
      "flux_gather", Policy3({1, 1, 1}, {p.nk + 1, p.nj + 1, p.ni + 1}),
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const int l = inddm(p, i, j, k);
        double f1, f2, f3, f4, f5;
        roe_o3_face(p, rop, ti, l, p.inci, f1, f2, f3, f4, f5);
        double d1 = f1, d2 = f2, d3 = f3, d4 = f4, d5 = f5;
        roe_o3_face(p, rop, ti, l + p.inci, p.inci, f1, f2, f3, f4, f5);
        d1 -= f1; d2 -= f2; d3 -= f3; d4 -= f4; d5 -= f5;
        roe_o3_face(p, rop, tj, l, p.incj, f1, f2, f3, f4, f5);
        d1 += f1; d2 += f2; d3 += f3; d4 += f4; d5 += f5;
        roe_o3_face(p, rop, tj, l + p.incj, p.incj, f1, f2, f3, f4, f5);
        d1 -= f1; d2 -= f2; d3 -= f3; d4 -= f4; d5 -= f5;
        roe_o3_face(p, rop, tk, l, p.inck, f1, f2, f3, f4, f5);
        d1 += f1; d2 += f2; d3 += f3; d4 += f4; d5 += f5;
        roe_o3_face(p, rop, tk, l + p.inck, p.inck, f1, f2, f3, f4, f5);
        drodm(l, 0) = d1 - f1;
        drodm(l, 1) = d2 - f2;
        drodm(l, 2) = d3 - f3;
        drodm(l, 3) = d4 - f4;
        drodm(l, 4) = d5 - f5;
      });
}

// Formulation "flux de face stockes" : un noyau par direction ecrit le flux
// de chaque face, puis un noyau de bilan. 1x le cout flottant, +3 tableaux.
inline void flux_faces(const Zone& z) {
  const Params p = z.p;
  const View2 rop = z.rop, ti = z.ti, tj = z.tj, tk = z.tk, drodm = z.drodm;
  const View2 fli = z.fli, flj = z.flj, flk = z.flk;
  Kokkos::parallel_for(
      "flux_faces_i", Policy3({1, 1, 1}, {p.nk + 1, p.nj + 1, p.ni + 2}),
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const int l = inddm(p, i, j, k);
        double f1, f2, f3, f4, f5;
        roe_o3_face(p, rop, ti, l, p.inci, f1, f2, f3, f4, f5);
        fli(l, 0) = f1; fli(l, 1) = f2; fli(l, 2) = f3; fli(l, 3) = f4; fli(l, 4) = f5;
      });
  Kokkos::parallel_for(
      "flux_faces_j", Policy3({1, 1, 1}, {p.nk + 1, p.nj + 2, p.ni + 1}),
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const int l = inddm(p, i, j, k);
        double f1, f2, f3, f4, f5;
        roe_o3_face(p, rop, tj, l, p.incj, f1, f2, f3, f4, f5);
        flj(l, 0) = f1; flj(l, 1) = f2; flj(l, 2) = f3; flj(l, 3) = f4; flj(l, 4) = f5;
      });
  Kokkos::parallel_for(
      "flux_faces_k", Policy3({1, 1, 1}, {p.nk + 2, p.nj + 1, p.ni + 1}),
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const int l = inddm(p, i, j, k);
        double f1, f2, f3, f4, f5;
        roe_o3_face(p, rop, tk, l, p.inck, f1, f2, f3, f4, f5);
        flk(l, 0) = f1; flk(l, 1) = f2; flk(l, 2) = f3; flk(l, 3) = f4; flk(l, 4) = f5;
      });
  Kokkos::parallel_for(
      "flux_balance", Policy3({1, 1, 1}, {p.nk + 1, p.nj + 1, p.ni + 1}),
      KOKKOS_LAMBDA(const int k, const int j, const int i) {
        const int l = inddm(p, i, j, k);
        for (int v = 0; v < 5; ++v)
          drodm(l, v) = (fli(l, v) - fli(l + p.inci, v)) + (flj(l, v) - flj(l + p.incj, v)) +
                        (flk(l, v) - flk(l + p.inck, v));
      });
}

// Variante CPU : boucles sur les lignes (k,j), boucle i interne vectorisee.
// Sert a mesurer ce que coute la perte de vectorisation du MDRange sur CPU ;
// sur GPU la boucle interne serait sequentielle par fil (a ne pas utiliser).
using Policy2 = Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;

inline void flux_faces_lines(const Zone& z) {
  const Params p = z.p;
  const View2 rop = z.rop, ti = z.ti, tj = z.tj, tk = z.tk, drodm = z.drodm;
  const View2 fli = z.fli, flj = z.flj, flk = z.flk;
  Kokkos::parallel_for(
      "flux_lines_i", Policy2({1, 1}, {p.nk + 1, p.nj + 1}), KOKKOS_LAMBDA(const int k, const int j) {
        const int l0 = inddm(p, 1, j, k);
#pragma omp simd
        for (int i = 0; i <= p.ni; ++i) {
          double f1, f2, f3, f4, f5;
          roe_o3_face(p, rop, ti, l0 + i, p.inci, f1, f2, f3, f4, f5);
          fli(l0 + i, 0) = f1; fli(l0 + i, 1) = f2; fli(l0 + i, 2) = f3; fli(l0 + i, 3) = f4; fli(l0 + i, 4) = f5;
        }
      });
  Kokkos::parallel_for(
      "flux_lines_j", Policy2({1, 1}, {p.nk + 1, p.nj + 2}), KOKKOS_LAMBDA(const int k, const int j) {
        const int l0 = inddm(p, 1, j, k);
#pragma omp simd
        for (int i = 0; i < p.ni; ++i) {
          double f1, f2, f3, f4, f5;
          roe_o3_face(p, rop, tj, l0 + i, p.incj, f1, f2, f3, f4, f5);
          flj(l0 + i, 0) = f1; flj(l0 + i, 1) = f2; flj(l0 + i, 2) = f3; flj(l0 + i, 3) = f4; flj(l0 + i, 4) = f5;
        }
      });
  Kokkos::parallel_for(
      "flux_lines_k", Policy2({1, 1}, {p.nk + 2, p.nj + 1}), KOKKOS_LAMBDA(const int k, const int j) {
        const int l0 = inddm(p, 1, j, k);
#pragma omp simd
        for (int i = 0; i < p.ni; ++i) {
          double f1, f2, f3, f4, f5;
          roe_o3_face(p, rop, tk, l0 + i, p.inck, f1, f2, f3, f4, f5);
          flk(l0 + i, 0) = f1; flk(l0 + i, 1) = f2; flk(l0 + i, 2) = f3; flk(l0 + i, 3) = f4; flk(l0 + i, 4) = f5;
        }
      });
  Kokkos::parallel_for(
      "flux_lines_balance", Policy2({1, 1}, {p.nk + 1, p.nj + 1}), KOKKOS_LAMBDA(const int k, const int j) {
        const int l0 = inddm(p, 1, j, k);
        for (int v = 0; v < 5; ++v)
#pragma omp simd
          for (int i = 0; i < p.ni; ++i) {
            const int l = l0 + i;
            drodm(l, v) = (fli(l, v) - fli(l + p.inci, v)) + (flj(l, v) - flj(l + p.incj, v)) +
                          (flk(l, v) - flk(l + p.inck, v));
          }
      });
}

}  // namespace mini
