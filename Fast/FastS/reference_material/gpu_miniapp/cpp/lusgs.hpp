// LU-SGS 5 equations 3dfull par hyperplans i+j+k = cte.
// Transcription de FastS/Compute/LU/{invlu_l,invlu_u,lu_i_3dfull,
// lu_dinv,mjr_newton}.for ; meme recurrence, ordre de parcours par plans.
#pragma once
#include "mini_zone.hpp"

namespace mini {

// [ +/- A(l1) + rho(l1) I ]/2 . dW(l1) pour une direction (lu_i_3dfull.for)
// idir : colonne de coe du rayon spectral (1:i, 2:j, 3:k, 0-based)
template <class V2>
KOKKOS_INLINE_FUNCTION void lu_mat(const Params& p, const V2& rop, const V2& coe, const V2& tn,
                                   const V2& dout, const int l1, const int lt1, const int lt2,
                                   const int idir, const double signe, double& b1, double& b2,
                                   double& b3, double& b4, double& b5) {
  const double gam1 = p.gamma;
  const double gamm1 = gam1 - 1.;
  const double cp = p.gamma * p.cv;
  const double gam2 = gamm1 - 1.;

  const double tcx = 0.5 * (tn(lt2, 0) + tn(lt1, 0));
  const double tcy = 0.5 * (tn(lt2, 1) + tn(lt1, 1));
  const double tcz = 0.5 * (tn(lt2, 2) + tn(lt1, 2));

  const double u = rop(l1, 1);
  const double v = rop(l1, 2);
  const double w = rop(l1, 3);
  const double t = rop(l1, 4);
  const double q2 = .5 * (u * u + v * v + w * w);
  const double h = cp * t + q2;
  const double ph2 = gamm1 * q2;
  const double qn = tcx * u + tcy * v + tcz * w;
  const double rad = coe(l1, idir) * signe;

  const double b11 = rad;
  const double b12 = tcx;
  const double b13 = tcy;
  const double b14 = tcz;
  const double b21 = (tcx * ph2 - u * qn);
  const double b22 = (qn - tcx * gam2 * u) + rad;
  const double b23 = (tcy * u - tcx * gamm1 * v);
  const double b24 = (tcz * u - tcx * gamm1 * w);
  const double b25 = tcx * gamm1;
  const double b31 = (tcy * ph2 - v * qn);
  const double b32 = (tcx * v - tcy * gamm1 * u);
  const double b33 = (qn - tcy * gam2 * v) + rad;
  const double b34 = (tcz * v - tcy * gamm1 * w);
  const double b35 = tcy * gamm1;
  const double b41 = (tcz * ph2 - w * qn);
  const double b42 = (tcx * w - tcz * gamm1 * u);
  const double b43 = (tcy * w - tcz * gamm1 * v);
  const double b44 = (qn - tcz * gam2 * w) + rad;
  const double b45 = tcz * gamm1;
  const double b51 = qn * (ph2 - h);
  const double b52 = (tcx * h - gamm1 * u * qn);
  const double b53 = (tcy * h - gamm1 * v * qn);
  const double b54 = (tcz * h - gamm1 * w * qn);
  const double b55 = (gam1 * qn) + rad;

  const double d1 = dout(l1, 0), d2 = dout(l1, 1), d3 = dout(l1, 2), d4 = dout(l1, 3),
               d5 = dout(l1, 4);
  b1 = b11 * d1 + b12 * d2 + b13 * d3 + b14 * d4;
  b2 = b21 * d1 + b22 * d2 + b23 * d3 + b24 * d4 + b25 * d5;
  b3 = b31 * d1 + b32 * d2 + b33 * d3 + b34 * d4 + b35 * d5;
  b4 = b41 * d1 + b42 * d2 + b43 * d3 + b44 * d4 + b45 * d5;
  b5 = b51 * d1 + b52 * d2 + b53 * d3 + b54 * d4 + b55 * d5;
}

KOKKOS_INLINE_FUNCTION void decode_ijk(const Params& p, const int l, int& i, int& j, int& k) {
  const int kk = l / (p.nx * p.ny);
  const int rem = l - kk * p.nx * p.ny;
  const int jj = rem / p.nx;
  const int ii = rem - jj * p.nx;
  i = ii - 1; j = jj - 1; k = kk - 1;
}

// balayage inferieur : (D+L) x = b, plans croissants
inline void lusgs_lower(const Zone& z, const View2& din, const View2& dout) {
  const Params p = z.p;
  const View2 rop = z.rop, coe = z.coe, ti = z.ti, tj = z.tj, tk = z.tk;
  const IView cells = z.plane_cell;
  const double signe = 0.5;
  for (int q = 0; q < z.nplane; ++q) {
    Kokkos::parallel_for(
        "lusgs_lower", Kokkos::RangePolicy<>(z.plane_ptr[q], z.plane_ptr[q + 1]),
        KOKKOS_LAMBDA(const int m) {
          const int l = cells(m);
          int i, j, k;
          decode_ijk(p, l, i, j, k);
          const double xal = coe(l, 0) * signe;
          double x1 = din(l, 0), x2 = din(l, 1), x3 = din(l, 2), x4 = din(l, 3), x5 = din(l, 4);
          double b1, b2, b3, b4, b5;
          if (k > 1) {
            lu_mat(p, rop, coe, tk, dout, l - p.inck, l - p.inck, l, 3, signe, b1, b2, b3, b4, b5);
            x1 += b1 * xal; x2 += b2 * xal; x3 += b3 * xal; x4 += b4 * xal; x5 += b5 * xal;
          }
          if (j > 1) {
            lu_mat(p, rop, coe, tj, dout, l - p.incj, l - p.incj, l, 2, signe, b1, b2, b3, b4, b5);
            x1 += b1 * xal; x2 += b2 * xal; x3 += b3 * xal; x4 += b4 * xal; x5 += b5 * xal;
          }
          if (i > 1) {
            lu_mat(p, rop, coe, ti, dout, l - p.inci, l - p.inci, l, 1, signe, b1, b2, b3, b4, b5);
            x1 += b1 * xal; x2 += b2 * xal; x3 += b3 * xal; x4 += b4 * xal; x5 += b5 * xal;
          }
          const double diag = 1. / coe(l, 4);
          dout(l, 0) = x1 * diag; dout(l, 1) = x2 * diag; dout(l, 2) = x3 * diag;
          dout(l, 3) = x4 * diag; dout(l, 4) = x5 * diag;
        });
  }
}

// balayage superieur : (D+U) x = D y, plans decroissants, puis Newton
inline void lusgs_upper(const Zone& z, const View2& dout, const View2& rop1) {
  const Params p = z.p;
  const View2 rop = z.rop, coe = z.coe, ti = z.ti, tj = z.tj, tk = z.tk;
  const IView cells = z.plane_cell;
  const double signe = -0.5;
  const double cvinv = 1. / p.cv;
  const double cvinv2 = 0.5 * cvinv;
  for (int q = z.nplane - 1; q >= 0; --q) {
    Kokkos::parallel_for(
        "lusgs_upper", Kokkos::RangePolicy<>(z.plane_ptr[q], z.plane_ptr[q + 1]),
        KOKKOS_LAMBDA(const int m) {
          const int l = cells(m);
          int i, j, k;
          decode_ijk(p, l, i, j, k);
          const double xal = coe(l, 0) * signe;
          double diag = coe(l, 4);
          double x1 = dout(l, 0) * diag, x2 = dout(l, 1) * diag, x3 = dout(l, 2) * diag,
                 x4 = dout(l, 3) * diag, x5 = dout(l, 4) * diag;
          double b1, b2, b3, b4, b5;
          if (k < p.nk) {
            lu_mat(p, rop, coe, tk, dout, l + p.inck, l + p.inck, l + 2 * p.inck, 3, signe, b1, b2, b3, b4, b5);
            x1 += b1 * xal; x2 += b2 * xal; x3 += b3 * xal; x4 += b4 * xal; x5 += b5 * xal;
          }
          if (j < p.nj) {
            lu_mat(p, rop, coe, tj, dout, l + p.incj, l + p.incj, l + 2 * p.incj, 2, signe, b1, b2, b3, b4, b5);
            x1 += b1 * xal; x2 += b2 * xal; x3 += b3 * xal; x4 += b4 * xal; x5 += b5 * xal;
          }
          if (i < p.ni) {
            lu_mat(p, rop, coe, ti, dout, l + p.inci, l + p.inci, l + 2 * p.inci, 1, signe, b1, b2, b3, b4, b5);
            x1 += b1 * xal; x2 += b2 * xal; x3 += b3 * xal; x4 += b4 * xal; x5 += b5 * xal;
          }
          diag = 1. / coe(l, 4);
          dout(l, 0) = x1 * diag; dout(l, 1) = x2 * diag; dout(l, 2) = x3 * diag;
          dout(l, 3) = x4 * diag; dout(l, 4) = x5 * diag;

          // mise a jour Newton (mjr_newton.for)
          const double ro_old = rop(l, 0), u_old = rop(l, 1), v_old = rop(l, 2),
                       w_old = rop(l, 3), t_old = rop(l, 4);
          rop1(l, 0) = rop(l, 0) + dout(l, 0);
          const double r_1 = 1. / rop1(l, 0);
          rop1(l, 1) = (ro_old * u_old + dout(l, 1)) * r_1;
          rop1(l, 2) = (ro_old * v_old + dout(l, 2)) * r_1;
          rop1(l, 3) = (ro_old * w_old + dout(l, 3)) * r_1;
          const double roe_old =
              ro_old * (p.cv * t_old + 0.5 * (u_old * u_old + v_old * v_old + w_old * w_old));
          rop1(l, 4) = (roe_old + dout(l, 4)) * r_1 * cvinv -
                       cvinv2 * (rop1(l, 1) * rop1(l, 1) + rop1(l, 2) * rop1(l, 2) +
                                 rop1(l, 3) * rop1(l, 3));
        });
  }
}

inline void lusgs(const Zone& z, const View2& din, const View2& dout, const View2& rop1) {
  lusgs_lower(z, din, dout);
  lusgs_upper(z, dout, rop1);
}

}  // namespace mini
