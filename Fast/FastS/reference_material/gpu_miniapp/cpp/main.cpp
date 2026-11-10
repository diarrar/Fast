// Mini-application FastS GPU, phase 0 (version C++ / Kokkos)
// usage : miniapp_kokkos [ni nj nk nrep]
// Mesure et verifie :
//   flux_gather : flux Roe o3, une cellule = 6 flux (pas de course)
//   flux_faces  : flux Roe o3 par face stocke + bilan
//   lusgs_hp    : LU-SGS par hyperplans + mise a jour Newton
#include <Kokkos_Core.hpp>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "init.hpp"
#include "lusgs.hpp"
#include "mini_zone.hpp"
#include "roe_flux.hpp"

using namespace mini;

static void checksum(const Zone& z, const char* name, const View2& f) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), f);
  const Params& p = z.p;
  for (int v = 0; v < 5; ++v) {
    double s = 0., sa = 0., smax = 0.;
    for (int k = 1; k <= p.nk; ++k)
      for (int j = 1; j <= p.nj; ++j)
        for (int i = 1; i <= p.ni; ++i) {
          const double x = h(inddm(p, i, j, k), v);
          s += x;
          sa += std::fabs(x);
          smax = std::fmax(smax, std::fabs(x));
        }
    std::printf("CHECK %s v%d sum= %24.16E abssum= %24.16E max= %24.16E\n", name, v + 1, s, sa, smax);
  }
}

static void maxdiff(const Zone& z, const char* name, const View2& f1, const View2& f2) {
  auto h1 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), f1);
  auto h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), f2);
  const Params& p = z.p;
  double dmax = 0., fmax = 0.;
  for (int v = 0; v < 5; ++v)
    for (int k = 1; k <= p.nk; ++k)
      for (int j = 1; j <= p.nj; ++j)
        for (int i = 1; i <= p.ni; ++i) {
          const int l = inddm(p, i, j, k);
          dmax = std::fmax(dmax, std::fabs(h1(l, v) - h2(l, v)));
          fmax = std::fmax(fmax, std::fabs(h1(l, v)));
        }
  std::printf("DIFF %s maxabs= %12.4E rel= %12.4E\n", name, dmax, dmax / std::fmax(fmax, 1e-300));
}

template <class F>
static double time_kernel(F&& f, int nrep) {
  f();  // chauffe (compilation JIT, premier lancement)
  Kokkos::fence();
  Kokkos::Timer t;
  for (int it = 0; it < nrep; ++it) f();
  Kokkos::fence();
  return t.seconds() / nrep;
}

int main(int argc, char** argv) {
  Kokkos::ScopeGuard guard(argc, argv);
  int ni = 64, nj = 64, nk = 64, nrep = 10;
  if (argc >= 4) { ni = std::atoi(argv[1]); nj = std::atoi(argv[2]); nk = std::atoi(argv[3]); }
  if (argc >= 5) nrep = std::atoi(argv[4]);

  Zone z = make_zone(ni, nj, nk);
  init_fields(z);
  build_planes(z);
  const double ncell = double(ni) * double(nj) * double(nk);
  std::printf("INFO impl=kokkos backend=%s ni=%d nj=%d nk=%d nrep=%d\n",
              Kokkos::DefaultExecutionSpace::name(), ni, nj, nk, nrep);

  // flux gather
  double t = time_kernel([&] { flux_gather(z); }, nrep);
  std::printf("TIME flux_gather ns_per_cell= %12.4E\n", t / ncell * 1e9);
  checksum(z, "flux_gather", z.drodm);
  Kokkos::deep_copy(z.drodm_ref, z.drodm);

  // flux par faces stockees
  t = time_kernel([&] { flux_faces(z); }, nrep);
  std::printf("TIME flux_faces ns_per_cell= %12.4E\n", t / ncell * 1e9);
  checksum(z, "flux_faces", z.drodm);
  maxdiff(z, "flux_faces_vs_gather", z.drodm_ref, z.drodm);

  // flux par faces, boucle i interne vectorisee (reference CPU)
  t = time_kernel([&] { flux_faces_lines(z); }, nrep);
  std::printf("TIME flux_faces_lines ns_per_cell= %12.4E\n", t / ncell * 1e9);
  checksum(z, "flux_faces_lines", z.drodm);
  maxdiff(z, "flux_faces_lines_vs_gather", z.drodm_ref, z.drodm);

  // LU-SGS par hyperplans (entree : resultat gather)
  t = time_kernel([&] { lusgs(z, z.drodm_ref, z.drodm_out, z.rop1); }, nrep);
  std::printf("TIME lusgs_hp ns_per_cell= %12.4E\n", t / ncell * 1e9);
  checksum(z, "lusgs_hp_dw", z.drodm_out);
  checksum(z, "lusgs_hp_rop1", z.rop1);
  return 0;
}
