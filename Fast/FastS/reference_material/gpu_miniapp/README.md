# Mini-application FastS GPU (phase 0)

Banc d'essai autonome, sans dépendance à Cassiopee, qui reproduit à
l'identique l'arithmétique de deux noyaux représentatifs de FastS :

* **flux de Roe o3 Euler 3dfull** (`fluroe_euler_o3_3dfull.for` et ses
  includes `fluFaceEuler_o3_3dfull_*`, `o3_slope_var`, `etat_roe_GD`,
  `fludiffer_3dfull_*`) ;
* **LU-SGS 5 équations 3dfull** avec mise à jour Newton (`invlu_l.for`,
  `invlu_u.for`, `lu_i_3dfull.for`, `lu_dinv.for`, `mjr_newton.for`).

Il sert à qualifier les chaînes de compilation GPU (NVIDIA, AMD), à
comparer les deux voies candidates (Fortran + OpenMP target, C++ + Kokkos)
et à trancher la formulation du bilan de flux (gather ou flux de face).

Toutes les versions partent des **mêmes données synthétiques** (champ
primitif régulier, métrique 3dfull déformée, coefficients implicites) et
impriment des sommes de contrôle comparables entre langages et machines.

## Arborescence

```
gpu_miniapp/
  fortran/        forme fixe .for + préprocesseur C, comme FastS
    flux_cpu.for      flux Roe, formulation "scatter" FastS (référence CPU)
    flux_gpu.for      flux Roe, formulation "gather", OpenMP target collapse(3)
    lusgs_cpu.for     LU-SGS lexicographique (référence CPU)
    lusgs_gpu.for     LU-SGS par hyperplans i+j+k, OpenMP target par plan
    flux_face_roe_o3.h, lu_mat.h, mjr_newton.h   includes partagés CPU/GPU
    init.for, main.for, Makefile
  cpp/            C++17 + Kokkos
    roe_flux.hpp      flux_gather, flux_faces (flux stockés + bilan), flux_faces_lines
    lusgs.hpp         LU-SGS par hyperplans + Newton
    init.hpp, mini_zone.hpp, main.cpp, CMakeLists.txt
  scripts/
    compare_checks.py comparaison des sommes de contrôle de deux sorties
    build_kokkos.sh   installation de Kokkos (openmp | cuda | hip)
    run_all.sh        compile, exécute et compare les deux versions
```

Conventions communes : zone `ni×nj×nk` avec 2 cellules fictives par côté,
tableaux `(ndimdx, nvar)` avec `l` contigu (disposition compacte FastS),
`param_int` / `param_real` avec les mêmes clés que `param_solver.h`.

## Compilation et exécution

### Fortran

```bash
cd fortran
make COMPILER=gnu                     # gfortran, régions target exécutées sur le host
make COMPILER=nvhpc GPU_ARCH=cc90     # nvfortran, offload NVIDIA (H100 : cc90, A100 : cc80)
make COMPILER=amd   GPU_ARCH=gfx942   # amdflang ROCm, offload AMD (MI300 : gfx942, MI250 : gfx90a)
make COMPILER=cray                    # Cray ftn (charger craype-accel-nvidia90 ou craype-accel-amd-gfx90a)
./miniapp_fortran 64 64 64 10         # ni nj nk nrep
```

Les réels par défaut sont en double précision (`-r8`), comme dans FastS.

### C++ / Kokkos

```bash
scripts/build_kokkos.sh $HOME/kokkos-openmp openmp
scripts/build_kokkos.sh $HOME/kokkos-cuda   cuda HOPPER90
scripts/build_kokkos.sh $HOME/kokkos-hip    hip  AMD_GFX942
cmake -S cpp -B build -DKokkos_ROOT=$HOME/kokkos-cuda -DCMAKE_BUILD_TYPE=Release
cmake --build build
./build/miniapp_kokkos 64 64 64 10
```

### Tout enchaîner et comparer

```bash
KOKKOS_ROOT=$HOME/kokkos-openmp FORTRAN_COMPILER=gnu scripts/run_all.sh 64 64 64 10
```

Le script compare les sommes de contrôle Fortran (référence CPU) et Kokkos
avec une tolérance relative de 1e-12 et renvoie un code d'erreur sinon.

## Sorties

* `TIME <noyau> ns_per_cell=` : temps par cellule et par appel, après un
  appel de chauffe, moyenné sur `nrep` appels (`fence` Kokkos inclus).
* `CHECK <noyau> v<n> sum= abssum= max=` : sommes de contrôle sur les
  cellules intérieures, par variable.
* `DIFF <a>_vs_<b> maxabs= rel=` : écart maximal entre deux formulations
  dans le même exécutable.

## Ce qu'il faut rapporter par machine

| Machine / compilateur | Compile ? | flux CPU réf. | flux gather GPU | flux faces GPU | LU-SGS CPU réf. | LU-SGS hyperplans GPU | écart max |
|---|---|---|---|---|---|---|---|
| ex. H100, nvfortran 25.x | | | | | | | |
| ex. H100, Kokkos CUDA | | | | | | | |
| ex. MI300A, amdflang ROCm 6.x | | | | | | | |
| ex. MI300A, Kokkos HIP | | | | | | | |
| ex. MI250X, Cray ftn | | | | | | | |

Compléter aussi : version du compilateur, options, et pour le Fortran les
diagnostics d'offload (`-Minfo=mp` chez NVIDIA) sur les fonctions
instruction (`inddm`, `indmtr`) et les listes `PRIVATE`.

## Premiers résultats (conteneur x86 4 cœurs, sans GPU)

gcc / gfortran 13.3, `-O3 -march=native`, zone 64³, 5 répétitions.
Les temps sont en ns par cellule et par appel ; le bruit de mesure est de
l'ordre de 10 %. Les régions `target` sont exécutées sur le host (aucun
GPU disponible ici) : leurs temps ne préjugent de rien sur GPU.

| Noyau | Fortran 1 thread | Kokkos OpenMP 1 thread | Kokkos OpenMP 4 threads |
|---|---|---|---|
| flux scatter FastS (référence) | 58 | — | — |
| flux gather | 339 (target sur host) | 325 à 340 | 91 à 113 |
| flux faces stockées + bilan | — | 148 à 188 | 50 à 68 |
| flux faces, boucle i vectorisée | — | 136 | 85 |
| LU-SGS lexicographique (référence) | 99 | — | — |
| LU-SGS hyperplans | 410 (target sur host) | 485 à 620 | 126 à 171 |

Validation :

| Comparaison | écart relatif max |
|---|---|
| gather (Fortran target) vs scatter (Fortran) | 4.7e-14 |
| hyperplans (Fortran target) vs lexicographique (Fortran) | 4.2e-16 |
| faces (Kokkos) vs gather (Kokkos) | 2.3e-15 |
| Fortran vs Kokkos, flux et LU-SGS, 64³ | 9.2e-15 |

## Observations de phase 0 (sans GPU)

1. **Tableaux à taille supposée.** FastS déclare partout `param_int(0:*)`
   et `param_real(0:*)`. OpenMP refuse de les mapper implicitement dans une
   région `target` : chaque noyau GPU devra recevoir des formes explicites
   ou des clauses `map` dimensionnées. C'est mécanique mais touche tous
   les noyaux générés, donc le générateur doit le faire.
2. **Les fonctions instruction** (`inddm`, `indmtr`) et les macros `TN`
   du préprocesseur passent avec gfortran dans une région `target`. À
   confirmer avec nvfortran et amdflang (point de vigilance du plan).
3. **Formulation du bilan.** Sur CPU, le gather coûte 2× les flux de face
   stockés, conformément au double calcul de flux. Sur GPU le rapport peut
   s'inverser (bande passante) : c'est la mesure attendue des machines.
4. **Noyaux Kokkos sur CPU.** Écrits naïvement, ils sont 2,5 à 6× plus
   lents que le Fortran FastS vectorisé par lignes ; forcer `omp simd`
   sur la boucle i n'a pas suffi. Remplacer le chemin CPU Fortran par les
   noyaux C++ n'est donc pas gratuit et ne fait pas partie du périmètre
   initial.
5. **LU-SGS par hyperplans.** À 64³ il y a 190 plans par balayage, soit
   380 lancements de noyau par appel. Sur GPU la latence de lancement
   dominera pour les petites zones : prévoir le regroupement de plans ou
   l'exécution asynchrone multi-zones.

## Suite immédiate

* Exécuter sur un nœud NVIDIA (nvfortran + Kokkos CUDA) et un nœud AMD
  (amdflang ou Cray ftn + Kokkos HIP), remplir le tableau ci-dessus.
* Ajouter, si le temps le permet, un flux visqueux laminaire et la
  variante SA (6 équations) pour mesurer la pression de registres.
