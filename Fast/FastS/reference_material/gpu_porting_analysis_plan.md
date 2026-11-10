# Portage GPU de Fast / FastS (NVIDIA + AMD) : analyse et plan

Date : 2026-11-10 (révision 2 : option C++ / Kokkos, résultats de la
mini-application de phase 0). Périmètre : modules `FastC` et `FastS` du
dépôt `Fast` (le cœur du solveur structuré est dans `FastS`, la couche
HPC, les transferts et la gestion mémoire compacte sont dans `FastC`).

---

## 1. Résumé exécutif

* Le solveur est un empilement **Python (pyTree) → C++ (orchestration) →
  Fortran 77 fixe préprocessé (noyaux)**. Les tableaux de calcul sont déjà
  **contigus et "compactés"** par zone (`FastC._compact`), ce qui est un atout
  majeur pour une résidence des données sur GPU.
* Le modèle parallèle CPU actuel (décomposition statique par thread OpenMP,
  blocage de cache, **verrous spin-lock** entre threads voisins, flux
  "scatter" sur deux cellules) est **incompatible tel quel** avec un GPU. Il
  faut une formulation "une cellule = un fil GPU" sans verrou ni écriture
  concurrente.
* La volumétrie brute est trompeuse : ~430 000 lignes de Fortran dans FastS,
  mais ~75 % sont **générées** par `generate_flu.py` à partir de quelques
  gabarits. Porter les gabarits + une cinquantaine de noyaux manuels suffit
  pour couvrir l'essentiel.
* La première tentative (OpenMP offload) a laissé des briques réutilisables
  mais **n'a jamais été raccordée** (aucune directive `target`, macro `E_GPU`
  jamais définie, gabarit GPU désynchronisé).
* **Recommandation révisée** : ne pas réécrire Fast, mais réécrire le
  **chemin chaud GPU en C++ généré, porté par Kokkos**, dans une
  bibliothèque séparée `libfastgpu` à API C appelée par `gsdr3.cpp`. Le
  Fortran CPU, l'orchestration C++ et Python restent en place. Le Fortran
  + OpenMP target est conservé comme **plan B** si la qualification de
  phase 0 sur machines réelles tournait à l'avantage du Fortran.
* La **mini-application de phase 0** (`reference_material/gpu_miniapp/`)
  existe dans les deux voies, reproduit l'arithmétique exacte du flux Roe
  o3 et du LU-SGS de FastS, et valide les formulations GPU (gather,
  hyperplans) à 1e-14 près contre les références CPU. Elle attend les
  mesures sur un nœud NVIDIA et un nœud AMD.
* Point dur n° 1 : la **phase implicite LU-SGS** (récurrence séquentielle
  en i, j, k), résolue par balayage en hyperplans, validé dans la mini-app.
  Point dur n° 2 : les **noyaux C++ sur CPU** ne remplacent pas le Fortran
  vectorisé sans travail dédié ; on ne remplace donc pas le chemin CPU.

---

## 2. Architecture actuelle (ce qu'il faut porter)

### 2.1 Chaîne d'appel d'une itération

```
FastS.PyTree._compute(t, metrics, it)                       Python
  └─ fasts._computePT(zones, metrics, ...)                  C++  Compute/computePT.cpp
       └─ boucle nit_c × boucle nstep (sous-itérations RK / Newton)
            └─ gsdr3(...)                                   C++  Compute/gsdr3.cpp
                 ├─ #pragma omp parallel  { rhs.cpp → navier_stokes_struct_()  ; lhs.cpp → invlu_() }
                 ├─ transfert_multiblock.cpp → FastC::setInterpTransfersFast()  (raccords, chimère, IBC)
                 └─ #pragma omp parallel  { bcs.cpp → BCzone() → bvbs_*() }     (conditions limites)
```

* En `layer="c"` (défaut) toute la boucle de sous-itérations tourne en C++
  sans retour à Python : c'est le mode à cibler pour le GPU.
* En `layer="python"` (explicite local, chimère instationnaire `ucData`), les
  transferts passent par `Converter.Mpi` côté Python à chaque sous-pas : hors
  périmètre GPU dans un premier temps.

### 2.2 Noyau RHS : `navier_stokes_struct.for`

Pour chaque tâche (zone, sous-zone, thread) et pour chaque bloc de cache :

1. viscosité laminaire / turbulente (`invist`, `vispalart`, `lesvist`) ;
2. pas de temps / CFL (`cptst3`, réductions min/max) ;
3. init du RHS + termes sources (`src_term`, SA, pulse, IBM) ;
4. flux convectifs + visqueux + bilan (`fluroe_select`, `fluausm_select`,
   `flusenseur_select` → noyaux générés) puis `correct_flux` ;
5. mise à jour : `core3ark3` (RK3 explicite), `core3as2*` (résidu Newton
   implicite), `core3_dtloc`, `core3as2_kry` (Krylov).

Entre ces étapes : `SYNCHRO_WAIT` / `SYNCHRO_GO` (verrous `lok`, routine
`verrou.for` = boucle d'attente active sur un entier avec `!$OMP FLUSH`).

### 2.3 Modèle de parallélisme CPU (à remplacer)

| Mécanisme | Fichiers | Problème sur GPU |
|---|---|---|
| Décomposition statique zone → sous-zones → threads (`ipt_omp`, `distributeThreads.cpp`, `souszones_list`) | FastC/Compute, FastS/Compute/rhs.cpp | Le GPU veut un seul lancement de noyau par zone, pas N sous-domaines de thread |
| Blocage de cache `ijkv_sdm`, `CACHEBLCK*` (`LOOP_CACHE_BEGIN.for`) | FastC/HPC_LAYER | Inutile ; à neutraliser (1 bloc) |
| Spin-locks inter-threads (`verrou.for`, `synchro_omp_scater.for`, `verrou_lhs.h`) | FastC/HPC_LAYER, FastS/Compute | Impossible sur GPU ; remplacés par la séparation en noyaux successifs |
| Bilan de flux "scatter" : `drodm(l) += flu ; drodm(l0) -= flu` (`assemble_drodm_plus_vec1.for`, `flu_send_scater_vec1.for`) | FastS/Compute | **Course d'écriture** entre fils GPU voisins |
| Boucles linéarisées sur `l` avec `!$OMP SIMD` (`loopI_begin.for`) | FastS/Compute | Il faut des boucles `i,j,k` collapsées (déjà prototypé : `loopGpu_begin.for`) |
| Balayages LU-SGS récursifs (`invlu_l.for`, `invlu_u.for`, variantes SA/ALE/overlap) | FastS/Compute/LU | Récurrence en i (et j, k) : parallélisme nul sans reformulation |

### 2.4 Données

* Variables primitives `rop` (ndimdx × neq, avec `SHIFTVAR` de padding),
  `rop_m1`, `rop_p1` (rotation de pointeurs côté host par
  `switchPointers__`), `drodm`, `coe`, `xmut`, `wig`, métriques `ti/tj/tk/vol`
  (`GET_TI`), vitesses d'entraînement `vent*`, coordonnées, `cellN`.
* Chaque zone porte `param_int` / `param_real` (267 clés dans
  `param_solver.h`) consommés par tous les noyaux, déclarés partout en
  taille supposée `(0:*)`.
* Placement NUMA par premier toucher (`initNuma`) : équivalent GPU = copie
  sur le device au `warmup`.

### 2.5 Volumétrie et taux de génération

| Répertoire | Fichiers .for | Lignes | Sous-programmes | Nature |
|---|---|---|---|---|
| Compute/ROE | 432 | 146 571 | 253 | générés (`generate_flu.py`) |
| Compute/AUSM | 409 | 79 245 | 231 | générés |
| Compute/SENSOR | 420 | 75 548 | 230 | générés |
| BC | 265 | 45 541 | 53 | partiellement générés (`generate_BC_d.py`) |
| Compute/LU | 117 | 23 004 | 40 | manuels, très répétitifs |
| Compute/SA | 97 | 15 709 | 21 | manuels |
| Compute/DTLOC | 48 | 7 264 | 47 | manuels (explicite local) |
| POST, STAT, DISPLAY, Metric, Init | ~115 | ~9 100 | — | post-traitement (host) |

758 sites de boucles passent déjà par des `#include` de la couche HPC
(`loop*_begin.for` / `loop_end.for`). Tout le Fortran est en **forme fixe
`.for`** (2 264 fichiers, aucun `.f90`), préprocessé par le préprocesseur
C (`#include`, `#ifdef`).

---

## 3. Bilan de la première tentative OpenMP offload

Éléments présents dans le dépôt (aucun commit dédié ; arrivés avec les
commits de nettoyage de 2024, e866fa4 / bf8d96c) :

| Brique | État | Verdict |
|---|---|---|
| `Compute/template_FluxAndBalance_gpu.for` | formulation **gather** (`assemble_drodm_mono.for`, flux recalculé des deux côtés avec `sens=±1`) | Idée correcte, reprise dans la mini-app ; gabarit désynchronisé du gabarit CPU |
| `Compute/loopGpu_begin.for`, `loopGpu_end.for`, `loop3dcartGpu_begin.for` | boucles `k,j,i` collapsables sous `#ifndef E_SCALAR_COMPUTER` | Réutilisable pour le plan B ; `E_SCALAR_COMPUTER` est toujours défini dans `setup.scons` |
| `BC/bvbs_wall_viscous_GPU_gpu.for`, `bvbs_wall_isothermal_GPU_gpu.for` | prototypes de CL vectorisables | Non référencés dans `srcs.py` ; modèle pour les CL |
| `FastC/Init/copy2gpu.cpp` | `#pragma omp target map` commenté, macro `E_GPU` | Squelette seulement |
| `generate_flu.py ... GPU` | option de génération avec le gabarit GPU | Mécanisme à conserver (il émettra du C++) |

Pourquoi elle n'a pas abouti : le build ne définit jamais `E_GPU`, aucune
directive `target` n'existe, la gestion des données (`enter data`,
`update`) n'a pas été écrite, les verrous/synchro n'ont pas été contournés,
et le gabarit GPU a divergé du gabarit CPU.

---

## 4. Choix technologique

### 4.1 Pourquoi pas une réécriture complète de Fast

L'orchestration (`computePT.cpp`, `gsdr3.cpp`, transferts) est déjà en
C++. Le Fortran CPU est tuné (blocage de cache, verrous fins) et validé :
le réécrire n'apporte rien sur CPU et coûte un à deux ans de validation.
Le post-traitement, les statistiques, l'explicite local et l'adjoint n'ont
rien à faire sur GPU au départ. Le périmètre de la réécriture est donc le
**chemin chaud GPU seulement**.

### 4.2 Pourquoi des noyaux C++ plutôt que du Fortran offloadé

Quel que soit le langage, les noyaux GPU doivent être restructurés
(gather ou flux de face, hyperplans LU, suppression des verrous) : il y
aura deux chemins de calcul. La différence est leur avenir :

* Fortran + OpenMP target : deux variantes Fortran qui divergent, et une
  dépendance durable aux compilateurs Fortran offload, dont `amdflang` sur
  du F77 forme fixe préprocessé est le maillon faible.
* C++ généré + Kokkos : compilateurs GPU C++ de premier rang chez les deux
  constructeurs, langue de Cassiopee, outillage (profileurs, sanitizers),
  possibilité à long terme d'un backend CPU unique. La mini-app montre
  toutefois que ce dernier point n'est **pas gratuit** (voir § 9) : on ne
  l'inscrit pas au périmètre.

Le vrai actif du code est `generate_flu.py` et ses gabarits. Étendre le
générateur pour émettre du C++ coûte quelques milliers de lignes ; les
parties manuelles du chemin chaud sont très répétitives (40 routines LU
d'un même motif, 21 CL quasi identiques) et se factorisent par templates.
Estimation du chemin chaud GPU en C++ : 15 000 à 25 000 lignes.

### 4.3 Modèles de programmation C++ comparés

| Modèle | NVIDIA | AMD | Avis |
|---|---|---|---|
| **Kokkos** | oui | oui | Premier choix : vues avec espace mémoire explicite, `MDRangePolicy` pour les stencils 3D, `RangePolicy` par hyperplan pour le LU, backend OpenMP pour tester sur CPU, écosystème large |
| CUDA/HIP natif derrière une fine couche de macros | oui | oui | Contrôle et performance maximaux, dépendance nulle, mais deux backends à entretenir et pas de CPU |
| stdpar (`nvc++ -stdpar`, `hipstdpar`) | oui | oui | Trop implicite : mémoire managée imposée, pas de contrôle du lancement ; à réserver au post-traitement |
| SYCL | via plugins | via plugins | Fonctionnel mais peu répandu dans les centres français ; Kokkos a un backend SYCL si besoin |
| Fortran + OpenMP target (plan B) | `nvfortran` mûr | `amdflang` / Cray `ftn` à qualifier | Garde le Fortran ; dépend des compilateurs Fortran offload |

### 4.4 Décision et critère de bascule

Voie principale : **C++ généré + Kokkos**. Plan B : Fortran + OpenMP
target. La mini-app de phase 0 fournit les deux implémentations ; la
décision finale est prise sur les mesures NVIDIA et AMD selon ce critère :
si le Fortran offloadé compile et tourne sur AMD avec moins de 30 % d'écart
de performance par rapport à Kokkos et sans contournement de compilateur,
le plan B redevient discutable au vu de l'effort moindre ; sinon la voie
principale est confirmée.

---

## 5. Points durs et décisions de conception

### 5.1 Architecture `libfastgpu`

* Bibliothèque C++ compilée avec le compilateur constructeur (nvcc /
  hipcc via Kokkos), exposant une **API C** : `fastgpu_init`,
  `fastgpu_map_zone(param_int, param_real, pointeurs host)`,
  `fastgpu_rhs(zone, nstep)`, `fastgpu_lhs(zone)`, `fastgpu_bc(zone)`,
  `fastgpu_update_host(zone, champ)`, `fastgpu_update_device(...)`,
  `fastgpu_finalize`.
* `gsdr3.cpp` l'appelle sous `E_GPU` ; le module Python et KCore restent
  compilés avec gcc. Cela isole les chaînes de compilation et évite deux
  runtimes OpenMP dans un même processus.
* Les vues Kokkos sont **non gérées** (`Kokkos::View<double**,
  LayoutLeft, Unmanaged>`) côté host sur les tableaux numpy compactés, et
  gérées côté device ; la disposition `(ndimdx, neq)` avec `l` contigu est
  conservée telle quelle (validé dans la mini-app).

### 5.2 Décomposition et lancement des noyaux

* **1 rang MPI = 1 GPU**, et dans un rang, **une tâche = une zone**
  entière. La boucle `navier_stokes_struct` devient une suite de noyaux
  (viscosité → dt → src → flux → bilan → mise à jour), chacun lancé sur
  toute la zone avec `MDRangePolicy<Rank<3>>`.
* Les appels `SYNCHRO_WAIT` / `SYNCHRO_GO` / `verrou_*` n'existent pas
  dans `libfastgpu` ; les dépendances sont assurées par l'ordre des
  noyaux et les `fence`.
* Recouvrement des petites zones et multi-zone : instances d'espace
  d'exécution Kokkos (flux CUDA/HIP) par zone, phase 4.

### 5.3 Bilan de flux sans course d'écriture

Deux formulations implémentées dans la mini-app, à départager sur GPU :

1. **Gather** : chaque cellule calcule ses 6 flux, 2× le coût flottant,
   zéro mémoire supplémentaire. Sur CPU : 2× plus lent que la formulation 2.
2. **Flux de face stockés** : un noyau par direction écrit le flux de
   chaque face, puis un noyau de bilan ; +3·neq·ndimdx réels.

Les atomiques sont exclus (FP64 atomique lent, surtout AMD).

### 5.4 Phase implicite (LU-SGS)

* Balayage par **hyperplans** `i + j + k = p` : listes de cellules par
  plan construites au `warmup`, boucle séquentielle sur `p`, parallèle à
  l'intérieur du plan. Même récurrence que le LU-SGS lexicographique de
  FastS : la mini-app le confirme à 4e-16 près.
* Coût : 3N plans par balayage (190 à 64³), donc latence de lancement à
  surveiller sur petites zones ; regroupement de plans ou exécution
  asynchrone multi-zones en phase 4.
* Alternatives présentes dans le code : relaxation SSOR (`NB_RELAX ≥ 2`),
  solveur Krylov (`IMPLICITSOLVER = 1`, `matvec`, GMRES).
* `LU_MATCH` (recouvrement entre sous-domaines de threads) n'a plus de
  sens avec une tâche par zone.

### 5.5 Résidence des données

* `warmup` : après `_compact`, `fastgpu_map_zone` copie sur le device les
  tableaux compactés de chaque zone (variables, métriques, `coe`,
  `drodm`, `xmut`, `wig`, coordonnées, `cellN`, `param_int`,
  `param_real`, tableaux de raccord).
* La rotation `switchPointers__` échange des vues numpy sur la **même**
  mémoire host : la table host → device reste valide si les trois
  niveaux de temps sont mappés.
* Nouvelles fonctions Python : `FastS._updateHost(t, vars)` et
  `FastS._updateDevice(t, vars)`, appelées par les sorties,
  `display_temporal_criteria`, `_computeStats`, la sauvegarde et les
  transferts Python. Sans appel explicite, tout reste sur le GPU.

### 5.6 Conditions limites et transferts

* CL : noyaux sur fenêtres (`ind_CL`) de petite taille, 21 types
  `bvbs_*`, factorisés par templates C++. Les CL "fichier" (`inflow_fich`,
  Lund) lisent des données host : à mapper au warmup.
* Raccords / chimère / IBC : `setInterpTransfersFast` est un gather
  indirect (listes de points et coefficients), parallèle par point
  receveur. Les lois de paroi IBC (Musker, linelets) ont des boucles de
  Newton locales : noyau dédié.
* MPI : 1 rang par GPU. Étape 1 : tampons recopiés sur le host autour
  des échanges existants ; étape 2 : MPI GPU-aware.

### 5.7 Réductions et diagnostics

CFL min/max/moy (`cptst3`), résidus Newton (`lssiter_verif`), normes
Krylov, débits/efforts : `parallel_reduce` Kokkos. Les chronomètres par
thread (`timer_omp`) et l'équilibrage dynamique sont neutralisés sous
`E_GPU`.

### 5.8 Chaîne de compilation

* `libfastgpu` : CMake + Kokkos (backends OpenMP pour les tests CPU, CUDA,
  HIP), installé comme dépendance externe ; `setup.scons` de FastS lie la
  bibliothèque et définit `E_GPU` quand `FAST_GPU=nvidia|amd`.
* Le reste (FastC, FastS Fortran, extension Python) reste compilé par
  `KCore.Dist` sans changement de chaîne.
* Images de référence : Docker `cassiopee486/cassiopee:main` + couche
  `nvhpc` ou `rocm` ; servent à la CI.

### 5.9 Risques

| Risque | Impact | Atténuation |
|---|---|---|
| Pression de registres des noyaux générés (Roe o3 + visqueux + SA) | perf | découper Euler / visqueux ; mesurer l'occupation ; `launch_bounds` |
| Non-reproductibilité bit à bit CPU/GPU (ordre des sommes) | validation | tolérances dans les tests (1e-12 relatif démontré sur la mini-app) |
| Deux chaînes de compilation (gcc pour Python, nvcc/hipcc pour `libfastgpu`) | intégration | API C, pas d'OpenMP host dans `libfastgpu` |
| Gestion mémoire des zones (`switchPointers`, réallocation Python) | corruption silencieuse | assertions en mode debug ; API `_updateHost/Device` obligatoire |
| Petites zones (multiblocs fins, chimère) | faible occupation | instances d'exécution par zone, fusion de noyaux (phase 4) |
| Générateur Python à étendre pour émettre du C++ | effort initial | commencer par le gabarit Roe/AUSM 3dfull, généraliser ensuite |
| Compétence C++/Kokkos de l'équipe | planning | formation courte ; mini-app comme modèle de code |

---

## 6. Plan de travail

Estimations en semaines-personne pour un développeur connaissant FastS,
hors validation applicative ONERA.

### Phase 0 — Qualification et infrastructure (3 à 4 semaines)

Fait dans ce dépôt (`reference_material/gpu_miniapp/`) :

1. Mini-application Fortran (forme fixe + cpp, comme FastS) : flux Roe o3
   scatter (référence), flux gather OpenMP target, LU-SGS lexicographique
   (référence), LU-SGS par hyperplans OpenMP target.
2. Mini-application C++ / Kokkos : flux gather, flux de face stockés,
   LU-SGS par hyperplans, mêmes données synthétiques.
3. Script de comparaison des sommes de contrôle, scripts de build
   (Kokkos openmp / cuda / hip) et d'exécution.
4. Validation croisée sur CPU : toutes les formulations concordent à
   mieux que 1e-13 (détails § 9).

Reste à faire :

5. Exécuter sur un nœud NVIDIA (nvfortran, Kokkos CUDA) et un nœud AMD
   (amdflang ou Cray ftn, Kokkos HIP) ; remplir le tableau du README.
6. Trancher gather vs flux de face, et confirmer la voie principale selon
   le critère du § 4.4.
7. Build : option `FAST_GPU` dans `setup.scons`, macro `E_GPU`, images
   Docker NVIDIA et AMD, job CI qui compile `libfastgpu`.
8. Baseline de performance CPU (1 nœud) sur `computePT.py`,
   `PplaneNasaSAMesh136x96_t1`, `_t2`/`_t4` ; métrique : µs par cellule et
   par sous-pas.

### Phase 1 — Explicite monozone sur GPU (5 à 6 semaines)

1. Squelette `libfastgpu` : API C, mapping d'une zone, vues non gérées sur
   les tableaux compactés, `update_host/device`, intégration dans
   `gsdr3.cpp` sous `E_GPU`.
2. Générateur : `generate_flu.py` émet des noyaux C++ (gabarit
   `template_FluxAndBalance.cpp.in`) pour ROE/AUSM/SENSOR, 3dfull/3dhomo/
   3dcart/2d, Euler/laminaire ; la mini-app sert de modèle.
3. Noyaux manuels : `invist`, `cptst3` (réductions), `src_term`,
   `core3ark3`, `correct_flux`, `extract_res`, `mjr_ale_3dhomocart`.
4. CL : farfield, wall inviscid/viscous (adia, isotherme), periodique,
   extrapolate, inflow/outflow/outpres ; `correct_coins`.
5. Couche HPC CPU : une tâche par zone et synchro no-op quand `E_GPU`.
6. Validation : `computePT`, `PplaneNasaSAMesh136x96_t1`, cas Euler 3D ;
   comparaison aux références `test.testT` avec tolérance.
7. Mesure : accélération vs 1 nœud CPU sur NVIDIA et AMD ; profil
   (`nsys` / `rocprof`) des noyaux dominants.

### Phase 2 — Multizone, RANS, transferts, MPI (5 à 6 semaines)

1. `setInterpTransfersFast` sur GPU (raccords coïncidents, nearmatch,
   chimère stationnaire), tampons MPI par recopie host puis GPU-aware.
2. SA : `vispalart`, termes sources SA, flux RANS 6 équations.
3. LES : `lesvist`, `computePT_mut`.
4. ALE : `mjr_ale`, vitesses d'entraînement, `core3as2_def`.
5. CL restantes : wallmodel, inflow fichier/Lund/newton, injection/débit,
   wallexchange.
6. IBC : transferts IBC et lois de paroi.
7. Validation : `cylindreChimerePT_t1..t4`, `FastIBMnoOverlap_t2..t4`, cas
   multiblocs MPI (2 et 4 GPU).

### Phase 3 — Implicite (6 à 8 semaines)

1. LU-SGS par hyperplans : variantes 5/6 équations, ALE, 3dfull/3dhomo/
   3dcart/2d, en templates C++ à partir de la mini-app.
2. Résidu Newton (`core3as2*`, `residus_navier.h`), `lssiter_verif`.
3. Relaxation SSOR et solveur Krylov (repli si les hyperplans déçoivent).
4. SFD, extraction de résidus, `mjro_newton`.
5. Validation : `PplaneNasaSAMesh136x96_t2..t7`, `computeSteadyNaca*`.

### Phase 4 — Optimisation, robustesse, CI (4 à 6 semaines, continu)

1. Fusion de noyaux, découpage Euler/visqueux si la pression de registres
   l'impose, réglage des tailles de tuiles par constructeur.
2. Instances d'exécution par zone (recouvrement multi-zones),
   recouvrement MPI/calcul, regroupement de plans LU.
3. Explicite local (`DTLOC`), statistiques, post-traitements fréquents.
4. CI nocturne : compilation + tests `FastS/test` sur un nœud NVIDIA et un
   nœud AMD, suivi des performances.
5. Documentation utilisateur : installation GPU, API
   `_updateHost/_updateDevice`, limitations.

### Hors périmètre initial

Couche Python (`layer="python"`), chimère instationnaire `ucData`,
`FastLBM`, adjoint Tapenade, `FastStrand`, wire model, raccord conservatif
nearmatch, remplacement du chemin CPU Fortran par les noyaux C++.

---

## 7. Ordre de bataille des fichiers (phase 1)

| Fichier | Action |
|---|---|
| `FastS/libfastgpu/` (nouveau) | CMake, API C, `zone.hpp`, noyaux générés et manuels, tests unitaires |
| `FastS/setup.scons` | option `FAST_GPU`, `-DE_GPU`, lien `libfastgpu` |
| `FastS/FastS/Compute/generate_flu.py` + gabarit C++ | émission des noyaux de flux C++ |
| `FastS/FastS/Compute/{gsdr3,rhs,lhs,bcs}.cpp`, `navier_stokes_struct.for` | appels `libfastgpu` sous `E_GPU`, une tâche par zone, verrous et timers neutralisés |
| `FastC/FastC/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for`, `LOOP_CACHE_*.for`, `SYNCHRO_*.for` | branchement `#ifdef E_GPU` |
| `FastC/FastC/PyTree.py` (`_compact`, `warmup`), `FastS/FastS/PyTree.py` | `fastgpu_map_zone`, `_updateHost/_updateDevice` |
| `FastS/test/*` | tolérances, variantes GPU des cas de non-régression |

---

## 8. Décisions à prendre avant la phase 1

1. Machines cibles prioritaires (NVIDIA A100/H100, AMD MI250X/MI300) et
   accès à au moins un nœud de chaque pour la phase 0 puis la CI.
2. Confirmation de la voie C++ / Kokkos après les mesures de phase 0
   (critère § 4.4).
3. Niveau de tolérance accepté entre résultats CPU et GPU dans les tests.
4. Formulation du bilan de flux (gather vs flux de face), après mesure
   sur GPU.

---

## 9. Résultats de la mini-application (phase 0, sans GPU)

Conteneur x86 4 cœurs, gcc / gfortran 13.3, `-O3 -march=native`, zone
64³, 5 répétitions ; temps en ns par cellule et par appel, bruit ~10 %.
Les régions `target` Fortran sont exécutées sur le host : leurs temps ne
préjugent de rien sur GPU.

| Noyau | Fortran 1 thread | Kokkos OpenMP 1 thread | Kokkos OpenMP 4 threads |
|---|---|---|---|
| flux scatter FastS (référence) | 58 | — | — |
| flux gather | 339 (target sur host) | 325 à 340 | 91 à 113 |
| flux faces stockées + bilan | — | 148 à 188 | 50 à 68 |
| LU-SGS lexicographique (référence) | 99 | — | — |
| LU-SGS hyperplans | 410 (target sur host) | 485 à 620 | 126 à 171 |

| Comparaison | écart relatif max |
|---|---|
| gather (Fortran target) vs scatter (Fortran) | 4.7e-14 |
| hyperplans (Fortran target) vs lexicographique (Fortran) | 4.2e-16 |
| faces (Kokkos) vs gather (Kokkos) | 2.3e-15 |
| Fortran vs Kokkos, flux et LU-SGS | 9.2e-15 |

Enseignements :

1. Les tableaux à taille supposée `param_int(0:*)` de FastS ne peuvent
   pas être mappés implicitement dans une région `target` : en plan B,
   chaque noyau devrait recevoir des formes explicites ou des clauses
   `map` dimensionnées (mécanique, mais à faire par le générateur).
2. Les fonctions instruction (`inddm`, `indmtr`) et les macros du
   préprocesseur passent avec gfortran dans une région `target` ; à
   confirmer avec nvfortran et amdflang.
3. Sur CPU le gather coûte 2× les flux de face stockés ; sur GPU le
   rapport peut s'inverser.
4. Les noyaux Kokkos écrits naïvement sont 2,5 à 6× plus lents sur CPU
   que le Fortran FastS vectorisé par lignes : le chemin CPU Fortran est
   conservé.
5. Le LU-SGS par hyperplans produit 380 lancements par appel à 64³ : la
   latence de lancement sera le facteur limitant sur petites zones.
