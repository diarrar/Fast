# Portage GPU de Fast / FastS (NVIDIA + AMD) : analyse et plan

Date : 2026-11-10. Périmètre : modules `FastC` et `FastS` du dépôt `Fast`
(le cœur du solveur structuré est dans `FastS`, la couche HPC, les transferts
et la gestion mémoire compacte sont dans `FastC`).

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
  gabarits (`template_FluxAndBalance.for`, `fluFaceEuler.for`, `fluvisq_*`).
  Porter les gabarits + une cinquantaine de noyaux manuels suffit pour couvrir
  l'essentiel.
* La première tentative (OpenMP offload) a laissé des briques réutilisables
  (gabarit de flux en formulation "gather", macros de boucles `loopGpu_*`,
  prototypes de CL `bvbs_*_GPU`) mais **n'a jamais été raccordée** : aucune
  directive `target` dans le code, macro `E_GPU` jamais définie dans le build,
  gabarit GPU désynchronisé du gabarit CPU.
* Recommandation : **OpenMP target offload en Fortran** (portable NVIDIA via
  `nvfortran`, AMD via `amdflang`/Cray `ftn`), avec les directives isolées
  derrière des macros de la couche HPC pour pouvoir basculer vers OpenACC ou
  un backend C++ si un compilateur fait défaut. Avancer par paliers :
  explicite monozone → multizone/transferts/RANS → implicite (LU-SGS par
  hyperplans ou GMRES) → optimisation et CI bi-constructeur.
* Point dur n° 1 : la **phase implicite LU-SGS** (récurrence séquentielle en
  i, j, k). Point dur n° 2 : la **maturité des compilateurs Fortran fixe +
  préprocesseur C + offload** sur AMD, à qualifier en phase 0 avec une
  mini-application avant tout engagement.

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
  `param_solver.h`) consommés par tous les noyaux.
* Placement NUMA par premier toucher (`initNuma`) : équivalent GPU = `omp
  target enter data` au `warmup`.

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
(`loop*_begin.for` / `loop_end.for`) : c'est le point d'injection naturel
des directives GPU. Les ~144 fichiers à boucles "en dur" sont majoritairement
dans LU, DTLOC, BC et POST.

Tout le Fortran est en **forme fixe `.for`** (2 264 fichiers, aucun `.f90`),
préprocessé par le préprocesseur C (`#include`, `#ifdef`). Cette combinaison
doit être validée explicitement avec chaque compilateur offload (voir § 5).

---

## 3. Bilan de la première tentative OpenMP offload

Éléments présents dans le dépôt (aucun commit dédié ; arrivés avec les
commits de nettoyage de 2024, e866fa4 / bf8d96c) :

| Brique | État | Verdict |
|---|---|---|
| `Compute/template_FluxAndBalance_gpu.for` | formulation **gather** (`assemble_drodm_mono.for`, flux recalculé des deux côtés avec `sens=±1`), boucles `DO i` + `!$OMP SIMD` | Idée correcte ; gabarit désynchronisé (pas de wiggle, `SENSORTYPE` vs `WIG_AMP`, signature différente). À réécrire à partir du gabarit CPU courant |
| `Compute/loopGpu_begin.for`, `loopGpu_end.for`, `loop3dcartGpu_begin.for` (+ copies dans FastC/HPC_LAYER) | boucles `k,j,i` collapsables, branchées sous `#ifndef E_SCALAR_COMPUTER` | Réutilisable ; mais `E_SCALAR_COMPUTER` est toujours défini dans `setup.scons`, et il faut ajouter la directive `target teams distribute parallel do collapse(3)` |
| `BC/bvbs_wall_viscous_GPU_gpu.for`, `bvbs_wall_isothermal_GPU_gpu.for` | prototypes de CL vectorisables | Non référencés dans `srcs.py` ; à reprendre comme modèle |
| `FastC/Init/copy2gpu.cpp` | fonction Python `copy2gpu`, `#pragma omp target map` commenté, macro `E_GPU` | Squelette seulement ; à remplacer par une vraie gestion de résidence |
| `generate_flu.py ... GPU` | option de génération avec le gabarit GPU | À conserver comme mécanisme |

Pourquoi elle n'a pas abouti (lecture du code) : le build ne définit jamais
`E_GPU`, aucune directive `target` n'existe, la gestion des données
(`enter data`, `update`) n'a pas été écrite, les verrous/synchro n'ont pas
été contournés, et le gabarit GPU a divergé du gabarit CPU. La tentative
s'est arrêtée au stade "boucles réécrites" sans jamais exécuter sur GPU.

---

## 4. Choix technologique

### 4.1 Options comparées

| Option | NVIDIA | AMD | Effort | Commentaire |
|---|---|---|---|---|
| **OpenMP target (Fortran)** | `nvfortran -mp=gpu` : mûr | `amdflang -fopenmp --offload-arch=gfx9xx` (ROCm ≥ 6.2, flang LLVM) ; Cray `ftn` sur HPE | moyen | Garde le Fortran et les gabarits ; cohérent avec la tentative précédente et le style `!$OMP SIMD` déjà en place |
| OpenACC (Fortran) | `nvfortran -acc` : très mûr | pas de support industriel AMD (GCC partiel, Clacc expérimental) | moyen | Non portable AMD → seulement comme backend alternatif NVIDIA |
| `do concurrent` + stdpar | `nvfortran -stdpar` | en cours chez AMD/flang | faible | Pas de contrôle fin de la mémoire ni des réductions ; trop immature AMD |
| Réécriture C++ (Kokkos / SYCL / HIP+CUDA) | oui | oui | élevé | Portabilité maximale, mais abandon du Fortran : à envisager seulement si la génération de code produit des noyaux C++ (voie de repli à long terme) |

### 4.2 Recommandation

OpenMP target offload en Fortran, **avec une couche d'abstraction** :

* toutes les directives GPU passent par des `#include` de `FastC/HPC_LAYER`
  (`GPU_LOOP_BEGIN.for`, `GPU_LOOP_END.for`, `GPU_REDUCTION.for`,
  `GPU_DATA_*.for`) ;
* une seule macro de build `E_GPU` (plus `E_GPU_OMP` / `E_GPU_ACC` pour le
  backend), et `E_SCALAR_COMPUTER` devient le mode CPU explicite ;
* le mode CPU reste strictement inchangé quand `E_GPU` n'est pas défini
  (non-régression garantie, comparaison bit à bit possible en CPU).

Compilateurs cibles de référence : NVIDIA HPC SDK (`nvfortran`, `nvc++`) et
ROCm (`amdflang`, `amdclang++`), avec Cray CCE comme second AMD si la
machine le fournit (Adastra MI250X/MI300A). `gfortran` offload sert de
compilateur de contrôle fonctionnel uniquement (performances faibles).

---

## 5. Points durs et décisions de conception

### 5.1 Décomposition et lancement des noyaux

* Sur GPU : **1 rang MPI = 1 GPU**, et dans un rang, **une tâche = une zone**
  entière (sous-zone unique, `ijkv_sdm = (1,1,1)`, `Nbre_thread_actif_loc =
  1`). La boucle `navier_stokes_struct` devient une suite de noyaux
  (viscosité → dt → src → flux_i → flux_j → flux_k → bilan → mise à jour),
  chacun lancé sur toute la zone avec `collapse(3)`.
* Les appels `SYNCHRO_WAIT` / `SYNCHRO_GO` / `verrou_*` sont compilés en
  no-op sous `E_GPU` (les dépendances sont assurées par l'ordre des noyaux).
* Pour recouvrir les petites zones et le multi-zone : lancements asynchrones
  (`nowait` + `depend`) depuis une région host, ou un thread host par flux
  GPU. Prévu en phase 4.

### 5.2 Bilan de flux sans course d'écriture

Trois formulations possibles ; **choix à trancher par mesure en phase 1**
sur un noyau Roe o3 3dfull :

1. **Gather** (gabarit GPU existant) : chaque cellule calcule ses 6 flux,
   2× le coût flottant, zéro mémoire supplémentaire, pas de course. Souvent
   gagnant sur GPU car le noyau est limité par la bande passante.
2. **Flux de face stockés** : noyau par direction écrivant `flu_i(ndimdx,
   neq)` puis noyau de bilan ; 1× le coût, +3·neq·ndimdx réels (≈ +60 % de
   mémoire en SA), validation plus simple face au CPU.
3. Atomics : à exclure (FP64 atomique lent, surtout AMD).

### 5.3 Phase implicite (LU-SGS)

* Les balayages `invlu_l` / `invlu_u` sont récursifs en i, j, k. Reformuler
  en **balayage par hyperplans** `i + j + k = p` : boucle séquentielle sur
  `p`, parallèle à l'intérieur du plan. Pour un bloc 100³ chaque plan a
  ~10⁴ cellules : occupation suffisante. La `lu_dinv` (inverse de la
  diagonale) est déjà locale.
* Alternatives déjà présentes dans le code et exploitables : relaxation SSOR
  (`NB_RELAX ≥ 2`, `invlussor_*`) qui se prête à une variante Jacobi par
  couleurs, et le solveur Krylov (`IMPLICITSOLVER = 1`, `matvec`, GMRES) dont
  les opérations sont naturellement parallèles.
* Effet sur la convergence : l'ordre des hyperplans ne change pas le
  résultat mathématique du LU-SGS classique (même récurrence) ; la variante
  colorée change le préconditionneur et sera validée sur les cas Newton
  (`PplaneNasaSAMesh136x96_t2..t7`).
* Le mode `LU_MATCH` (recouvrement entre sous-domaines de threads) n'a plus
  de sens avec une tâche par zone : désactivé sous `E_GPU`.

### 5.4 Résidence des données

* `warmup` : après `_compact`, `omp target enter data map(to:)` des tableaux
  compactés de chaque zone (variables, métriques, `coe`, `drodm`, `xmut`,
  `wig`, coordonnées, `cellN`, `param_int`, `param_real`, tableaux de
  raccord `param_int_tc` / `param_real_tc`).
* La rotation `switchPointers__` échange des vues numpy sur la **même**
  mémoire host déjà mappée : aucune copie. Il faut que la table
  host → device reste valide (les trois niveaux de temps mappés).
* Nouvelles fonctions Python : `FastS._updateHost(t, vars)` et
  `FastS._updateDevice(t, vars)` (`omp target update`), appelées par les
  sorties, `display_temporal_criteria`, `_computeStats`, la sauvegarde, et
  les transferts Python. Sans appel explicite, tout reste sur le GPU.
* Nettoyage à la destruction de l'arbre (`exit data`) ; détection d'accès
  host non synchronisé en mode debug (`__DEBUG__`).

### 5.5 Conditions limites et transferts

* CL : noyaux sur fenêtres (`ind_CL`) de petite taille, à lancer tels quels
  sur GPU (21 types `bvbs_*`, les variantes `_d` adjointes exclues). Les CL
  "fichier" (`inflow_fich`, Lund) lisent des données host : à mapper au
  warmup. `correct_coins` idem.
* Raccords / chimère / IBC : `setInterpTransfersFast` est un gather indirect
  (listes de points donneurs/receveurs + coefficients) : parallèle par point
  receveur, portable en `target teams distribute`. Les lois de paroi IBC
  (`setIBCTransfersCommonVar2`, Musker, linelets) demandent une attention
  particulière (boucles de Newton locales).
* MPI : 1 rang par GPU. Étape 1 : tampons de communication recopiés sur le
  host (`update`) autour des échanges existants ; étape 2 : MPI GPU-aware
  (UCX / ROCm-aware) en passant directement les pointeurs device.

### 5.6 Réductions et diagnostics

CFL min/max/moy (`cptst3`), résidus Newton (`lssiter_verif`,
`residus_navier.h`), normes Krylov, débits/efforts (`compute_effort`,
`cp_debitIBM`), senseur de wiggle : `reduction(min/max/+:)` sur les
noyaux concernés. Les chronomètres par thread (`timer_omp`) et
l'équilibrage dynamique n'ont pas de sens sur GPU : neutralisés.

### 5.7 Chaîne de compilation

* `setup.scons` (FastC, FastS) : variable d'environnement `FAST_GPU`
  (`nvidia` / `amd` / `off`) + `FAST_GPU_ARCH` (`cc90`, `gfx942`, …) →
  `-DE_GPU`, options `-mp=gpu -gpu=ccXX` ou `-fopenmp --offload-arch=gfxXXX`,
  bibliothèques d'exécution (`libnvomp`, `libflang_rt`, `libomptarget`).
* Dépendance Cassiopee : `KCore.Dist` fournit compilateurs et drapeaux
  (`getDistUtilsCompilers`, `createFortranBuilder`). Il faudra soit un
  correctif amont dans Cassiopee (ajout d'une option GPU), soit surcharger
  localement `f77compiler` et les drapeaux dans `setup.scons`.
* Mélange d'exécutables OpenMP : l'extension Python et KCore sont compilés
  avec gcc/libgomp ; les noyaux Fortran offload avec le compilateur
  constructeur. Deux runtimes OpenMP dans le même processus est un risque
  réel (déjà rencontré avec Intel). Règle : compiler **tout** FastC/FastS
  (C++ inclus) avec le compilateur constructeur quand `FAST_GPU` est actif.
* Images de référence : Docker `cassiopee486/cassiopee:main` (base actuelle)
  + couche `nvhpc` ou `rocm` ; servent à la CI.

### 5.8 Risques

| Risque | Impact | Atténuation |
|---|---|---|
| `amdflang` sur Fortran fixe + préprocesseur + offload (bugs, perf) | bloquant AMD | mini-app en phase 0 ; repli Cray `ftn` ; repli lointain : génération C++ |
| Pression de registres des noyaux générés (Roe o3 + visqueux + SA inlinés) | perf | découper flux Euler / visqueux en noyaux séparés ; mesurer occupation ; `-gpu=maxregcount` |
| Non-reproductibilité bit à bit CPU/GPU (ordre des sommes) | validation | tolérances dans les tests ; comparaison sur résidus et grandeurs intégrées |
| Deux runtimes OpenMP dans un processus Python | plantages | toolchain unique sous `FAST_GPU` |
| Gestion mémoire des zones (`switchPointers`, réallocation Python) | corruption silencieuse | assertions en mode debug ; API `_updateHost/Device` obligatoire ; tests dédiés |
| Petites zones (multiblocs fins, chimère) | faible occupation | flux asynchrones multi-zones, fusion de noyaux (phase 4) |

---

## 6. Plan de travail

Estimations en semaines-personne pour un développeur connaissant FastS,
hors validation applicative ONERA.

### Phase 0 — Qualification et infrastructure (3 à 4 semaines)

1. **Mini-application** : programme Fortran autonome dans
   `FastS/reference_material/gpu_miniapp/` appelant un noyau généré
   (`fluroe_lamin_o3_3dfull`) et `invlu_l/u` sur des données synthétiques,
   compilable sans Cassiopee. Sert de banc de test compilateurs.
2. Matrice compilateurs : `nvfortran` (HPC SDK courant), `amdflang` (ROCm
   courant), Cray `ftn` si disponible, `gfortran` offload en contrôle. Points
   à vérifier : forme fixe, `#include` imbriqués, `collapse(3)`, réductions,
   `declare target` pour les petites fonctions (`shape_tab_mtr`, `inddm`).
3. Build : options `FAST_GPU`/`FAST_GPU_ARCH` dans les deux `setup.scons`,
   macro `E_GPU`, images Docker NVIDIA et AMD, job CI qui compile.
4. Baseline de performance CPU (1 nœud) sur 3 cas : `computePT.py`
   (Euler cartésien), `PplaneNasaSAMesh136x96_t1` (laminaire explicite),
   `_t2`/`_t4` (implicite, SA) ; métrique : µs par cellule et par sous-pas.
5. Décision formelle sur la formulation du bilan de flux (gather vs flux de
   face) sur la mini-app.

Livrable : rapport de faisabilité par compilateur, build GPU fonctionnel
(sans noyau porté), baseline.

### Phase 1 — Explicite monozone sur GPU (5 à 6 semaines)

1. Couche HPC : nouvelles macros `GPU_LOOP_*`, `SYNCHRO_*` no-op,
   `LOOP_CACHE_*` à 1 bloc, `WORK_DISTRIBUTION` à 1 tâche par zone sous
   `E_GPU`.
2. Gabarits : réécriture de `template_FluxAndBalance_gpu.for` à partir du
   gabarit CPU courant (options wiggle/senseur incluses), option `GPU` de
   `generate_flu.py` branchée dans `srcs.py` ; régénération ROE/AUSM/SENSOR.
3. Noyaux manuels : `invist`, `cptst3` (réductions), `src_term`,
   `core3ark3`, `correct_flux`, `extract_res`, `mjr_ale_3dhomocart`.
4. CL : farfield, wall inviscid/viscous (adia, isotherme), periodique,
   extrapolate, inflow/outflow/outpres ; `correct_coins`.
5. Résidence des données au `warmup`, `_updateHost/_updateDevice`,
   `exit data`.
6. Validation : `computePT`, `PplaneNasaSAMesh136x96_t1`, cas Euler
   cartésien 3D ; comparaison aux références `test.testT` avec tolérance.
7. Mesure : accélération vs 1 nœud CPU sur NVIDIA et AMD ; profil
   (`nsys`/`rocprof`) des 5 noyaux dominants.

Livrable : explicite RK3 Euler/laminaire 2D/3D sur une zone, sur les deux
constructeurs.

### Phase 2 — Multizone, RANS, transferts, MPI (5 à 6 semaines)

1. `setInterpTransfersFast` sur GPU (raccords coïncidents, nearmatch,
   chimère stationnaire), tampons MPI par recopie host puis GPU-aware.
2. SA : `vispalart`, termes sources SA, flux RANS (`fluFaceSA`, `fluvisq`
   6 équations), `assemble_drodm_plus_SA`.
3. LES : `lesvist`, `computePT_mut`.
4. ALE : `mjr_ale`, vitesses d'entraînement, `core3as2_def`.
5. CL restantes : wallmodel, inflow fichier/Lund/newton, injection/débit,
   wallexchange.
6. IBC : transferts IBC et lois de paroi (`setIBCTransfersCommonVar2`).
7. Validation : `cylindreChimerePT_t1..t4`, `FastIBMnoOverlap_t2..t4`, cas
   multiblocs MPI (2 et 4 GPU).

Livrable : explicite multizone RANS/LES, multi-GPU.

### Phase 3 — Implicite (6 à 8 semaines)

1. LU-SGS par hyperplans : `invlu_l`, `invlu_u`, variantes SA et ALE
   (3dfull, 3dhomo, 3dcart, 2d) ; `lu_dinv`, `extrap_coe`.
2. Résidu Newton (`core3as2*`, `residus_navier.h`), `lssiter_verif`,
   affichage des critères temporels depuis le device.
3. Relaxation SSOR (`NB_RELAX ≥ 2`) et solveur Krylov (`matvec`, GMRES,
   préconditionneur) : le chemin Krylov est le repli si les hyperplans
   déçoivent.
4. SFD, extraction de résidus, `mjro_newton`.
5. Validation : `PplaneNasaSAMesh136x96_t2..t7` (implicite SA),
   `computeSteadyNaca*`.

Livrable : implicite Newton/LU-SGS sur GPU.

### Phase 4 — Optimisation, robustesse, CI (4 à 6 semaines, continu)

1. Fusion de noyaux (viscosité + dt + src), découpage Euler/visqueux si
   la pression de registres l'impose, réglage `collapse`/taille de teams
   par constructeur.
2. Recouvrement multi-zones (flux asynchrones), recouvrement MPI/calcul.
3. Explicite local (`DTLOC`), statistiques (`STAT`), post-traitements
   fréquents (`POST`) sur GPU si le profil le justifie.
4. CI nocturne : compilation + tests `FastS/test` sur un nœud NVIDIA et un
   nœud AMD, suivi des performances.
5. Documentation utilisateur (`docs/doc`) : installation GPU, API
   `_updateHost/_updateDevice`, limitations.

### Hors périmètre initial

Couche Python (`layer="python"`), chimère instationnaire `ucData`,
`FastLBM` (en quarantaine), adjoint Tapenade (`ADJOINT`, `*_d.for`),
`FastStrand`, wire model, raccord conservatif nearmatch, outils Pointwise.

---

## 7. Ordre de bataille des fichiers (phase 1)

| Fichier | Action |
|---|---|
| `FastC/setup.scons`, `FastS/setup.scons` | options `FAST_GPU`, `-DE_GPU`, drapeaux et libs constructeur |
| `FastC/FastC/HPC_LAYER/loopGpu_*.for`, nouveaux `GPU_*.for` | macros de boucles, réductions et no-op synchro |
| `FastC/FastC/HPC_LAYER/WORK_DISTRIBUTION_BEGIN.for`, `LOOP_CACHE_*.for`, `SYNCHRO_*.for` | branchement `#ifdef E_GPU` |
| `FastC/FastC/PyTree.py` (`_compact`, `warmup`), `FastC/FastC/Init/copy2gpu.cpp` → `gpu_data.cpp` | `enter data`, `update`, `exit data`, API Python |
| `FastS/FastS/Compute/template_FluxAndBalance_gpu.for`, `generate_flu.py`, `srcs.py` | gabarit gather/face-flux, génération GPU |
| `FastS/FastS/Compute/{invist,cptst3,src_term,core3ark3,correct_flux}.for` | boucles `GPU_LOOP_*`, réductions |
| `FastS/FastS/Compute/{gsdr3,rhs,lhs,bcs}.cpp`, `navier_stokes_struct.for` | une tâche par zone, suppression verrous et timers sous `E_GPU` |
| `FastS/FastS/BC/bvbs_*.for`, `BCzone.cpp` | boucles fenêtre sur GPU (modèle `bvbs_*_GPU_gpu.for`) |
| `FastS/test/*` | tolérances, variantes GPU des cas de non-régression |

---

## 8. Décisions à prendre avant la phase 1

1. Machines cibles prioritaires (NVIDIA A100/H100, AMD MI250X/MI300) et
   accès à au moins un nœud de chaque pour la CI.
2. Niveau de tolérance accepté entre résultats CPU et GPU dans les tests.
3. Accord avec Cassiopee pour l'ajout d'une option GPU dans `KCore.Dist`
   (sinon surcharge locale dans `setup.scons`).
4. Formulation du bilan de flux (gather vs flux de face), après mesure en
   phase 0.
