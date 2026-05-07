# Fast RAG

Couche RAG locale pour la racine du dépôt Fast, pensée pour s'appuyer sur un backend LLM déjà déployé.

## Objectifs

- priorité aux sources Fortran issues des listes officielles `srcs.py`
- récupération hybride avec surpondération du noyau solveur
- assemblage d'un contexte sourcé pour API Python, noyau numérique, build/perf/GPU et validation
- aucun besoin Docker/Apptainer

## Modèles par défaut (mai 2026)

- dense embedding général: `Qwen/Qwen3-Embedding-8B`
- dense embedding code: `nomic-ai/nomic-embed-code`
- reranker: `BAAI/bge-reranker-v2-m3`

Ces choix sont configurables dans `tools/fast_rag/fast_rag/config.py`.

## Corpus indexé

- Fortran officiel via:
  - `/home/runner/work/Fast/Fast/Fast/FastS/srcs.py`
  - `/home/runner/work/Fast/Fast/Fast/FastC/srcs.py`
  - `/home/runner/work/Fast/Fast/Fast/Fast/srcs.py`
- Pont Python:
  - `/home/runner/work/Fast/Fast/Fast/FastS/FastS/**/*.py`
  - `/home/runner/work/Fast/Fast/Fast/FastC/FastC/**/*.py`
  - `/home/runner/work/Fast/Fast/Fast/Fast/Fast/**/*.py`
- Documentation:
  - `/home/runner/work/Fast/Fast/docs/**/*.html`

## Exemples

```bash
cd <repo_root>
python -m tools.fast_rag.fast_rag.cli manifest
python -m tools.fast_rag.fast_rag.cli search "farfield boundary condition" --top-k 5
python -m tools.fast_rag.fast_rag.cli prompt "où est calculée la métrique ?" --top-k 4
python -m tools.fast_rag.fast_rag.cli benchmark /path/to/benchmark.jsonl --top-k 5
```

## Format benchmark

Chaque ligne JSONL doit contenir:

```json
{"question":"où est calculée la métrique ?","expected_paths":["Fast/FastC/FastC/Metric/skmtr.for"],"intent":"noyau_numerique"}
```

## Limites actuelles

- le score dense est un hook optionnel: l'intégration effective à Qdrant/embeddings locaux reste branchable sans imposer de dépendance Python ici
- le reranking par défaut est heuristique, avec un emplacement dédié pour brancher un cross-encoder local
- l'assemblage final du prompt est optimisé pour être envoyé à un backend externe déjà disponible
