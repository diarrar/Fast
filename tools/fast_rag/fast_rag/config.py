from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


DEFAULT_EMBEDDING_MODEL = 'Qwen/Qwen3-Embedding-8B'
DEFAULT_CODE_EMBEDDING_MODEL = 'nomic-ai/nomic-embed-code'
DEFAULT_RERANKER_MODEL = 'BAAI/bge-reranker-v2-m3'
DEFAULT_SYSTEM_PROMPT = (
    'Réponds uniquement à partir du contexte fourni. '
    'Sépare clairement les faits observés dans le code et les hypothèses. '
    'Si la preuve est insuffisante, dis "non confirmé dans les sources indexées".'
)


@dataclass(frozen=True)
class ModelCatalog:
    dense_embedding_model: str = DEFAULT_EMBEDDING_MODEL
    code_embedding_model: str = DEFAULT_CODE_EMBEDDING_MODEL
    reranker_model: str = DEFAULT_RERANKER_MODEL


@dataclass(frozen=True)
class RetrievalWeights:
    lexical_weight: float = 1.0
    dense_weight: float = 0.35
    rerank_weight: float = 0.25
    fortran_priority_boost: float = 1.35
    python_priority_boost: float = 1.1
    docs_penalty: float = 0.9


@dataclass(frozen=True)
class RAGConfig:
    repo_root: Path = field(default_factory=lambda: Path(__file__).resolve().parents[3])
    chunk_size_lines: int = 120
    chunk_stride_lines: int = 90
    max_results: int = 8
    max_chunk_iterations: int = 100000
    retrieval_weights: RetrievalWeights = field(default_factory=RetrievalWeights)
    model_catalog: ModelCatalog = field(default_factory=ModelCatalog)
    system_prompt: str = DEFAULT_SYSTEM_PROMPT

    @property
    def docs_root(self) -> Path:
        return self.repo_root / 'docs'

    @property
    def fast_root(self) -> Path:
        return self.repo_root / 'Fast'
