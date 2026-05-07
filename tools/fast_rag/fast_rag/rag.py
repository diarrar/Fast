from __future__ import annotations

import math
import re
import unicodedata
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

from .config import RAGConfig
from .corpus import BenchmarkEntry, discover_corpus, load_benchmark_entries
from .fortran import Chunk, chunk_fortran_source, chunk_text_source


TOKEN_RE = re.compile(r'[A-Za-z_][A-Za-z0-9_./+-]*')


class DenseScorer(Protocol):
    def score(self, query: str, chunk: Chunk) -> float:
        ...


class Reranker(Protocol):
    def score(self, query: str, chunk: Chunk) -> float:
        ...


class NullDenseScorer:
    def __init__(self, model_name: str):
        self.model_name = model_name

    def score(self, query: str, chunk: Chunk) -> float:
        return 0.0


class HeuristicReranker:
    def __init__(self, model_name: str):
        self.model_name = model_name

    def score(self, query: str, chunk: Chunk) -> float:
        query_tokens = set(tokenize(query))
        path_tokens = set(tokenize(str(chunk.path)))
        routine_tokens = set(tokenize(chunk.routine_name or ''))
        return 0.6 * overlap_score(query_tokens, path_tokens) + 0.4 * overlap_score(query_tokens, routine_tokens)


@dataclass(frozen=True)
class SearchResult:
    chunk: Chunk
    score: float
    lexical_score: float
    dense_score: float
    rerank_score: float
    intent: str


@dataclass(frozen=True)
class PromptPayload:
    intent: str
    query_variants: tuple[str, ...]
    system_prompt: str
    context: str
    citations: tuple[str, ...]


def tokenize(text: str) -> list[str]:
    normalized = unicodedata.normalize('NFKD', text)
    ascii_text = normalized.encode('ascii', 'ignore').decode('ascii')
    return [token.lower() for token in TOKEN_RE.findall(ascii_text)]


def overlap_score(left: set[str], right: set[str]) -> float:
    if not left or not right:
        return 0.0
    return len(left & right) / len(left)


def classify_intent(query: str) -> str:
    lowered = query.lower()
    if any(token in lowered for token in ('pytree', 'api', 'python', 'module', 'wrapper')):
        return 'api_python'
    if any(token in lowered for token in ('gpu', 'openmp', 'perf', 'performance', 'build', 'scons', 'compile')):
        return 'build_perf_gpu'
    if any(token in lowered for token in ('test', 'validation', 'benchmark', 'accuracy', 'recall')):
        return 'validation'
    return 'noyau_numerique'


def expand_query(query: str, intent: str) -> tuple[str, ...]:
    variants = [query.strip()]
    normalized = ' '.join(tokenize(query))
    synonym_hints = {
        'metrique': 'metric skmtr',
        'metric': 'metric skmtr',
        'farfield': 'bc farfield bvbs_farfield',
        'boundary': 'bc boundary condition',
        'condition': 'bc boundary condition',
        'gpu': 'gpu openmp offload',
    }
    for keyword, hint in synonym_hints.items():
        if keyword in normalized:
            variants.append(f'{hint} {query}')
    if intent == 'api_python':
        variants.append(f'python wrapper {query}')
    elif intent == 'build_perf_gpu':
        variants.append(f'fortran gpu build {query}')
    elif intent == 'validation':
        variants.append(f'benchmark evaluation {query}')
    else:
        variants.append(f'fortran routine {query}')
        variants.append(f'cfd solver {query}')
    return tuple(dict.fromkeys(variant for variant in variants if variant))


def lexical_score(query: str, chunk: Chunk, inverse_document_frequency: dict[str, float]) -> float:
    query_tokens = tokenize(query)
    chunk_tokens = tokenize(chunk.text)
    if not query_tokens or not chunk_tokens:
        return 0.0
    frequencies = Counter(chunk_tokens)
    score = 0.0
    for token in query_tokens:
        score += frequencies[token] * inverse_document_frequency.get(token, 1.0)
    normalization = (len(chunk_tokens) ** 0.75) or 1.0
    return score / normalization


def metadata_score(query: str, chunk: Chunk, inverse_document_frequency: dict[str, float]) -> float:
    query_tokens = set(tokenize(query))
    metadata_tokens = set(tokenize(
        f'{chunk.title} {chunk.path} {chunk.family} {chunk.routine_name or ""} '
        f'{" ".join(chunk.calls)} {" ".join(chunk.includes)}'
    ))
    if not query_tokens or not metadata_tokens:
        return 0.0
    score = sum(
        inverse_document_frequency.get(token, 1.0)
        for token in query_tokens
        if token in metadata_tokens
    )
    path_tokens = set(tokenize(chunk.path.name))
    routine_tokens = set(tokenize(chunk.routine_name or ''))
    score += sum(
        1.5 * inverse_document_frequency.get(token, 1.0)
        for token in query_tokens
        if token in path_tokens or token in routine_tokens
    )
    return score / len(query_tokens)


class FastRAGPipeline:
    def __init__(
        self,
        config: RAGConfig | None = None,
        dense_scorer: DenseScorer | None = None,
        reranker: Reranker | None = None,
    ):
        self.config = config or RAGConfig()
        self.dense_scorer = dense_scorer or NullDenseScorer(self.config.model_catalog.code_embedding_model)
        self.reranker = reranker or HeuristicReranker(self.config.model_catalog.reranker_model)
        self.documents = discover_corpus(self.config)
        self.chunks = self._build_chunks()
        self.inverse_document_frequency = self._build_idf()

    def _build_chunks(self) -> list[Chunk]:
        chunks: list[Chunk] = []
        for document in self.documents:
            if document.language == 'fortran':
                chunks.extend(chunk_fortran_source(document, self.config))
            else:
                chunks.extend(chunk_text_source(document, self.config))
        return chunks

    def _build_idf(self) -> dict[str, float]:
        documents_by_token: defaultdict[str, int] = defaultdict(int)
        for chunk in self.chunks:
            for token in set(tokenize(chunk.text)):
                documents_by_token[token] += 1
        total = max(len(self.chunks), 1)
        return {
            token: math.log((1 + total) / (1 + frequency)) + 1.0
            for token, frequency in documents_by_token.items()
        }

    def _language_boost(self, chunk: Chunk, intent: str) -> float:
        weights = self.config.retrieval_weights
        boost = 1.0
        if chunk.language == 'fortran':
            boost *= weights.fortran_priority_boost
            if intent in ('noyau_numerique', 'build_perf_gpu'):
                boost *= 1.15
        elif chunk.language == 'python':
            boost *= weights.python_priority_boost
            if intent == 'api_python':
                boost *= 1.15
        elif chunk.language == 'html':
            boost *= weights.docs_penalty
        elif chunk.language == 'cpp' and intent == 'noyau_numerique':
            boost *= 1.1
        if chunk.source_kind == 'srcs.py':
            boost *= 1.05
        return boost / chunk.priority

    def search(self, query: str, top_k: int | None = None) -> list[SearchResult]:
        limit = top_k or self.config.max_results
        intent = classify_intent(query)
        query_variants = expand_query(query, intent)
        scored: list[SearchResult] = []
        weights = self.config.retrieval_weights
        for chunk in self.chunks:
            lexical = max(lexical_score(variant, chunk, self.inverse_document_frequency) for variant in query_variants)
            lexical += 2.5 * max(metadata_score(variant, chunk, self.inverse_document_frequency) for variant in query_variants)
            dense = max(self.dense_scorer.score(variant, chunk) for variant in query_variants)
            rerank = self.reranker.score(query, chunk)
            base_score = (
                lexical * weights.lexical_weight +
                dense * weights.dense_weight +
                rerank * weights.rerank_weight
            )
            score = base_score * self._language_boost(chunk, intent)
            if score > 0:
                scored.append(
                    SearchResult(
                        chunk=chunk,
                        score=score,
                        lexical_score=lexical,
                        dense_score=dense,
                        rerank_score=rerank,
                        intent=intent,
                    )
                )
        scored.sort(key=lambda item: item.score, reverse=True)
        return scored[:limit]

    def build_prompt_payload(self, query: str, top_k: int | None = None) -> PromptPayload:
        intent = classify_intent(query)
        query_variants = expand_query(query, intent)
        results = self.search(query, top_k=top_k)
        ordered_languages = ('fortran', 'cpp', 'python', 'html')
        grouped: defaultdict[str, list[SearchResult]] = defaultdict(list)
        for result in results:
            grouped[result.chunk.language].append(result)
        sections: list[str] = []
        citations: list[str] = []
        for language in ordered_languages:
            if not grouped.get(language):
                continue
            lines = [language.upper()]
            for result in grouped[language]:
                chunk = result.chunk
                citation = f'{chunk.path}:{chunk.line_start}-{chunk.line_end}'
                citations.append(citation)
                lines.append(f'- {citation} :: {chunk.title}')
                lines.append(chunk.text)
            sections.append('\n'.join(lines))
        system_prompt = (
            'Réponds uniquement à partir du contexte fourni. '
            'Sépare clairement les faits observés dans le code et les hypothèses. '
            'Si la preuve est insuffisante, dis "non confirmé dans les sources indexées".'
        )
        return PromptPayload(
            intent=intent,
            query_variants=query_variants,
            system_prompt=system_prompt,
            context='\n\n'.join(sections),
            citations=tuple(citations),
        )

    def benchmark_recall_at_k(self, benchmark_path: Path, top_k: int = 5) -> dict[str, float]:
        entries = load_benchmark_entries(benchmark_path)
        if not entries:
            return {'recall_at_k': 0.0, 'questions': 0}
        hits = 0
        for entry in entries:
            results = self.search(entry.question, top_k=top_k)
            returned_paths = {str(result.chunk.path) for result in results}
            expected_paths = {
                str(Path(expected))
                if Path(expected).is_absolute()
                else str((self.config.repo_root / expected).resolve())
                for expected in entry.expected_paths
            }
            if any(expected in returned_paths for expected in expected_paths):
                hits += 1
        return {'recall_at_k': hits / len(entries), 'questions': len(entries)}
