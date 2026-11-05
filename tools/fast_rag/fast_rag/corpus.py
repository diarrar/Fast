from __future__ import annotations

import ast
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .config import RAGConfig


OFFICIAL_SOURCE_LISTS = {
    'FastS': 'Fast/FastS/srcs.py',
    'FastC': 'Fast/FastC/srcs.py',
    'Fast': 'Fast/Fast/srcs.py',
}

PYTHON_GLOBS = (
    'Fast/FastS/FastS/**/*.py',
    'Fast/FastC/FastC/**/*.py',
    'Fast/Fast/Fast/**/*.py',
)

DOC_GLOBS = ('docs/**/*.html',)

CPP_SUFFIXES = {'.cpp', '.cxx', '.cc', '.c'}
FORTRAN_SUFFIXES = {'.for', '.f', '.f90', '.f95', '.F90', '.F'}


@dataclass(frozen=True)
class SourceDocument:
    path: Path
    module_name: str
    language: str
    priority: int
    source_kind: str


@dataclass(frozen=True)
class BenchmarkEntry:
    question: str
    expected_paths: tuple[str, ...]
    intent: str | None = None


def _literal_list(path: Path, variable_name: str) -> list[str]:
    tree = ast.parse(path.read_text(encoding='utf-8'))
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id == variable_name:
                    values: list[str] = []
                    for value in getattr(node.value, 'elts', []):
                        if isinstance(value, ast.Constant) and isinstance(value.value, str):
                            values.append(value.value)
                        elif isinstance(value, ast.Str):
                            values.append(value.s)
                    return values
    return []


def _language_from_path(path: Path) -> str:
    suffix = path.suffix
    if suffix in FORTRAN_SUFFIXES:
        return 'fortran'
    if suffix in CPP_SUFFIXES:
        return 'cpp'
    if suffix == '.py':
        return 'python'
    if suffix == '.html':
        return 'html'
    return suffix.lstrip('.') or 'text'


def _priority_for_language(language: str) -> int:
    if language == 'fortran':
        return 1
    if language == 'python':
        return 2
    if language == 'html':
        return 3
    if language == 'cpp':
        return 2
    return 4


def _module_from_path(relative_path: Path) -> str:
    parts = relative_path.parts
    for candidate in ('FastS', 'FastC', 'Fast'):
        if candidate in parts:
            return candidate
    return 'Fast'


def _official_sources(repo_root: Path) -> Iterable[SourceDocument]:
    for module_name, relative_srcs in OFFICIAL_SOURCE_LISTS.items():
        srcs_path = repo_root / relative_srcs
        source_root = srcs_path.parent
        for variable_name in ('for_srcs', 'cpp_srcs'):
            for relative_file in _literal_list(srcs_path, variable_name):
                absolute_path = source_root / relative_file
                if absolute_path.exists():
                    language = _language_from_path(absolute_path)
                    yield SourceDocument(
                        path=absolute_path,
                        module_name=module_name,
                        language=language,
                        priority=_priority_for_language(language),
                        source_kind='srcs.py',
                    )


def _glob_sources(config: RAGConfig, patterns: Iterable[str]) -> Iterable[SourceDocument]:
    seen: set[Path] = set()
    for pattern in patterns:
        for path in config.repo_root.glob(pattern):
            if path.is_file() and path not in seen:
                seen.add(path)
                relative_path = path.relative_to(config.repo_root)
                language = _language_from_path(path)
                yield SourceDocument(
                    path=path,
                    module_name=_module_from_path(relative_path),
                    language=language,
                    priority=_priority_for_language(language),
                    source_kind='glob',
                )


def discover_corpus(config: RAGConfig) -> list[SourceDocument]:
    documents: dict[Path, SourceDocument] = {}
    for document in _official_sources(config.repo_root):
        documents[document.path] = document
    for document in _glob_sources(config, PYTHON_GLOBS + DOC_GLOBS):
        documents.setdefault(document.path, document)
    return sorted(documents.values(), key=lambda item: (item.priority, item.module_name, str(item.path)))


def load_benchmark_entries(path: Path) -> list[BenchmarkEntry]:
    entries: list[BenchmarkEntry] = []
    for raw_line in path.read_text(encoding='utf-8').splitlines():
        line = raw_line.strip()
        if not line:
            continue
        payload = json.loads(line)
        entries.append(
            BenchmarkEntry(
                question=payload['question'],
                expected_paths=tuple(payload.get('expected_paths', [])),
                intent=payload.get('intent'),
            )
        )
    return entries
