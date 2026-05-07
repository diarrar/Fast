from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .config import RAGConfig
from .corpus import SourceDocument


ROUTINE_START_RE = re.compile(r'^\s*(subroutine|function)\s+([a-z_][a-z0-9_]*)', re.IGNORECASE)
ROUTINE_END_RE = re.compile(r'^\s*end(\s+(subroutine|function))?\b', re.IGNORECASE)
CALL_RE = re.compile(r'\bcall\s+([a-z_][a-z0-9_]*)', re.IGNORECASE)
INCLUDE_RE = re.compile(r'^(?:\s*#include|\s*include)\s+["\']([^"\']+)["\']', re.IGNORECASE)
COMMENT_RE = re.compile(r'^\s*[c*!]')


@dataclass(frozen=True)
class Chunk:
    chunk_id: str
    path: Path
    module_name: str
    language: str
    priority: int
    text: str
    title: str
    line_start: int
    line_end: int
    routine_name: str | None
    family: str
    calls: tuple[str, ...]
    includes: tuple[str, ...]
    source_kind: str


def _family_for_path(path: Path) -> str:
    families = ('BC', 'ROE', 'LES', 'POST', 'ADJOINT', 'Metric', 'Compute', 'FILTER', 'INTERP', 'STAT', 'Init', 'ALE')
    parts = {part.upper(): part for part in path.parts}
    for family in families:
        if family.upper() in parts:
            return family
    return path.parent.name or 'root'


def _window_chunks(document: SourceDocument, lines: list[str], config: RAGConfig) -> list[Chunk]:
    chunks: list[Chunk] = []
    start = 0
    chunk_index = 0
    while start < len(lines):
        end = min(start + config.chunk_size_lines, len(lines))
        text = '\n'.join(lines[start:end]).strip()
        if text:
            chunks.append(
                Chunk(
                    chunk_id=f'{document.path}:{start + 1}-{end}',
                    path=document.path,
                    module_name=document.module_name,
                    language=document.language,
                    priority=document.priority,
                    text=text,
                    title=f'{document.path.name}:{start + 1}-{end}',
                    line_start=start + 1,
                    line_end=end,
                    routine_name=None,
                    family=_family_for_path(document.path),
                    calls=tuple(sorted({match.group(1).lower() for match in CALL_RE.finditer(text)})),
                    includes=tuple(sorted({match.group(1) for match in INCLUDE_RE.finditer(text)})),
                    source_kind=document.source_kind,
                )
            )
        if end == len(lines):
            break
        start += config.chunk_stride_lines
        chunk_index += 1
        if chunk_index > 100000:
            break
    return chunks


def chunk_fortran_source(document: SourceDocument, config: RAGConfig) -> list[Chunk]:
    text = document.path.read_text(encoding='utf-8', errors='ignore')
    lines = text.splitlines()
    chunks: list[Chunk] = []
    routine_name: str | None = None
    routine_start = 0
    buffer: list[str] = []
    for line_number, line in enumerate(lines, start=1):
        if routine_name is None:
            match = ROUTINE_START_RE.match(line)
            if match:
                routine_name = match.group(2).lower()
                routine_start = line_number
                buffer = [line]
            continue
        buffer.append(line)
        if ROUTINE_END_RE.match(line):
            routine_text = '\n'.join(buffer).strip()
            chunks.append(
                Chunk(
                    chunk_id=f'{document.path}:{routine_name}:{routine_start}-{line_number}',
                    path=document.path,
                    module_name=document.module_name,
                    language=document.language,
                    priority=document.priority,
                    text=routine_text,
                    title=f'{routine_name} [{document.path.name}]',
                    line_start=routine_start,
                    line_end=line_number,
                    routine_name=routine_name,
                    family=_family_for_path(document.path),
                    calls=tuple(sorted({match.group(1).lower() for match in CALL_RE.finditer(routine_text)})),
                    includes=tuple(sorted({match.group(1) for match in INCLUDE_RE.finditer(routine_text)})),
                    source_kind=document.source_kind,
                )
            )
            routine_name = None
            buffer = []
    if chunks:
        return chunks
    return _window_chunks(document, lines, config)


def chunk_text_source(document: SourceDocument, config: RAGConfig) -> list[Chunk]:
    lines = document.path.read_text(encoding='utf-8', errors='ignore').splitlines()
    return _window_chunks(document, lines, config)
