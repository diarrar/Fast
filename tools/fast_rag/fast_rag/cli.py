from __future__ import annotations

import argparse
import json
from pathlib import Path

from .config import RAGConfig
from .rag import FastRAGPipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Fast RAG tooling focused on Fortran-first retrieval.')
    parser.add_argument('--repo-root', type=Path, default=RAGConfig().repo_root)
    subparsers = parser.add_subparsers(dest='command', required=True)

    manifest = subparsers.add_parser('manifest', help='Summarize the discovered corpus.')
    manifest.set_defaults(command='manifest')

    search = subparsers.add_parser('search', help='Search the hybrid index.')
    search.add_argument('query')
    search.add_argument('--top-k', type=int, default=8)

    prompt = subparsers.add_parser('prompt', help='Build a citation-oriented prompt payload.')
    prompt.add_argument('query')
    prompt.add_argument('--top-k', type=int, default=8)

    benchmark = subparsers.add_parser('benchmark', help='Evaluate recall@k from a JSONL benchmark file.')
    benchmark.add_argument('path', type=Path)
    benchmark.add_argument('--top-k', type=int, default=5)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    pipeline = FastRAGPipeline(RAGConfig(repo_root=args.repo_root))

    if args.command == 'manifest':
        counts: dict[str, int] = {}
        for document in pipeline.documents:
            key = f'{document.module_name}:{document.language}'
            counts[key] = counts.get(key, 0) + 1
        print(json.dumps(counts, indent=2, sort_keys=True))
        return 0

    if args.command == 'search':
        results = pipeline.search(args.query, top_k=args.top_k)
        payload = [
            {
                'score': round(result.score, 6),
                'path': str(result.chunk.path),
                'language': result.chunk.language,
                'module': result.chunk.module_name,
                'title': result.chunk.title,
                'line_start': result.chunk.line_start,
                'line_end': result.chunk.line_end,
                'intent': result.intent,
                'calls': result.chunk.calls,
            }
            for result in results
        ]
        print(json.dumps(payload, indent=2))
        return 0

    if args.command == 'prompt':
        payload = pipeline.build_prompt_payload(args.query, top_k=args.top_k)
        print(json.dumps({
            'intent': payload.intent,
            'query_variants': payload.query_variants,
            'system_prompt': payload.system_prompt,
            'context': payload.context,
            'citations': payload.citations,
        }, indent=2))
        return 0

    if args.command == 'benchmark':
        metrics = pipeline.benchmark_recall_at_k(args.path, top_k=args.top_k)
        print(json.dumps(metrics, indent=2))
        return 0

    parser.error(f'Unsupported command: {args.command}')
    return 2


if __name__ == '__main__':
    raise SystemExit(main())
