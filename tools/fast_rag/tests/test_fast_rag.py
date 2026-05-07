from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from tools.fast_rag.fast_rag.config import RAGConfig
from tools.fast_rag.fast_rag.corpus import BenchmarkEntry, SourceDocument, discover_corpus, load_benchmark_entries
from tools.fast_rag.fast_rag.fortran import chunk_fortran_source
from tools.fast_rag.fast_rag.rag import FastRAGPipeline, classify_intent


REPO_ROOT = Path(__file__).resolve().parents[3]


class FastRAGTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = RAGConfig(repo_root=REPO_ROOT)

    def test_discover_corpus_uses_official_fortran_lists(self):
        documents = discover_corpus(self.config)
        fasts_farfield = REPO_ROOT / 'Fast/FastS/FastS/BC/bvbs_farfield.for'
        self.assertIn(fasts_farfield, {document.path for document in documents})
        selected = next(document for document in documents if document.path == fasts_farfield)
        self.assertEqual(selected.language, 'fortran')
        self.assertEqual(selected.source_kind, 'srcs.py')

    def test_chunk_fortran_source_extracts_routine_and_calls(self):
        document = SourceDocument(
            path=REPO_ROOT / 'Fast/FastC/FastC/Metric/skmtr.for',
            module_name='FastC',
            language='fortran',
            priority=1,
            source_kind='srcs.py',
        )
        chunks = chunk_fortran_source(document, self.config)
        self.assertTrue(chunks)
        self.assertEqual(chunks[0].routine_name, 'skmtr')
        self.assertIn('cp_tijk', chunks[0].calls)
        self.assertEqual(chunks[0].family, 'Metric')

    def test_search_prioritizes_fortran_for_solver_queries(self):
        pipeline = FastRAGPipeline(self.config)
        results = pipeline.search('farfield boundary condition routine', top_k=5)
        self.assertTrue(results)
        self.assertEqual(results[0].chunk.language, 'fortran')
        self.assertIn('bvbs_farfield.for', str(results[0].chunk.path))

    def test_prompt_payload_orders_fortran_before_docs(self):
        pipeline = FastRAGPipeline(self.config)
        payload = pipeline.build_prompt_payload('farfield boundary condition', top_k=5)
        self.assertIn('FORTRAN', payload.context)
        if 'HTML' in payload.context:
            self.assertLess(payload.context.index('FORTRAN'), payload.context.index('HTML'))

    def test_benchmark_loading_and_recall(self):
        with tempfile.NamedTemporaryFile('w+', suffix='.jsonl', delete=False) as handle:
            path = Path(handle.name)
            handle.write(json.dumps({
                'question': 'où est calculée la métrique ?',
                'expected_paths': ['Fast/FastC/FastC/Metric/skmtr.for'],
                'intent': 'noyau_numerique',
            }) + '\n')
        try:
            entries = load_benchmark_entries(path)
            self.assertEqual(entries, [BenchmarkEntry(
                question='où est calculée la métrique ?',
                expected_paths=('Fast/FastC/FastC/Metric/skmtr.for',),
                intent='noyau_numerique',
            )])
            metrics = FastRAGPipeline(self.config).benchmark_recall_at_k(path, top_k=5)
            self.assertEqual(metrics['questions'], 1)
            self.assertGreaterEqual(metrics['recall_at_k'], 1.0)
        finally:
            path.unlink(missing_ok=True)

    def test_intent_classifier(self):
        self.assertEqual(classify_intent('show python api wrapper for warmup'), 'api_python')
        self.assertEqual(classify_intent('gpu build flags'), 'build_perf_gpu')
        self.assertEqual(classify_intent('validation benchmark'), 'validation')
        self.assertEqual(classify_intent('roe flux routine'), 'noyau_numerique')


if __name__ == '__main__':
    unittest.main()
