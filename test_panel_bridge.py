import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import panel_bridge as pb


class PanelBridgeTests(unittest.TestCase):
    def fixture(self, root):
        source = b'A\nB\nC\nD\n'
        target = b'D\nB\nA\n'
        specs = {}
        for key, filename, raw in [('source', 'source.genes.txt', source), ('target', 'target.genes.txt', target)]:
            specs[key] = dict(genes_file=filename, n_genes=len(raw.splitlines()),
                genes_sha256=pb.digest('\n'.join(raw.decode().splitlines()).encode())[:16],
                min_cells=2, max_cells=5, obsm_required=['spatial_3D'])
            (root / filename).write_bytes(raw)
        (root / 'index.json').write_text(json.dumps(specs))
        files = {f.name: pb.digest(f.read_bytes()) for f in root.iterdir()}
        (root / 'receipt.json').write_text(json.dumps({'files': files}))
        return specs

    def test_exact_order_projection_and_explicit_drop(self):
        report = pb.compare(['A', 'B', 'C', 'D'], ['D', 'B', 'A'])
        self.assertEqual(report['source_column_indices_in_target_order'], [3, 1, 0])
        self.assertEqual(report['excluded_source_genes'], ['C'])
        self.assertTrue(report['order_changes'])
        self.assertFalse(report['exact_match'])

    def test_missing_never_produces_a_partial_mapping(self):
        report = pb.compare(['A'], ['A', 'Z'])
        self.assertEqual(report['missing_from_source'], ['Z'])
        self.assertIsNone(report['source_column_indices_in_target_order'])
        self.assertFalse(report['can_project_existing_columns'])

    def test_duplicates_whitespace_and_empty_rejected(self):
        for raw in [b'A\nA\n', b'A\n B\n', b'A\n\n', b'']:
            with self.assertRaises(ValueError): pb.labels(raw)
        with self.assertRaises(ValueError): pb.compare(['A', 'A'], ['A'])
        self.assertEqual(pb.labels(b'A\r\nB\r\n'), ['A', 'B'])

    def test_case_sensitive_no_symbol_guessing(self):
        self.assertEqual(pb.compare(['gata4'], ['Gata4'])['missing_from_source'], ['Gata4'])

    def test_offline_audit_cell_edges_and_immutability(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); self.fixture(root)
            before = {p.name: p.read_bytes() for p in root.iterdir()}
            for cells in [2, 5]:
                r = pb.audit(root, 'source', None, 'target', cells)
                self.assertTrue(r['cell_count_in_range'])
                self.assertEqual(r['dense_float32_bytes_arithmetic_only'], cells * 3 * 4)
            for cells in [1, 6]:
                self.assertFalse(pb.audit(root, 'source', None, 'target', cells)['cell_count_in_range'])
            self.assertEqual(before, {p.name: p.read_bytes() for p in root.iterdir()})

    def test_published_digest_and_receipt_detect_different_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); specs = self.fixture(root)
            specs['target']['genes_sha256'] = '0' * 16
            with self.assertRaisesRegex(ValueError, 'published index'):
                pb.panel(root, specs, 'target')
            self.fixture(root)
            raw = (root / 'target.genes.txt').read_bytes().replace(b'\n', b'\r\n')
            (root / 'target.genes.txt').write_bytes(raw)
            with self.assertRaisesRegex(ValueError, 'snapshot changed'):
                pb.panel(root, pb.read_index(root), 'target')

    def test_malformed_metadata_is_a_value_error(self):
        for index in [[], {'x': []}, {'x': {'genes_file': '../file'}}]:
            with self.assertRaises(ValueError): pb.board_spec(index, 'x')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); specs = self.fixture(root)
            (root / 'receipt.json').write_text('{"files": []}')
            with self.assertRaises(ValueError): pb.panel(root, specs, 'target')

    def test_fetch_checks_every_panel_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); specs = self.fixture(root)
            destination = root / 'snapshot'
            index_raw = json.dumps(specs).encode()
            def fake(name):
                return index_raw if name == 'index.json' else b'WRONG\n'
            with patch.object(pb, 'download', side_effect=fake):
                with self.assertRaises(ValueError): pb.fetch(destination, ['target'])
            self.assertFalse(destination.exists())
            with self.assertRaises(ValueError): pb.fetch(root, ['target'])


if __name__ == '__main__': unittest.main()
