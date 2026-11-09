#!/usr/bin/env python3
"""Compare ordered Virtual Embryo board panels without opening expression data."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import re
import sys
from datetime import datetime, timezone
from urllib.request import Request, urlopen

BASE = 'https://virtualembryo.ai/challenge/panels/'
MAX_BYTES = 2_000_000


def digest(data):
    return hashlib.sha256(data).hexdigest()


def labels(data):
    text = data.decode('utf-8-sig')
    result = text.splitlines()
    if not result or any(not x or x.strip() != x for x in result):
        raise ValueError('Gene list must contain one unpadded, nonempty name per line')
    if len(set(result)) != len(result):
        raise ValueError('Duplicate gene names make an index mapping ambiguous')
    return result


def read_index(cache):
    index = json.loads((cache / 'index.json').read_text())
    if not isinstance(index, dict):
        raise ValueError('Board index must be an object')
    return index


def board_spec(index, key):
    if not isinstance(index, dict):
        raise ValueError('Board index must be an object')
    if key not in index:
        raise ValueError('Unknown board: ' + key)
    spec = index[key]
    if not isinstance(spec, dict):
        raise ValueError('Board specification must be an object')
    filename = spec.get('genes_file', '')
    if not re.fullmatch(r'[A-Za-z0-9_]+\.genes\.txt', filename):
        raise ValueError('Unexpected panel filename')
    for field in ('n_genes', 'min_cells', 'max_cells'):
        if type(spec.get(field)) is not int or spec[field] <= 0:
            raise ValueError('Invalid board field: ' + field)
    if spec['min_cells'] > spec['max_cells']:
        raise ValueError('Invalid cell range')
    if not isinstance(spec.get('obsm_required'), list) or not all(isinstance(x, str) for x in spec['obsm_required']):
        raise ValueError('Invalid required coordinate-key list')
    return spec


def check_panel(data, spec):
    names = labels(data)
    if len(names) != spec['n_genes']:
        raise ValueError('Panel count differs from index')
    expected = spec.get('genes_sha256', '')
    if not isinstance(expected, str) or not re.fullmatch(r'[0-9a-f]{16}', expected):
        raise ValueError('Expected published 16-character panel digest')
    # The three real downloaded panels match the digest of names joined with LF,
    # without a trailing LF; raw HTTP bytes have a separate receipt digest.
    if digest('\n'.join(names).encode())[:16] != expected:
        raise ValueError('Ordered panel digest differs from published index')
    return names


def panel(cache, index, key):
    spec = board_spec(index, key)
    data = (cache / spec['genes_file']).read_bytes()
    names = check_panel(data, spec)
    receipt = json.loads((cache / 'receipt.json').read_text())
    if not isinstance(receipt, dict) or not isinstance(receipt.get('files'), dict):
        raise ValueError('Snapshot receipt must contain a file-digest object')
    expected = receipt['files'].get(spec['genes_file'])
    if expected != digest(data) or receipt['files'].get('index.json') != digest((cache / 'index.json').read_bytes()):
        raise ValueError('Cached snapshot changed since acquisition; fetch a new directory')
    return names, spec, digest(data)


def download(name):
    req = Request(BASE + name, headers={'User-Agent': 'VirtualEmbryoPanelBridge/1.0 (public-panel-metadata-only)'})
    with urlopen(req, timeout=30) as response:
        data = response.read(MAX_BYTES + 1)
    if len(data) > MAX_BYTES:
        raise ValueError('Unexpectedly large metadata response')
    return data


def fetch(cache, keys):
    if cache.exists():
        raise ValueError('Snapshot destination already exists; use a new directory')
    raw = download('index.json')
    index = json.loads(raw)
    files = {'index.json': raw}
    for key in keys:
        spec = board_spec(index, key)
        name = spec['genes_file']
        if name not in files:
            files[name] = download(name)
        check_panel(files[name], spec)
    # Only commit a snapshot once all requested public metadata has passed checks.
    cache.mkdir(parents=True)
    for name, data in files.items():
        (cache / name).write_bytes(data)
    receipt = {'retrieved_utc': datetime.now(timezone.utc).isoformat(), 'base_url': BASE,
               'boards': keys, 'files': {name: digest(data) for name, data in files.items()},
               'bytes': sum(map(len, files.values())),
               'scope': 'Public submission-axis metadata only; not expression, labels or a data-use approval'}
    (cache / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    return receipt


def compare(source, target):
    if len(set(source)) != len(source) or len(set(target)) != len(target):
        raise ValueError('Duplicate names are ambiguous')
    lookup = {gene: i for i, gene in enumerate(source)}
    target_set = set(target)
    missing = [gene for gene in target if gene not in lookup]
    dropped = [gene for gene in source if gene not in target_set]
    indices = None if missing else [lookup[gene] for gene in target]
    return {'source_genes': len(source), 'target_genes': len(target),
            'exact_match': source == target, 'can_project_existing_columns': not missing,
            'missing_from_source': missing, 'excluded_source_genes': dropped,
            'source_column_indices_in_target_order': indices,
            'order_changes': False if indices is None else indices != sorted(indices)}


def audit(cache, source_board, source_genes, target_board, cells):
    index = read_index(cache)
    target, spec, target_hash = panel(cache, index, target_board)
    if source_board:
        source, source_spec, source_hash = panel(cache, index, source_board)
        origin = {'board': source_board, 'sha256': source_hash}
    else:
        data = source_genes.read_bytes()
        source = labels(data)
        origin = {'local_gene_list': str(source_genes), 'sha256': digest(data)}
    result = compare(source, target)
    result.update(source=origin, target_board=target_board, target_panel_sha256=target_hash,
                  required_obsm=spec['obsm_required'], accepted_cells=[spec['min_cells'], spec['max_cells']],
                  scope='Column-axis planning only; not a valid prediction, score, biological transfer, or permission to reuse data')
    if cells is not None:
        if cells <= 0:
            raise ValueError('Cell count must be positive')
        result['requested_cells'] = cells
        result['cell_count_in_range'] = spec['min_cells'] <= cells <= spec['max_cells']
        result['dense_float32_bytes_arithmetic_only'] = cells * len(target) * 4
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch', help='Fetch small public panel metadata to a NEW local snapshot directory')
    f.add_argument('--cache', type=Path, required=True)
    f.add_argument('--boards', nargs='+', required=True)
    a = sub.add_parser('audit', help='Plan axis projection; never reads or rewrites expression')
    a.add_argument('--cache', type=Path, required=True)
    source = a.add_mutually_exclusive_group(required=True)
    source.add_argument('--source-board')
    source.add_argument('--source-genes', type=Path)
    a.add_argument('--target-board', required=True)
    a.add_argument('--cells', type=int)
    a.add_argument('--output', type=Path, help='New JSON report path; refuses overwrite')
    args = parser.parse_args()
    try:
        if args.command == 'fetch':
            result = fetch(args.cache, args.boards)
        else:
            result = audit(args.cache, args.source_board, args.source_genes, args.target_board, args.cells)
            if args.output:
                with args.output.open('x') as f:
                    f.write(json.dumps(result, indent=2) + '\n')
        print(json.dumps(result, indent=2))
        if args.command == 'audit' and (not result['can_project_existing_columns'] or result.get('cell_count_in_range') is False):
            return 2
        return 0
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(json.dumps({'error': str(exc)}), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
