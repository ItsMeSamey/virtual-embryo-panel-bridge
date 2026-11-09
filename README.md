# Panel Bridge

Panel Bridge checks an ordered gene list against a public Virtual Embryo panel before inference. It returns source-column indices in target order, lists excluded or missing genes, verifies the panel digest, and reports cell bounds and required coordinate keys.

Use Python 3.10 or newer. The tool uses the standard library.

## Fetch the panel specifications

```sh
python panel_bridge.py fetch --cache .cache/example \
  --boards T2:embryo:val_interp T2:embryo:test_interp T2:heart:test_extrap
```

Use a new cache directory for each fetch. Each directory keeps a fixed snapshot for later audits. The [official data page](https://virtualembryo.ai/challenge/data) and [panel index](https://virtualembryo.ai/challenge/panels/index.json) describe the board specifications.

## Check an axis

Export your model's gene names in their current order, one per line. Preserve spelling and case. Then run:

```sh
python panel_bridge.py audit --cache .cache/example \
  --source-genes model_genes.txt --target-board T2:embryo:test_interp \
  --cells 2000 --output model-axis-report.json
```

Choose a new output path for each report. To compare two public panels, replace `--source-genes model_genes.txt` with `--source-board T2:embryo:val_interp`.

Read `source_column_indices_in_target_order` before selecting prediction columns. Missing target genes set this mapping to `null`. Excluded source genes are listed separately. Duplicate, blank or padded names produce an error.

`order_changes` reports whether retained columns need a new order. `exact_match` also requires every gene to be present on both axes. `dense_float32_bytes_arithmetic_only` reports `cells * genes * 4`; size the full inference process separately.

## Exit codes

For a completed audit:

- `0`: report created; all target genes are available and the optional cell count is in range.
- `2`: report created; target genes are missing or the cell count is out of range.

A successful fetch returns `0`. Invalid command-line arguments return `2` with usage text. Input, fetch or snapshot processing errors return `1` with a JSON error on stderr.

## Tests

```sh
python -m unittest -v test_panel_bridge
```

Eight synthetic tests cover ordering, exclusion, missing genes, duplicate labels, case sensitivity, cell boundaries, fixed snapshots, digest checks and malformed metadata.

## Licence and credits

Original code and synthetic tests use the MIT licence. Official panel specifications are published by the Virtual Embryo Challenge, Qiu Lab at Stanford, and remain subject to the [Challenge Terms](https://virtualembryo.ai/challenge/terms). The data-generating team is Dr. Neil Chi's group at UCSD. Downloaded specifications stay in the chosen cache directory.
