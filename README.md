# DESeq2 Pipeline

Pseudobulk + DESeq2 pipeline for single-cell perturbation screens. It takes an AnnData (with raw counts, not log-transformed), splits cells into a user-specified number of pseudobulk replicates per perturbation, runs DESeq2 against control, and returns differential expression results.

---

## Setup (only need to run once)

Requires conda or mamba. No sudo needed. R and all dependencies install into the conda environment.

```bash
bash setup.sh
conda activate deseq2_pipeline
```

This `setup.sh` installs Python 3.11, R 4.3, and all required packages (DESeq2, glmGamPoi, etc.) as pre-built conda binaries — no compilation.

---

## Verify installation

```bash
pytest test_pipeline.py -v   # all 14 tests should pass
```

To speed things up, the main script employs parallelization, both across pseudobulking and DESeq implementation. This `test_pipeline.py` tests whether parallelization matches sequential results.

---

## Running the pipeline

### CLI

```bash
python deseq2_pipeline.py \
    --h5ad /path/to/data.h5ad \
    --pert-col target_gene \
    --ctrl-label non-targeting \
    --outdir results \
    --n-threads 4 \
    --n-workers-r 50
```

Test on a subset first with `--n-sample 50`. This would take a random sample of 50 perturbations and return results on that. Note the time, and then accordingly experiment with `n-threads` and `n-workers`. `n-threads` decides the number of threads for pseudobulking (in Python), while `n-workers-r` decides parallelization in R for actual DESeq2 implementation.

### Python API

```python
from deseq2_pipeline import run_pipeline

results_path = run_pipeline(
    h5ad_path   = "/path/to/data.h5ad",
    pert_col    = "target_gene",
    ctrl_label  = "non-targeting",
    n_threads   = 4,
    n_workers_r = 50,
)

# Optional: if you want to see the results
import pandas as pd
results = pd.read_csv(results_path)
```

---

## Output

Each run creates a unique timestamped directory under `--outdir` so reruns never overwrite each other:

```
results/
  run_20260221_143502_mydata/
    run_config.json          # all parameters used
    pseudobulk.csv.gz        # pseudobulk count matrix
    deseq2_results.csv.gz    # final results
    skipped_perts.csv        # perts excluded (too few cells)
```

### Result columns

| Column | Description |
|---|---|
| `gene` | Gene name |
| `perturbation` | Perturbation name |
| `baseMean` | Mean normalized counts across all samples |
| `log2FoldChange` | Log2 fold change (pert vs ctrl) |
| `fold_change` | Linear fold change (2^log2FoldChange) |
| `lfcSE` | Standard error of log2 fold change |
| `stat` | Wald statistic |
| `pvalue` | Nominal p-value |
| `fdr` | BH-adjusted p-value (Benjamini-Hochberg) |

---

## Reproducibility

**Pseudobulking is random** — control cells are sampled without replacement for each perturbation replicate. This means two runs on the same data will produce slightly different pseudobulk matrices and therefore slightly different DESeq2 results unless the same seed is used.

To fully reproduce a run, use the same `--random-seed` (default: 42) and the same parameters. All parameters are saved to `run_config.json` in the output directory for reference.

```bash
# reproduce a previous run exactly
python deseq2_pipeline.py \
    --h5ad /path/to/data.h5ad \
    --pert-col target_gene \
    --ctrl-label non-targeting \
    --random-seed 42          # must match original run
    --n-sample 100            # must match original run
```

---

## Key parameters

| Parameter | Default | Notes |
|---|---|---|
| `--n-threads` | 4 | Pseudobulk parallelism. 4 is the sweet spot — beyond 4 threads scipy sparse hits the GIL and gains plateau |
| `--n-workers-r` | 50 | Parallel R processes. Safe up to ~100 on a machine with 200GB RAM since adata is freed before R launches |
| `--min-cells` | 10 | Perts with fewer cells than this are skipped. Names of such skipped perturbations are logged to `skipped_perts.csv` |
| `--n-sample` | None | Subsample N random perts. Useful for quick tests before running all perts |
| `--random-seed` | 42 | Controls pseudobulk cell sampling. Fix this to reproduce results exactly |
| `--chunk-size` | 500 | Rows per chunk in sparse row sum. Keep at 500 — lower values avoid RSS spikes on systems with tight memory limits |
| `--n-reps` | 2 | Number of pseudobulk replicates per perturbation. Each replicate uses `rep-frac` of the pert's cells. Total cells used = `n-reps * rep-frac`|
| `--rep-frac` | 0.5 | Fraction of each pert's cells used per replicate. Default 0.5 means each replicate uses 50% of cells. Must satisfy `n-reps * rep-frac ≤ 1.0` — e.g. 3 replicates requires `rep-frac ≤ 0.33` |
