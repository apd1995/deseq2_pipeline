# deseq2_pipeline

Pseudobulk + DESeq2 pipeline for single-cell perturbation screens. Takes an AnnData, splits cells into pseudobulk replicates per perturbation, runs DESeq2 against control, and returns differential expression results.

---

## Setup (run once)

Requires conda or mamba. No sudo needed — R and all dependencies install into the conda environment.

```bash
bash setup.sh
conda activate deseq2_pipeline
```

That's it. `setup.sh` installs Python 3.11, R 4.3, and all required packages (DESeq2, glmGamPoi, etc.) as pre-built conda binaries — no compilation.

---

## Verify installation

```bash
pytest test_pipeline.py -v   # all 14 tests should pass
```

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

Test on a subset first with `--n-sample 50`.

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

import pandas as pd
results = pd.read_csv(results_path)
```

### SLURM

```bash
#!/bin/bash
#SBATCH --mem=300G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err

source $(conda info --base)/etc/profile.d/conda.sh
conda activate deseq2_pipeline

python deseq2_pipeline.py \
    --h5ad /path/to/data.h5ad \
    --pert-col target_gene \
    --ctrl-label non-targeting \
    --n-threads 8 \
    --n-workers-r 50
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
| `baseMean` | Mean normalized counts |
| `log2FoldChange` | Log2 fold change (pert vs ctrl) |
| `lfcSE` | Standard error of LFC |
| `stat` | Wald statistic |
| `pvalue` | Nominal p-value |
| `padj` | BH-adjusted p-value |

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
| `--min-cells` | 10 | Perts with fewer cells are skipped and logged to `skipped_perts.csv` |
| `--n-sample` | None | Subsample N random perts. Useful for quick tests before running all perts |
| `--random-seed` | 42 | Controls pseudobulk cell sampling. Fix this to reproduce results exactly |
| `--chunk-size` | 500 | Rows per chunk in sparse row sum. Keep at 500 — lower values avoid RSS spikes on systems with tight memory limits |
