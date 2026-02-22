"""
deseq2_pipeline.py

Pseudobulk + DESeq2 pipeline for single-cell perturbation data.

Usage (CLI):
    python deseq2_pipeline.py \\
        --h5ad /path/to/data.h5ad \\
        --pert-col target_gene \\
        --ctrl-label non-targeting \\
        --outdir results \\
        --n-threads 4 \\
        --n-workers-r 50 \\
        --n-sample 100          # optional: subsample perts for testing

Usage (Python API):
    from deseq2_pipeline import run_pipeline

    results_path = run_pipeline(
        h5ad_path   = "/path/to/data.h5ad",
        pert_col    = "target_gene",
        ctrl_label  = "non-targeting",
        n_threads   = 4,
        n_workers_r = 50,
    )
"""

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import subprocess
import threading
import os
import gc
import json
import time
import argparse
import psutil
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path


# ── Helpers ──────────────────────────────────────────────────

def _ram() -> str:
    proc = psutil.Process(os.getpid())
    return f"process: {proc.memory_info().rss/1e9:.1f}GB"


def _make_run_dir(outdir: str, h5ad_path: str) -> str:
    """Create a unique timestamped run directory to prevent overwrites."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dataset   = Path(h5ad_path).stem[:20]
    run_dir   = os.path.join(outdir, f"run_{timestamp}_{dataset}")
    os.makedirs(os.path.join(run_dir, "chunks"), exist_ok=True)
    return run_dir


def _save_config(run_dir: str, config: dict):
    """Save run parameters to JSON for reproducibility."""
    path = os.path.join(run_dir, "run_config.json")
    with open(path, "w") as f:
        json.dump(config, f, indent=2, default=str)
    print(f"Config saved: {path}", flush=True)


# ── Core functions ───────────────────────────────────────────

def load_adata(h5ad_path: str) -> ad.AnnData:
    """
    Load AnnData from h5ad and ensure CSR sparse format.
    CSR is required for efficient row slicing during pseudobulking.
    """
    print(f"Loading AnnData from {h5ad_path}...", flush=True)
    t0    = time.time()
    adata = ad.read_h5ad(h5ad_path)
    if sp.issparse(adata.X):
        if not isinstance(adata.X, sp.csr_matrix):
            print("Converting to CSR...", flush=True)
            adata.X = adata.X.tocsr()
        print(
            f"X is sparse CSR  "
            f"({adata.X.data.nbytes/1e9:.1f}GB data, "
            f"{adata.X.indices.nbytes/1e9:.1f}GB indices)",
            flush=True,
        )
    else:
        print("X is dense numpy array", flush=True)
    print(f"Loaded: {adata.shape}  {_ram()}  ({time.time()-t0:.1f}s)", flush=True)
    return adata


def build_group_structure(
    adata:       ad.AnnData,
    pert_col:    str,
    ctrl_label:  str,
    rep_frac:    float,
    n_reps:      int,
    min_cells:   int,
    n_sample:    int | None,
    random_seed: int,
    run_dir:     str,
) -> dict:
    """
    Infer perturbation group sizes from obs, filter invalid perturbations,
    and optionally subsample. Returns a group structure dict used downstream.

    Saves skipped_perts.csv to run_dir if any perts are filtered out.
    """
    t0         = time.time()
    pert_col_v = adata.obs[pert_col].values
    gene_names = adata.var_names.tolist()

    unique_perts, counts = np.unique(pert_col_v, return_counts=True)
    value_counts = dict(zip(unique_perts, counts))

    n_ctrl_cells       = int(value_counts[ctrl_label])
    cells_per_ctrl_rep = int(n_ctrl_cells * rep_frac)
    ctrl_idx           = np.where(pert_col_v == ctrl_label)[0]

    pert_valid = []
    skipped    = {}
    for p in unique_perts:
        if p == ctrl_label:
            continue
        n        = int(value_counts[p])
        rep_size = int(n * rep_frac)
        if n < min_cells:
            skipped[p] = f"only {n} cells (< min_cells={min_cells})"
        elif rep_size * n_reps > n:
            skipped[p] = (
                f"not enough cells for {n_reps} non-overlapping reps "
                f"(n={n}, rep_size={rep_size})"
            )
        else:
            pert_valid.append(p)

    if skipped:
        skip_path = os.path.join(run_dir, "skipped_perts.csv")
        pd.DataFrame({
            "perturbation": list(skipped.keys()),
            "reason":       list(skipped.values()),
        }).to_csv(skip_path, index=False)
        print(f"Skipped {len(skipped)} perts → {skip_path}", flush=True)

    np.random.seed(random_seed)
    if n_sample is not None and n_sample < len(pert_valid):
        pert_names = list(
            np.random.choice(pert_valid, size=n_sample, replace=False)
        )
    else:
        pert_names = list(pert_valid)

    pert_rep_sizes = {p: int(value_counts[p] * rep_frac) for p in pert_names}

    print(f"Group structure done ({time.time()-t0:.1f}s):", flush=True)
    print(f"  Valid perts        : {len(pert_valid):,}  (skipped {len(skipped):,})", flush=True)
    print(f"  Running            : {len(pert_names):,}", flush=True)
    print(f"  Genes              : {len(gene_names):,}", flush=True)
    print(f"  Control cells      : {n_ctrl_cells:,}", flush=True)
    print(f"  Ctrl cells/rep     : {cells_per_ctrl_rep:,}  (fixed)", flush=True)
    print(
        f"  Pert rep size range: {min(pert_rep_sizes.values()):,} "
        f"– {max(pert_rep_sizes.values()):,}",
        flush=True,
    )

    return dict(
        pert_col_v         = pert_col_v,
        gene_names         = gene_names,
        ctrl_idx           = ctrl_idx,
        cells_per_ctrl_rep = cells_per_ctrl_rep,
        pert_names         = pert_names,
        pert_rep_sizes     = pert_rep_sizes,
        n_ctrl_cells       = n_ctrl_cells,
    )


def draw_cell_assignments(
    group:       dict,
    n_reps:      int,
    random_seed: int,
) -> dict:
    """
    Pre-draw all cell index assignments for all perturbations and replicates
    using a fixed random seed. This ensures pseudobulk is identical regardless
    of whether it runs sequentially or in parallel.

    Returns:
        dict: pert_name -> list of (p_idx, c_idx) tuples, one per replicate
    """
    np.random.seed(random_seed)

    pert_col_v         = group["pert_col_v"]
    ctrl_idx           = group["ctrl_idx"]
    cells_per_ctrl_rep = group["cells_per_ctrl_rep"]
    pert_rep_sizes     = group["pert_rep_sizes"]
    pert_names         = group["pert_names"]

    assignments = {}
    for pert in pert_names:
        pert_idx    = np.where(pert_col_v == pert)[0]
        ctrl_sample = np.random.choice(
            ctrl_idx, cells_per_ctrl_rep * n_reps, replace=False
        )
        reps = []
        for rep in range(n_reps):
            p_idx = pert_idx[
                rep * pert_rep_sizes[pert] : (rep + 1) * pert_rep_sizes[pert]
            ]
            c_idx = ctrl_sample[
                rep * cells_per_ctrl_rep : (rep + 1) * cells_per_ctrl_rep
            ]
            reps.append((p_idx.copy(), c_idx.copy()))
        assignments[pert] = reps
    return assignments


def sparse_row_sum(
    X:          sp.csr_matrix,
    idx:        np.ndarray,
    chunk_size: int = 500,
) -> np.ndarray:
    """
    Sum rows of a sparse CSR matrix in chunks to avoid large RSS spikes.

    Naive slicing of e.g. 100k rows at once can spike RSS by ~8GB even on
    a sparse matrix with 38% density. chunk_size=500 keeps each spike to
    ~36MB at 18k genes, well within typical memory budgets.

    Args:
        X:          Sparse CSR matrix (cells x genes)
        idx:        Row indices to sum
        chunk_size: Rows per chunk (default 500)

    Returns:
        1D array of shape (n_genes,)
    """
    result = np.zeros(X.shape[1], dtype=np.float64)
    for start in range(0, len(idx), chunk_size):
        chunk_idx = np.sort(idx[start:start + chunk_size])
        chunk     = X[chunk_idx, :]
        result   += np.asarray(chunk.sum(axis=0)).flatten()
        del chunk
    return result


def pseudobulk(
    X:           sp.csr_matrix,
    group:       dict,
    assignments: dict,
    n_reps:      int,
    n_threads:   int,
    chunk_size:  int = 500,
) -> tuple[pd.DataFrame, np.ndarray]:
    """
    Compute pseudobulk counts for all perturbations using pre-drawn
    cell assignments.

    Uses ThreadPoolExecutor for parallelism. Threads share X with zero
    memory copy (unlike multiprocessing). scipy sparse .sum() releases
    the GIL so threads run truly in parallel.

    Args:
        X:           Sparse CSR matrix (cells x genes)
        group:       Output of build_group_structure()
        assignments: Output of draw_cell_assignments()
        n_reps:      Number of replicates
        n_threads:   Number of parallel threads
        chunk_size:  Passed to sparse_row_sum()

    Returns:
        (meta_df, counts_np) where counts_np is (n_samples x n_genes)
    """
    pert_names = group["pert_names"]
    lock       = threading.Lock()
    meta_rows  = []
    count_rows = []
    completed  = [0]
    t_window   = [time.time()]
    t0         = time.time()

    print(
        f"\nPseudobulking {len(pert_names)} perts "
        f"with {n_threads} threads...  {_ram()}",
        flush=True,
    )

    def _pseudobulk_one(pert: str):
        local_meta   = []
        local_counts = []
        for rep, (p_idx, c_idx) in enumerate(assignments[pert]):
            pert_counts = sparse_row_sum(X, p_idx, chunk_size)
            ctrl_counts = sparse_row_sum(X, c_idx, chunk_size)
            local_meta.append({
                "perturbation": pert,
                "replicate":    rep + 1,
                "condition":    "pert",
                "n_pert_cells": len(p_idx),
                "n_ctrl_cells": len(c_idx),
            })
            local_meta.append({
                "perturbation": pert,
                "replicate":    rep + 1,
                "condition":    "ctrl",
                "n_pert_cells": len(p_idx),
                "n_ctrl_cells": len(c_idx),
            })
            local_counts.append(pert_counts)
            local_counts.append(ctrl_counts)

        with lock:
            meta_rows.extend(local_meta)
            count_rows.extend(local_counts)
            completed[0] += 1
            if completed[0] % 10 == 0:
                now = time.time()
                print(
                    f"  [{completed[0]:>4} / {len(pert_names)}]  "
                    f"last 10: {now - t_window[0]:.2f}s, "
                    f"{(now - t_window[0])/10:.2f}s/pert  {_ram()}",
                    flush=True,
                )
                t_window[0] = now

    with ThreadPoolExecutor(max_workers=n_threads) as executor:
        futures = [executor.submit(_pseudobulk_one, pert) for pert in pert_names]
        for f in futures:
            f.result()  # re-raises exceptions from worker threads

    meta_df   = pd.DataFrame(meta_rows)
    counts_np = np.vstack(count_rows)

    # DESeq2 requires integer counts. Pseudobulk sums from float32 sparse
    # matrices can have tiny floating point residuals — round to nearest int.
    if not np.issubdtype(counts_np.dtype, np.integer):
        counts_np = np.round(counts_np).astype(np.int32)

    print(
        f"Pseudobulk done ({time.time()-t0:.1f}s) — "
        f"matrix: {counts_np.shape}  dtype: {counts_np.dtype}  {_ram()}",
        flush=True,
    )
    return meta_df, counts_np


def save_pseudobulk(
    meta_df:    pd.DataFrame,
    counts_np:  np.ndarray,
    gene_names: list,
    run_dir:    str,
) -> tuple[str, pd.DataFrame]:
    """
    Save full pseudobulk matrix (metadata + counts) to run_dir/pseudobulk.csv.gz.

    Returns:
        (path, full_df)
    """
    t0      = time.time()
    full_df = pd.concat(
        [meta_df.reset_index(drop=True),
         pd.DataFrame(counts_np, columns=gene_names)],
        axis=1,
    )
    path = os.path.join(run_dir, "pseudobulk.csv.gz")
    full_df.to_csv(path, index=False)
    print(f"Pseudobulk saved → {path}  ({time.time()-t0:.1f}s)", flush=True)
    return path, full_df


def write_chunks(
    full_df:    pd.DataFrame,
    pert_names: list,
    run_dir:    str,
) -> tuple[list, list]:
    """
    Write one CSV per perturbation into run_dir/chunks/.
    Zero-padded filenames (chunk_0000.csv) to keep directory sorted.

    Returns:
        chunk_paths:   list of (chunk_input, chunk_output, worker_id, pert)
        results_paths: list of chunk_output paths
    """
    t0         = time.time()
    chunks_dir = os.path.join(run_dir, "chunks")
    chunk_paths   = []
    results_paths = []

    for wi, pert in enumerate(pert_names):
        chunk_input  = os.path.join(chunks_dir, f"chunk_{wi:04d}.csv")
        chunk_output = os.path.join(chunks_dir, f"result_{wi:04d}.csv.gz")
        full_df[full_df["perturbation"] == pert].to_csv(chunk_input, index=False)
        chunk_paths.append((chunk_input, chunk_output, str(wi), pert))
        results_paths.append(chunk_output)

    print(
        f"Chunks written → {chunks_dir}  ({time.time()-t0:.1f}s)",
        flush=True,
    )
    return chunk_paths, results_paths


def run_deseq2_batched(
    chunk_paths:  list,
    r_script:     str,
    n_workers_r:  int,
) -> list:
    """
    Launch R workers in batches of n_workers_r. Each worker reads one
    chunk CSV, runs DESeq2, and writes a result CSV. Stdout is streamed
    live to the notebook/terminal via background threads.

    Args:
        chunk_paths:  Output of write_chunks()
        r_script:     Absolute path to deseq2_worker.R
        n_workers_r:  Max parallel R processes per batch

    Returns:
        List of failed worker IDs (empty if all succeeded)
    """
    def _stream(proc, worker_id):
        for line in proc.stdout:
            print(line, end="", flush=True)

    batches = [
        chunk_paths[i:i + n_workers_r]
        for i in range(0, len(chunk_paths), n_workers_r)
    ]
    failed = []
    t0     = time.time()

    print(
        f"\nLaunching R workers: {len(batches)} batches of {n_workers_r}  "
        f"({len(chunk_paths)} perts total)",
        flush=True,
    )

    for bi, batch in enumerate(batches):
        print(
            f"  Batch {bi+1} / {len(batches)}: launching {len(batch)} workers...",
            flush=True,
        )
        procs = []
        for chunk_input, chunk_output, worker_id, pert in batch:
            p = subprocess.Popen(
                ["Rscript", r_script, chunk_input, chunk_output, worker_id],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            procs.append((p, worker_id, pert))

        threads = []
        for proc, worker_id, pert in procs:
            t = threading.Thread(target=_stream, args=(proc, worker_id))
            t.daemon = True
            t.start()
            threads.append(t)
        for t in threads:
            t.join()

        for proc, worker_id, pert in procs:
            proc.wait()
            if proc.returncode != 0:
                err         = proc.stderr.read()
                error_lines = [
                    l for l in err.split("\n")
                    if l.startswith("Error") or l.startswith("Warning")
                ]
                print(
                    f"  Worker {worker_id} ({pert}) FAILED:\n"
                    + "\n".join(error_lines),
                    flush=True,
                )
                failed.append(worker_id)

        print(f"  Batch {bi+1} done", flush=True)

    print(
        f"All R batches finished ({time.time()-t0:.1f}s)  "
        f"{len(chunk_paths)-len(failed)} succeeded, {len(failed)} failed",
        flush=True,
    )
    return failed


def merge_results(
    results_paths: list,
    run_dir:       str,
    cleanup:       bool       = True,
    chunk_paths:   list       = None,
) -> pd.DataFrame:
    """
    Merge per-perturbation result CSVs into run_dir/deseq2_results.csv.gz.
    Optionally clean up chunk files.

    Returns:
        Merged results DataFrame
    """
    t0         = time.time()
    result_dfs = []
    for p in results_paths:
        if os.path.exists(p):
            result_dfs.append(pd.read_csv(p))
        else:
            print(f"WARNING: missing result file {p}", flush=True)

    final      = pd.concat(result_dfs, ignore_index=True)
    final_path = os.path.join(run_dir, "deseq2_results.csv.gz")
    final.to_csv(final_path, index=False)
    print(f"Results merged → {final_path}  ({time.time()-t0:.1f}s)", flush=True)

    if cleanup and chunk_paths:
        for chunk_input, chunk_output, _, _ in chunk_paths:
            for f in [chunk_input, chunk_output]:
                if os.path.exists(f):
                    os.remove(f)
        print("Chunk files cleaned up", flush=True)

    return final


def run_pipeline(
    h5ad_path:   str,
    pert_col:    str,
    ctrl_label:  str,
    outdir:      str        = "results",
    n_threads:   int        = 4,
    n_workers_r: int        = 50,
    rep_frac:    float      = 0.5,
    n_reps:      int        = 2,
    min_cells:   int        = 10,
    n_sample:    int | None = None,
    random_seed: int        = 42,
    chunk_size:  int        = 500,
    r_script:    str        = "deseq2_worker.R",
    cleanup:     bool       = True,
) -> str:
    """
    Full pipeline: load AnnData → pseudobulk → DESeq2 → merge results.

    Each run creates a unique timestamped subdirectory under outdir so
    repeated runs never overwrite each other.

    Args:
        h5ad_path:   Path to input .h5ad file
        pert_col:    obs column containing perturbation labels
        ctrl_label:  Value in pert_col that identifies control cells
        outdir:      Parent output directory (default: "results")
        n_threads:   Threads for pseudobulking (default: 4)
        n_workers_r: Parallel R processes per batch (default: 50)
        rep_frac:    Fraction of cells per replicate (default: 0.5)
        n_reps:      Number of replicates (default: 2)
        min_cells:   Minimum cells for a pert to be included (default: 10)
        n_sample:    Subsample N random perts; None = run all (default: None)
        random_seed: Random seed for reproducibility (default: 42)
        chunk_size:  Rows per chunk in sparse_row_sum (default: 500)
        r_script:    Path to deseq2_worker.R (default: "deseq2_worker.R")
        cleanup:     Delete chunk files after merging (default: True)

    Returns:
        Path to deseq2_results.csv.gz
    """
    # prevent R/BLAS from spawning extra threads and competing with Python
    os.environ["OMP_NUM_THREADS"]      = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"]      = "1"

    t_start = time.time()
    run_dir = _make_run_dir(outdir, h5ad_path)
    print(f"Run directory: {run_dir}", flush=True)

    _save_config(run_dir, dict(
        h5ad_path   = h5ad_path,
        pert_col    = pert_col,
        ctrl_label  = ctrl_label,
        n_threads   = n_threads,
        n_workers_r = n_workers_r,
        rep_frac    = rep_frac,
        n_reps      = n_reps,
        min_cells   = min_cells,
        n_sample    = n_sample,
        random_seed = random_seed,
        chunk_size  = chunk_size,
        r_script    = r_script,
    ))

    r_script_abs = os.path.abspath(r_script)
    assert os.path.exists(r_script_abs), \
        f"R script not found: {r_script_abs}\n" \
        f"Make sure deseq2_worker.R is in the same directory or pass --r-script."

    # load
    adata = load_adata(h5ad_path)

    # group structure
    group = build_group_structure(
        adata, pert_col, ctrl_label, rep_frac, n_reps,
        min_cells, n_sample, random_seed, run_dir,
    )

    # pre-draw cell assignments with fixed seed — ensures reproducibility
    assignments = draw_cell_assignments(group, n_reps, random_seed)

    # pseudobulk
    meta_df, counts_np = pseudobulk(
        adata.X, group, assignments, n_reps, n_threads, chunk_size
    )

    # free adata — no longer needed after pseudobulking
    del adata
    gc.collect()
    print(f"adata freed  {_ram()}", flush=True)

    # save
    _, full_df = save_pseudobulk(
        meta_df, counts_np, group["gene_names"], run_dir
    )
    chunk_paths, results_paths = write_chunks(
        full_df, group["pert_names"], run_dir
    )
    del full_df, counts_np, meta_df
    gc.collect()

    # DESeq2
    failed = run_deseq2_batched(chunk_paths, r_script_abs, n_workers_r)
    if failed:
        print(f"WARNING: {len(failed)} workers failed: {failed}", flush=True)

    # merge
    merge_results(results_paths, run_dir, cleanup, chunk_paths)

    t_total = time.time() - t_start
    print(f"\n{'='*55}", flush=True)
    print(f"  TOTAL   : {t_total:.1f}s  ({t_total/60:.1f} min)", flush=True)
    print(f"  Output  : {run_dir}", flush=True)
    print(f"{'='*55}", flush=True)

    return os.path.join(run_dir, "deseq2_results.csv.gz")


# ── CLI entry point ──────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Pseudobulk + DESeq2 pipeline for scRNA-seq perturbation data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--h5ad",        required=True,
                        help="Path to input .h5ad file")
    parser.add_argument("--pert-col",    required=True,
                        help="obs column with perturbation labels")
    parser.add_argument("--ctrl-label",  required=True,
                        help="Value in pert-col that identifies control cells")
    parser.add_argument("--outdir",      default="results",
                        help="Parent output directory")
    parser.add_argument("--n-threads",   type=int,   default=4,
                        help="Threads for pseudobulking")
    parser.add_argument("--n-workers-r", type=int,   default=50,
                        help="Parallel R processes per batch")
    parser.add_argument("--rep-frac",    type=float, default=0.5,
                        help="Fraction of cells per replicate")
    parser.add_argument("--n-reps",      type=int,   default=2,
                        help="Number of replicates")
    parser.add_argument("--min-cells",   type=int,   default=10,
                        help="Min cells per pert to be included")
    parser.add_argument("--n-sample",    type=int,   default=None,
                        help="Subsample N perts (default: run all)")
    parser.add_argument("--random-seed", type=int,   default=42,
                        help="Random seed")
    parser.add_argument("--chunk-size",  type=int,   default=500,
                        help="Rows per chunk in sparse row sum")
    parser.add_argument("--r-script",    default="deseq2_worker.R",
                        help="Path to deseq2_worker.R")
    parser.add_argument("--no-cleanup",  action="store_true",
                        help="Keep chunk files after run (useful for debugging)")

    args = parser.parse_args()

    results_path = run_pipeline(
        h5ad_path   = args.h5ad,
        pert_col    = args.pert_col,
        ctrl_label  = args.ctrl_label,
        outdir      = args.outdir,
        n_threads   = args.n_threads,
        n_workers_r = args.n_workers_r,
        rep_frac    = args.rep_frac,
        n_reps      = args.n_reps,
        min_cells   = args.min_cells,
        n_sample    = args.n_sample,
        random_seed = args.random_seed,
        chunk_size  = args.chunk_size,
        r_script    = args.r_script,
        cleanup     = not args.no_cleanup,
    )
    print(f"\nResults: {results_path}")
