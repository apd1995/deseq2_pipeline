"""
test_pipeline.py

Unit tests for deseq2_pipeline.py.

Simulates a small dataset (30 perts x 500 genes) and verifies:
  1. sparse_row_sum matches dense sum exactly
  2. sparse_row_sum is chunk-size invariant
  3. Parallel pseudobulk matches sequential pseudobulk exactly
  4. DESeq2 results are numerically identical for seq vs parallel pseudobulk

Run:
    pytest test_pipeline.py -v
"""

import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import os
import shutil
import pytest

from deseq2_pipeline import (
    build_group_structure,
    draw_cell_assignments,
    sparse_row_sum,
    pseudobulk,
    save_pseudobulk,
    write_chunks,
    run_deseq2_batched,
    merge_results,
)

# ── Simulation parameters ────────────────────────────────────
N_PERTS        = 30
N_GENES        = 500
CELLS_PER_PERT = 100
CELLS_CTRL     = 10_000
REP_FRAC       = 0.5
N_REPS         = 2
MIN_CELLS      = 10
RANDOM_SEED    = 42
R_SCRIPT       = "deseq2_worker.R"
CTRL_LABEL     = "control"
PERT_COL       = "perturbation"
TEST_OUTDIR    = "test_results"
N_WORKERS_R    = 5   # small for tests


# ── Fixtures ─────────────────────────────────────────────────

def make_test_adata() -> ad.AnnData:
    """Simulate a small sparse AnnData for testing."""
    np.random.seed(RANDOM_SEED)
    n_cells   = N_PERTS * CELLS_PER_PERT + CELLS_CTRL
    n_nonzero = int(n_cells * N_GENES * 0.1)

    row_idx = np.random.randint(0, n_cells, size=n_nonzero)
    col_idx = np.random.randint(0, N_GENES,  size=n_nonzero)
    values  = np.random.negative_binomial(1, 0.5, size=n_nonzero).astype(np.float32)
    X       = sp.csr_matrix(
        (values, (row_idx, col_idx)), shape=(n_cells, N_GENES)
    )
    pert_labels = (
        [CTRL_LABEL] * CELLS_CTRL
        + [f"pert_{i:02d}" for i in range(1, N_PERTS + 1)
           for _ in range(CELLS_PER_PERT)]
    )
    obs = pd.DataFrame({"perturbation": pert_labels})
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(1, N_GENES + 1)])
    return ad.AnnData(X=X, obs=obs, var=var)


def sort_pseudobulk(
    meta_df: pd.DataFrame, counts_np: np.ndarray
) -> tuple[pd.DataFrame, np.ndarray]:
    """Sort pseudobulk rows by perturbation + replicate + condition."""
    idx = meta_df.sort_values(
        ["perturbation", "replicate", "condition"]
    ).index
    return meta_df.loc[idx].reset_index(drop=True), counts_np[idx]


# ── Tests: sparse_row_sum ────────────────────────────────────

class TestSparseRowSum:

    def test_matches_dense_sum(self):
        """sparse_row_sum should match naive dense sum."""
        np.random.seed(0)
        X_dense  = np.random.randint(0, 10, size=(1000, 200)).astype(np.float32)
        X_sparse = sp.csr_matrix(X_dense)
        idx      = np.random.choice(1000, size=300, replace=False)

        expected = X_dense[idx].sum(axis=0)
        result   = sparse_row_sum(X_sparse, idx, chunk_size=50)

        assert np.allclose(result, expected, rtol=1e-5, atol=1e-5), \
            "sparse_row_sum does not match dense sum"

    def test_chunk_size_invariant(self):
        """Result must be identical regardless of chunk_size."""
        np.random.seed(1)
        X   = sp.csr_matrix(
            np.random.randint(0, 5, size=(500, 100)).astype(np.float32)
        )
        idx = np.arange(200)

        r1 = sparse_row_sum(X, idx, chunk_size=10)
        r2 = sparse_row_sum(X, idx, chunk_size=100)
        r3 = sparse_row_sum(X, idx, chunk_size=500)

        assert np.allclose(r1, r2, rtol=1e-10), \
            "chunk_size=10 vs 100 give different results"
        assert np.allclose(r1, r3, rtol=1e-10), \
            "chunk_size=10 vs 500 give different results"

    def test_empty_idx(self):
        """Empty index should return zero vector."""
        X      = sp.csr_matrix(np.ones((100, 50), dtype=np.float32))
        result = sparse_row_sum(X, np.array([], dtype=int))
        assert result.shape == (50,)
        assert np.all(result == 0)

    def test_single_row(self):
        """Single row index should return that row's values."""
        np.random.seed(2)
        X_dense = np.random.rand(10, 20).astype(np.float32)
        X       = sp.csr_matrix(X_dense)
        idx     = np.array([3])
        result  = sparse_row_sum(X, idx, chunk_size=500)
        assert np.allclose(result, X_dense[3], rtol=1e-5)


# ── Tests: pseudobulk reproducibility ───────────────────────

class TestPseudobulkReproducibility:
    """Parallel pseudobulk must be identical to sequential."""

    def setup_method(self):
        self.adata   = make_test_adata()
        self.run_dir = os.path.join(TEST_OUTDIR, "test_pseudobulk")
        os.makedirs(os.path.join(self.run_dir, "chunks"), exist_ok=True)

        self.group = build_group_structure(
            self.adata, PERT_COL, CTRL_LABEL,
            REP_FRAC, N_REPS, MIN_CELLS, None, RANDOM_SEED,
            self.run_dir,
        )
        self.assignments = draw_cell_assignments(
            self.group, N_REPS, RANDOM_SEED
        )

    def teardown_method(self):
        if os.path.exists(self.run_dir):
            shutil.rmtree(self.run_dir)

    def _run(self, n_threads: int) -> tuple[pd.DataFrame, np.ndarray]:
        return pseudobulk(
            self.adata.X, self.group, self.assignments,
            N_REPS, n_threads=n_threads,
        )

    def _assert_equal(
        self,
        meta1: pd.DataFrame, counts1: np.ndarray,
        meta2: pd.DataFrame, counts2: np.ndarray,
        label: str,
    ):
        meta1_s, counts1_s = sort_pseudobulk(meta1, counts1)
        meta2_s, counts2_s = sort_pseudobulk(meta2, counts2)

        assert meta1_s.equals(meta2_s), \
            f"{label}: metadata mismatch"
        assert np.allclose(counts1_s, counts2_s, rtol=1e-10, atol=1e-10), \
            (
                f"{label}: count matrix mismatch  "
                f"max_diff={np.abs(counts1_s - counts2_s).max():.2e}"
            )

    def test_sequential_vs_2_threads(self):
        seq_meta, seq_counts = self._run(n_threads=1)
        par_meta, par_counts = self._run(n_threads=2)
        self._assert_equal(
            seq_meta, seq_counts, par_meta, par_counts,
            "sequential vs 2 threads",
        )

    def test_sequential_vs_4_threads(self):
        seq_meta, seq_counts = self._run(n_threads=1)
        par_meta, par_counts = self._run(n_threads=4)
        self._assert_equal(
            seq_meta, seq_counts, par_meta, par_counts,
            "sequential vs 4 threads",
        )

    def test_parallel_runs_identical(self):
        """Two parallel runs with same assignments must be identical."""
        meta1, counts1 = self._run(n_threads=4)
        meta2, counts2 = self._run(n_threads=4)
        self._assert_equal(
            meta1, counts1, meta2, counts2,
            "parallel run 1 vs run 2",
        )

    def test_correct_shape(self):
        """Output shape must be (n_perts * n_reps * 2) x n_genes."""
        meta_df, counts_np = self._run(n_threads=1)
        expected_rows = len(self.group["pert_names"]) * N_REPS * 2
        assert counts_np.shape == (expected_rows, N_GENES), \
            f"Expected shape ({expected_rows}, {N_GENES}), got {counts_np.shape}"
        assert len(meta_df) == expected_rows

    def test_all_perts_present(self):
        """All selected perts must appear in the output."""
        meta_df, _ = self._run(n_threads=1)
        output_perts = set(meta_df["perturbation"].unique())
        expected     = set(self.group["pert_names"])
        assert output_perts == expected, \
            f"Missing perts: {expected - output_perts}"

    def test_replicate_balance(self):
        """Each pert must have exactly n_reps pert rows and n_reps ctrl rows."""
        meta_df, _ = self._run(n_threads=1)
        for pert in self.group["pert_names"]:
            subset    = meta_df[meta_df["perturbation"] == pert]
            n_pert    = (subset["condition"] == "pert").sum()
            n_ctrl    = (subset["condition"] == "ctrl").sum()
            assert n_pert == N_REPS, \
                f"{pert}: expected {N_REPS} pert rows, got {n_pert}"
            assert n_ctrl == N_REPS, \
                f"{pert}: expected {N_REPS} ctrl rows, got {n_ctrl}"


# ── Tests: DESeq2 results ────────────────────────────────────

class TestDESeq2Results:
    """DESeq2 results from parallel pseudobulk must match sequential."""

    def setup_method(self):
        self.adata   = make_test_adata()
        self.run_seq = os.path.join(TEST_OUTDIR, "test_deseq2_seq")
        self.run_par = os.path.join(TEST_OUTDIR, "test_deseq2_par")
        os.makedirs(os.path.join(self.run_seq, "chunks"), exist_ok=True)
        os.makedirs(os.path.join(self.run_par, "chunks"), exist_ok=True)

        self.group = build_group_structure(
            self.adata, PERT_COL, CTRL_LABEL,
            REP_FRAC, N_REPS, MIN_CELLS, None, RANDOM_SEED,
            self.run_seq,
        )
        self.assignments = draw_cell_assignments(
            self.group, N_REPS, RANDOM_SEED
        )

    def teardown_method(self):
        for d in [self.run_seq, self.run_par]:
            if os.path.exists(d):
                shutil.rmtree(d)

    def _run_and_collect(
        self, run_dir: str, n_threads: int
    ) -> pd.DataFrame:
        """Run pseudobulk + DESeq2 and return sorted results DataFrame."""
        meta_df, counts_np = pseudobulk(
            self.adata.X, self.group, self.assignments,
            N_REPS, n_threads=n_threads,
        )
        meta_s, counts_s = sort_pseudobulk(meta_df, counts_np)

        _, full_df = save_pseudobulk(
            meta_s, counts_s, self.group["gene_names"], run_dir
        )
        chunk_paths, results_paths = write_chunks(
            full_df, self.group["pert_names"], run_dir
        )

        r_script_abs = os.path.abspath(R_SCRIPT)
        failed = run_deseq2_batched(chunk_paths, r_script_abs, N_WORKERS_R)
        assert not failed, f"DESeq2 workers failed: {failed}"

        results = merge_results(
            results_paths, run_dir, cleanup=True, chunk_paths=chunk_paths
        )
        return (
            results
            .sort_values(["perturbation", "gene"])
            .reset_index(drop=True)
        )

    def test_deseq2_sequential_vs_parallel(self):
        """DESeq2 results must be numerically identical for seq vs parallel."""
        seq_results = self._run_and_collect(self.run_seq, n_threads=1)
        par_results = self._run_and_collect(self.run_par, n_threads=4)

        assert seq_results.shape == par_results.shape, \
            f"Shape mismatch: {seq_results.shape} vs {par_results.shape}"

        for col in ["gene", "perturbation"]:
            assert seq_results[col].equals(par_results[col]), \
                f"Column '{col}' mismatch"

        numeric_cols = [
            "baseMean", "log2FoldChange", "fold_change", "lfcSE", "stat", "pvalue", "fdr"
        ]
        for col in numeric_cols:
            seq_v = seq_results[col].values
            par_v = par_results[col].values

            nan_match = np.all(np.isnan(seq_v) == np.isnan(par_v))
            val_match = np.allclose(
                np.nan_to_num(seq_v),
                np.nan_to_num(par_v),
                rtol=1e-6,
                atol=1e-6,
            )
            assert nan_match and val_match, (
                f"Column '{col}' mismatch  "
                f"max_diff="
                f"{np.abs(np.nan_to_num(seq_v)-np.nan_to_num(par_v)).max():.2e}"
            )

        print(
            f"\nPASS: all {len(numeric_cols)} DESeq2 columns match "
            f"across {len(seq_results['perturbation'].unique())} perts",
            flush=True,
        )

    def test_output_columns(self):
        """Results must contain all expected DESeq2 output columns."""
        results = self._run_and_collect(self.run_seq, n_threads=1)
        expected_cols = {
            "gene", "perturbation", "baseMean",
            "log2FoldChange", "fold_change", "lfcSE", "stat", "pvalue", "fdr",
        }
        assert expected_cols.issubset(set(results.columns)), \
            f"Missing columns: {expected_cols - set(results.columns)}"

    def test_all_perts_in_results(self):
        """All valid perts must appear in the final results."""
        results      = self._run_and_collect(self.run_seq, n_threads=1)
        output_perts = set(results["perturbation"].unique())
        expected     = set(self.group["pert_names"])
        assert output_perts == expected, \
            f"Missing perts in results: {expected - output_perts}"

    def test_results_row_count(self):
        """Results must have n_perts * n_genes rows."""
        results    = self._run_and_collect(self.run_seq, n_threads=1)
        n_perts    = len(self.group["pert_names"])
        expected   = n_perts * N_GENES
        assert len(results) == expected, \
            f"Expected {expected} rows, got {len(results)}"


# ── Cleanup after all tests ───────────────────────────────────

def teardown_module(module):
    if os.path.exists(TEST_OUTDIR):
        shutil.rmtree(TEST_OUTDIR)


# ── Run directly ─────────────────────────────────────────────
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
