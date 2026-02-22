#!/bin/bash
# setup.sh
#
# Fully automated setup for deseq2_pipeline on HPC SLURM clusters.
# No sudo required. Installs R 4.5.2 + all dependencies via conda + pak.
#
# Usage:
#   bash setup.sh
#
# What it does:
#   1. Finds conda/mamba
#   2. Creates or updates the conda environment
#   3. Installs Python 3.11 + R 4.5.2 + compiler tools via conda
#   4. Installs R packages (DESeq2, glmGamPoi, data.table) via pak inside R
#   5. Verifies everything works
#   6. Prints a ready-to-use SLURM job template

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

info()    { echo -e "${GREEN}[INFO]${NC}  $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }

echo ""
echo "=================================================="
echo "  deseq2_pipeline HPC setup"
echo "=================================================="
echo ""

# ── Step 1: Find conda/mamba ─────────────────────────────────
info "Checking for conda/mamba..."

CONDA_CMD=""
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
    success "Found mamba: $(mamba --version | head -1)"
elif command -v conda &>/dev/null; then
    CONDA_CMD="conda"
    success "Found conda: $(conda --version)"
else
    for candidate in \
        "$HOME/miniconda3/bin/conda" \
        "$HOME/anaconda3/bin/conda" \
        "/opt/conda/bin/conda" \
        "/usr/local/anaconda/bin/conda" \
        "/apps/conda/bin/conda"
    do
        if [ -f "$candidate" ]; then
            source "$(dirname "$candidate")/../etc/profile.d/conda.sh"
            CONDA_CMD="conda"
            success "Found conda at $candidate"
            break
        fi
    done
fi

if [ -z "$CONDA_CMD" ]; then
    error "conda not found. Install Miniconda first:
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p \$HOME/miniconda3
    source \$HOME/miniconda3/etc/profile.d/conda.sh
    conda init bash
Then re-run this script."
fi

# ── Step 2: Create or update conda environment ───────────────
ENV_NAME="deseq2_pipeline"
info "Setting up conda environment: $ENV_NAME..."

get_prefix() {
    conda env list | grep "^${ENV_NAME} " | awk '{print $NF}'
}

if conda env list | grep -q "^${ENV_NAME} "; then
    warn "Environment '$ENV_NAME' already exists — checking for missing packages..."
else
    info "Creating new environment '$ENV_NAME'..."
    $CONDA_CMD create -n "$ENV_NAME" \
        -c conda-forge \
        python=3.11 \
        pip \
        -y
    success "Base environment created"
fi

CONDA_PREFIX=$(get_prefix)
[ -z "$CONDA_PREFIX" ] && error "Could not find conda prefix for $ENV_NAME"

PYTHON_BIN="$CONDA_PREFIX/bin/python"
RSCRIPT_BIN="$CONDA_PREFIX/bin/Rscript"

# ── Step 3: Clean conda cache ────────────────────────────────
info "Cleaning conda package cache..."
$CONDA_CMD clean --packages --tarballs -y
success "Conda cache cleaned"

# ── Step 4: Install all conda packages ───────────────────────
info "Installing Python + R 4.5.2 + compiler tools via conda..."
info "This may take 10-20 min on first run..."

$CONDA_CMD install -n "$ENV_NAME" \
    -c conda-forge \
    --no-update-deps \
    -y \
    python=3.11 \
    pip \
    numpy \
    pandas \
    scipy \
    anndata \
    h5py \
    psutil \
    pytest \
    r-base=4.5.2 \
    compilers \
    make \
    libcurl \
    openssl \
    zlib

success "Python + R 4.5.2 + compiler tools installed"

# ── Step 5: Verify R is present ──────────────────────────────
if [ ! -f "$RSCRIPT_BIN" ]; then
    error "Rscript not found at $RSCRIPT_BIN after installation"
fi

info "Environment prefix : $CONDA_PREFIX"
info "Python             : $($PYTHON_BIN --version)"
info "R                  : $($RSCRIPT_BIN --version 2>&1 | head -1)"

# ── Step 6: Install R packages via pak ───────────────────────
# conda-forge does not have Bioconductor packages for R 4.5.2 yet
# pak installs from source with full dependency resolution
# minimal list — only what deseq2_worker.R actually needs
info "Installing R packages via pak (10-20 min first run)..."

# explicitly set compiler paths so pak subprocess can find gcc/g++
export PATH="$CONDA_PREFIX/bin:$PATH"
export CC="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gcc"
export CXX="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-g++"
export FC="$CONDA_PREFIX/bin/x86_64-conda-linux-gnu-gfortran"

$RSCRIPT_BIN - <<'REOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

# bootstrap pak using stable binary repo
cat("==> Bootstrapping pak...\n")
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak",
                   repos = "https://r-lib.github.io/p/pak/stable/",
                   quiet = FALSE)
}
cat(sprintf("  pak version: %s\n", as.character(packageVersion("pak"))))

# install minimal set — pak resolves all true transitive dependencies
# (DESeq2 requires SummarizedExperiment, HDF5Array etc — all handled automatically)
cat("==> Installing DESeq2, glmGamPoi, data.table, R.utils...\n")
pak::pkg_install(
  c(
    "data.table",
    "R.utils",
    "bioc::DESeq2",
    "bioc::glmGamPoi"
  ),
  ask     = FALSE,
  upgrade = FALSE
)

# verify all four packages load
cat("\n==> Verifying packages load correctly...\n")
pkgs   <- c("DESeq2", "glmGamPoi", "data.table", "R.utils")
failed <- c()
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
    cat(sprintf("  FAIL  %s\n", pkg))
  } else {
    cat(sprintf("  OK    %-15s %s\n", pkg, as.character(packageVersion(pkg))))
  }
}
if (length(failed) > 0) {
  cat(sprintf("\nFailed: %s\n", paste(failed, collapse = ", ")))
  quit(status = 1)
}
cat("\nAll R packages installed and verified.\n")
REOF

success "R packages installed and verified"

# ── Step 7: Verify Python imports ────────────────────────────
info "Verifying Python imports..."
$PYTHON_BIN - <<'PYEOF'
import importlib, sys
pkgs   = ["numpy", "pandas", "scipy", "anndata", "h5py", "psutil"]
failed = []
for pkg in pkgs:
    try:
        m   = importlib.import_module(pkg)
        ver = getattr(m, "__version__", "?")
        print(f"  OK    {pkg:<15} {ver}")
    except ImportError as e:
        print(f"  FAIL  {pkg}: {e}")
        failed.append(pkg)
if failed:
    print(f"\nFailed: {failed}")
    sys.exit(1)
print("\nAll Python packages verified.")
PYEOF

success "Python imports verified"

# ── Step 8: Verify Rscript callable from Python ──────────────
info "Verifying Rscript callable from Python subprocess..."
$PYTHON_BIN - <<PYEOF
import subprocess, sys
result = subprocess.run(
    ["$RSCRIPT_BIN", "--version"],
    capture_output=True, text=True
)
if result.returncode != 0:
    print("FAIL: Rscript not callable via subprocess")
    sys.exit(1)
ver = (result.stdout or result.stderr).strip().splitlines()[0]
print(f"  OK    {ver}")
PYEOF

success "Rscript callable from Python"

# ── Done ─────────────────────────────────────────────────────
echo ""
echo "=================================================="
echo -e "${GREEN}  Setup complete!${NC}"
echo "=================================================="
echo ""
echo "Activate the environment:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Run tests:"
echo "  conda activate $ENV_NAME && pytest test_pipeline.py -v"
echo ""
echo "Run the pipeline:"
echo "  conda activate $ENV_NAME"
echo "  python deseq2_pipeline.py \\"
echo "      --h5ad /path/to/data.h5ad \\"
echo "      --pert-col target_gene \\"
echo "      --ctrl-label non-targeting"
echo ""

# ── SLURM template ────────────────────────────────────────────
echo "=================================================="
echo "  SLURM job template (save as submit_job.sh)"
echo "=================================================="
cat << SLURMEOF

#!/bin/bash
#SBATCH --job-name=deseq2_pipeline
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20
#SBATCH --time=04:00:00
#SBATCH --output=logs/%j.out
#SBATCH --error=logs/%j.err

mkdir -p logs
source \$(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

python run_analysis.py

SLURMEOF
echo "Run with: sbatch submit_job.sh"
echo ""
