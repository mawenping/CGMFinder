# CGMFinder

Correlated gene module identification from scRNA-seq data.

## Requirements

- Python >= 3.12

## Installation

```bash
# Step 1: create conda environment
conda create -n cgmfinder python=3.12
conda activate cgmfinder

# Step 2: install locked dependencies (exact versions for reproducibility)
pip install -r requirements_locked.txt

# Step 3: install CGMFinder
pip install -e .
```

If you do not need exact reproducibility, simply run `pip install -e .` after creating the environment.
The results will be biologically equivalent but may differ slightly due to BLAS library differences.

### Quick check

```bash
cgmfinder --help
```

## BLAS dependency

The Pearson correlation matrix is computed via `numpy.corrcoef()`, which relies on the BLAS library
(Basic Linear Algebra Subprograms) for matrix operations. Different BLAS implementations
(e.g., OpenBLAS vs MKL) or different BLAS versions can produce results that differ at the
machine-precision level (≈ 1 × 10⁻¹⁵). While individually negligible, these tiny differences can
accumulate and affect the final edge-binning results.

This is a known property of numerical computing and is not a bug in CGMFinder. To reproduce
bit-identical results, use the exact same environment (see `requirements_locked.txt`).

## Usage

```bash
cgmfinder input_name.feather output_dir/ --cpu_cores 5
```

### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--top_neighbor_fraction` | `0.005` | Fraction of top correlated genes retained per gene for network construction. If >1, treated as absolute neighbour count. |
| `--core_size` | `4` | Minimum number of genes required for a module to be recognised after consecutive network refinement. |
| `--min_edges_per_bin` | `300` | Minimum number of edges accumulated in each bin during sequential edge binning. Smaller values improve accuracy but increase computation time. |
| `--dominance_factor` | `2.0` | Ratio requiring the strongest module connection to be at least this many times the second-strongest before assignment. |
| `--triangle_confirmation_fraction` | `0.3` | Minimum fraction of connections a gene must maintain to members of a candidate module in the triangle network to confirm assignment. |
| `--min_cell_fraction` | `0.1` | Minimum fraction of cells in which a gene must be detected to be retained. Controls the number of genes entering network construction. |
| `--cpu_cores` | `1` | Number of CPU cores (recommended ≤ 5). More cores do not improve performance. |
| `--ercc_correction` | `True` | Perform ERCC spike-in correction when ERCC genes are found. Set to `False` to skip (e.g. for toy data without ERCC). |

### Input format

The expression matrix can be provided in three formats:

| Format | Extension | Notes |
|--------|-----------|-------|
| Feather | `.feather` | Fastest, recommended. Columns = genes, rows = cells. |
| CSV | `.csv` | First column is treated as row names (`index_col=0`). Columns = genes, rows = cells. |
| AnnData | `.h5ad` | Directly loaded via `scanpy.read_h5ad()`. Orientation is checked automatically. |

All formats require: **rows = cells, columns = genes**.

## Output files

All output files are written to the specified output directory. File names follow the pattern:

```
input_name_0.005_4_300_2.0_0.3_0.1.run.pars
input_name_0.005_4_300_2.0_0.3_0.1.0.1.h5ad
input_name_0.005_4_300_2.0_0.3_0.1.pearson.feather
input_name_0.005_4_300_2.0_0.3_0.1.module.pickle
input_name_0.005_4_300_2.0_0.3_0.1.geneattr.pickle
input_name_0.005_4_300_2.0_0.3_0.1.dismtx.feather
```

The filename prefix encodes the run parameters:
`{basename}_{top_neighbor_fraction}_{core_size}_{min_edges_per_bin}_{dominance_factor}_{triangle_confirmation_fraction}_{min_cell_fraction}`.

For example, `input_name_0.005_4_300_2.0_0.3_0.1` means
`top_neighbor_fraction=0.005`, `core_size=4`, `min_edges_per_bin=300`, `dominance_factor=2.0`,
`triangle_confirmation_fraction=0.3`, `min_cell_fraction=0.1`.

| File | Format | Description |
|------|--------|-------------|
| `*.run.pars` | plain text | Run parameters recorded for reproducibility |
| `*.h5ad` | AnnData (h5ad) | Filtered, normalised, log-transformed expression matrix (checkpoint) |
| `*.pearson.feather` | Feather (DataFrame) | Pairwise Pearson correlation matrix for genes passing `min_cell_fraction` filter (checkpoint) |
| `*.module.pickle` | pickle (dict) | **Identified modules.** Dict of `{module_name: [gene_list]}` |
| `*.geneattr.pickle` | pickle (dict) | Per-gene attributes including module assignment, expression statistics, and network properties |
| `*.dismtx.feather` | Feather (DataFrame) | Module distance matrix for downstream module-level analysis |

> **Note:** `.h5ad` and `.pearson.feather` are checkpoint files. If they already exist from a previous run,
> they are not recomputed. Delete them to force a fresh computation.

## Examples

### Run with Jupyter notebook

Open `analysis/run_cgmfinder.ipynb` in Jupyter Lab and run the cells sequentially.
The notebook uses the toy dataset (500 cells × 3000 genes) in `analysis/data/`.

### Slurm cluster submission

An example script is provided at `analysis/run_cgmfinder.slurm`:

```bash
sbatch analysis/run_cgmfinder.slurm
```

The script runs CGMFinder on the toy dataset with default parameters.
Adjust `--cpu_cores` and `--mem` as needed for your data.

## Downstream analysis

### Module clustering

The `analysis/module_clustering.py` script demonstrates hierarchical clustering on the module distance matrix (`*.dismtx.feather`):

```bash
python analysis/module_clustering.py
```
## License

MIT
