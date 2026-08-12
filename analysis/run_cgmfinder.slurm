#!/bin/bash
#SBATCH -p compute_pro
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 6
#SBATCH --time=7-00:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=cgmfinder
#SBATCH --output=cgmfinder_%j.out

# ====== Environment ======
source /path/to/miniconda3/etc/profile.d/conda.sh
conda activate cgmfinder

# ====== Paths ======
# Toy dataset: 500 cells x 3000 genes
umi=analysis/data/toy.feather
out_dir=data/output

# ====== Run CGMFinder ======
cgmfinder $umi $out_dir \
    --top_neighbor_fraction 0.005 \
    --core_size 4 \
    --min_edges_per_bin 300 \
    --dominance_factor 2.0 \
    --triangle_confirmation_fraction 0.3 \
    --min_cell_fraction 0.1 \
    --cpu_cores 5 \
    --ercc_correction False
