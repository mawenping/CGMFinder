"""CGMFinder - correlated gene module identification from scRNA-seq data."""

from cgmfinder.cgmfinder import (
    run_analysis,
    pre_processing,
    ercc_correction,
    compute_corr_mr,
    main,
)

__all__ = [
    "run_analysis",
    "pre_processing",
    "ercc_correction",
    "compute_corr_mr",
    "main",
]
