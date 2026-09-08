"""CGMFinder: Correlated Gene Module Identification from scRNA-seq Data.

This module identifies correlated gene modules (CGMs) from scRNA-seq
expression data using a graph-based filtering approach. It first identifies
CGM cores from highly-expressed genes, then links lowly-expressed genes
to these cores. A distance matrix between modules is also computed for
module-level functional analysis.

Typical usage:

    >>> import cgmfinder as cgmf
    >>> cgmf.run_analysis("expression_matrix.feather", "output_dir/")
"""

from __future__ import annotations

import argparse
import copy
import importlib.resources
import inspect
import itertools
import logging
import os
import pickle
import time
from collections import defaultdict
from multiprocessing import Pool
from typing import Dict, List, Optional, Tuple

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

logger = logging.getLogger(__name__)

# Global state shared across the CGM calling workflow.
# These are set by run_analysis() and read/written by downstream functions.
# Note: the module is not re-entrant - only one analysis at a time.
modules: Dict[str, list] = {}
module_dis: Dict[int, pd.DataFrame] = {}
gene_attr: Dict[str, dict] = {}
mr_rank: pd.DataFrame = None  # type: ignore[assignment]

# ---------------------------------------------------------------------------
# ERCC correction
# ---------------------------------------------------------------------------


def _load_cell_cycle_genes() -> Optional[pd.DataFrame]:
    """Load the cell-cycle gene reference shipped with the package."""
    try:
        package = importlib.resources.files("cgmfinder")
        resource = package / "data" / "cell_cycle_genes.csv"
        with resource.open("r", encoding="utf-8") as fh:
            return pd.read_csv(fh, header=0)
    except (FileNotFoundError, ModuleNotFoundError):
        logger.warning("Cell-cycle gene file not found; skipping ERCC correction.")
        return None


def ercc_correction(
    edata: ad.AnnData, ercc_genes: List[str]
) -> ad.AnnData:
    """Correct ERCC spike-in expression bias across cell-cycle phases.

    Reads the bundled cell-cycle gene list, scores cell-cycle phases via
    Scanpy, and scales ERCC counts in G1 cells to match the mean expression
    level of S/G2M cells.
    """
    ccg = _load_cell_cycle_genes()
    if ccg is None:
        return edata

    s_genes = ccg[ccg.phase == "S"]["ensembl_id"].tolist()
    g2m_genes = ccg[ccg.phase == "G2M"]["ensembl_id"].tolist()
    sc.tl.score_genes_cell_cycle(edata, s_genes=s_genes, g2m_genes=g2m_genes)

    if ercc_genes:
        g1_mean = np.mean(
            edata[edata.obs["phase"] == "G1", ercc_genes].X.sum(axis=1)
        ).item()
        sg2_mean = np.mean(
            edata[edata.obs["phase"] != "G1", ercc_genes].X.sum(axis=1)
        ).item()
        if g1_mean > 0 and sg2_mean > 0:
            # Scale G1 ERCC counts directly on the underlying .X array.
            # Do not assign to the .X of an AnnData view: anndata >= 0.13
            # "initializes view as actual" on view .X writes, which changed
            # the resulting edge-binning (anndata 0.12 -> 0.13 comparison).
            g1_idx = np.flatnonzero((edata.obs["phase"] == "G1").values)
            ercc_idx = edata.var_names.get_indexer(ercc_genes)
            scale = sg2_mean / g1_mean
            ix = np.ix_(g1_idx, ercc_idx)
            if hasattr(edata.X, "toarray"):
                X = edata.X.tocsr().tolil()
                X[ix] = X[ix] * scale
                edata.X = X.tocsr()
            else:
                edata.X[ix] = edata.X[ix] * scale
    return edata


# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------


def pre_processing(exp_mtx: str, min_cell_fraction: float, flag: bool = True) -> ad.AnnData:
    """Load, filter, and normalise an expression matrix.

    Parameters
    ----------
    exp_mtx:
        Path to a ``.h5ad``, ``.feather``, or ``.csv`` file.
        Rows = cells, columns = genes.  For ``.csv`` / ``.feather`` the matrix is
        loaded into a new AnnData object; ``.h5ad`` is read directly.
    min_cell_fraction:
        Minimum fraction of cells in which a gene must be detected to be retained.
        Controls the number of genes entering network construction.
    flag:
        If *True*, attempt ERCC spike-in correction when ERCC genes are found.

    Returns
    -------
    AnnData with normalised expression, log-transformed, and cell-cycle scores.
    """
    logger.info("Preprocessing started at %s", time.ctime())

    if exp_mtx.endswith(".feather"):
        umi_matrix = pd.read_feather(exp_mtx)
    elif exp_mtx.endswith(".csv"):
        umi_matrix = pd.read_csv(exp_mtx, index_col=0)
    elif exp_mtx.endswith(".h5ad"):
        edata = ad.read_h5ad(exp_mtx)
        # Check orientation: rows must be cells, columns must be genes
        if edata.shape[0] > edata.shape[1]:
            logger.warning(
                f"Input .h5ad has {edata.shape[0]} rows and {edata.shape[1]} columns. "
                f"Expected format: rows = cells, columns = genes."
            )
        # AnnData already built; skip re-initialisation below
    else:
        raise ValueError(
            "Unsupported file format: " + exp_mtx + ". "
            "Only .h5ad, .feather, and .csv are allowed."
        )

    if not exp_mtx.endswith(".h5ad"):
        edata = ad.AnnData(umi_matrix)
        # Check orientation: rows must be cells, columns must be genes
        if edata.shape[0] > edata.shape[1]:
            logger.warning(
                f"Input {os.path.splitext(exp_mtx)[1]} has {edata.shape[0]} rows "
                f"and {edata.shape[1]} columns. "
                f"Expected format: rows = cells, columns = genes."
            )
    sc.pp.filter_genes(edata, min_cells=1)
    edata.obs["umi_counts"] = edata.X.sum(axis=1)
    edata.var["mean_exp_unnor"] = edata.X.mean(axis=0)
    sc.pp.normalize_total(edata, target_sum=1e4)
    edata.var["mean_exp_nor"] = edata.X.mean(axis=0)
    sc.pp.log1p(edata)
    edata.var["mean_exp_log"] = edata.X.mean(axis=0)

    ercc_genes = [g for g in edata.var_names if str(g).startswith("ERCC-")]
    logger.info("ERCC gene count: %d, flag: %s", len(ercc_genes), flag)
    if ercc_genes and flag:
        edata = ercc_correction(edata, ercc_genes)
        logger.info("ERCC correction finished")
    else:
        logger.info("ERCC correction not performed")

    logger.info("Preprocessing finished at %s", time.ctime())
    logger.info("*" * 70)
    return edata
def compute_corr_mr(
    edata: ad.AnnData,
    min_cell_fraction: float,
    out_dir: str,
    f_base_name: str,
) -> pd.DataFrame:
    """Compute pairwise gene Pearson correlation and mutual rank matrix.

    The correlation matrix is saved as a Feather file, and the geometric
    mutual-rank matrix (geometric mean of bidirectional ranks) is returned.
    """
    logger.info("Compute correlation and mutual ranks started at %s", time.ctime())

    min_cells = int(edata.shape[0] * min_cell_fraction)
    logger.info("Genes expressed in less than %d cells are removed", min_cells)

    logger.info("Use Pearson correlation")
    sc.pp.filter_genes(edata, min_cells=min_cells)
    corr = np.corrcoef(edata.to_df().values, rowvar=False)
    target_gene_list = list(edata.var_names)

    edata_name = os.path.join(out_dir, f"{f_base_name}_{min_cell_fraction}.h5ad")
    if not os.path.exists(edata_name):
        edata.write(edata_name)

    corr_name = os.path.join(
        out_dir, f"{f_base_name}_{min_cell_fraction}.pearson.feather"
    )
    if not os.path.exists(corr_name):
        pd.DataFrame(corr, index=target_gene_list, columns=target_gene_list).to_feather(
            corr_name
        )

    np.fill_diagonal(corr, -1.0)
    df_corr = pd.DataFrame(corr, index=target_gene_list, columns=target_gene_list)
    df_rank = df_corr.rank(ascending=False, method="min")

    logger.info("Use geometric mean for mutual ranking calculation")
    mr_result = np.sqrt(df_rank * df_rank.T)

    logger.info("Compute correlation and mutual ranks finished at %s", time.ctime())
    logger.info("%d genes remained for CGM identification", corr.shape[0])
    logger.info("*" * 70)
    return mr_result

# ---------------------------------------------------------------------------
# Network binning
# ---------------------------------------------------------------------------


def get_network_bins(
    mr_rank: pd.DataFrame, top_neighbor_fraction: float, min_edges_per_bin: int
) -> Tuple[List[List[Tuple]], int]:
    """Bin edges by mutual-rank value.

    For each gene, the closest neighbours (rank-based) are
    collected, deduplicated, sorted by mutual rank, and split into bins of
    at least *min_edges_per_bin* edges each.
    """
    logger.info("Edge binning started at %s", time.ctime())
    if top_neighbor_fraction < 1:
        top_neighbor_fraction = mr_rank.shape[0] * top_neighbor_fraction
    top_neighbor_fraction = int(top_neighbor_fraction)

    # Collect unique edges (soft threshold: include all ties at boundary)
    edge_map: Dict[frozenset, float] = {}
    for col in mr_rank.columns:
        sorted_vals = mr_rank[col].sort_values()
        cutoff_idx = min(top_neighbor_fraction - 1, len(sorted_vals) - 1)
        cutoff_val = sorted_vals.iloc[cutoff_idx]
        for gene, mr_val in sorted_vals[sorted_vals <= cutoff_val].items():
            key = frozenset((col, gene))
            edge_map[key] = mr_val

    edges_df = pd.DataFrame(
        [(a, b, c) for (a, b), c in edge_map.items()],
        columns=["n1", "n2", "mr"],
    ).sort_values(by="mr")

    # Bin
    mr_bins: List[List[Tuple]] = []
    current_bin: List[Tuple] = []
    for mr_val in sorted(edges_df["mr"].unique()):
        chunk = list(edges_df[edges_df["mr"] == mr_val].itertuples(index=False))
        current_bin.extend(chunk)
        if len(current_bin) >= min_edges_per_bin:
            mr_bins.append(current_bin)
            current_bin = []
    if current_bin:
        mr_bins.append(current_bin)

    logger.info("Top neighbor fraction: %s", str(top_neighbor_fraction))
    logger.info("Total bins: %d", len(mr_bins))
    logger.info("Min edges per bin: %d", min_edges_per_bin)
    logger.info("Edge binning finished at %s", time.ctime())
    logger.info("*" * 70)
    return mr_bins, top_neighbor_fraction

# ---------------------------------------------------------------------------
# Graph helpers
# ---------------------------------------------------------------------------


def _safe_cpu_count() -> int:
    """Return the number of usable CPUs, defaulting to 1."""
    return max(1, (os.cpu_count() or 1) - 1)


def add_weighted_edges_undirected(
    adj: Dict[str, Dict[str, float]], edges: List[Tuple[str, str, float]]
) -> Dict[str, Dict[str, float]]:
    for n1, n2, w in edges:
        adj.setdefault(n1, {})[n2] = w
        adj.setdefault(n2, {})[n1] = w
    return adj


def non_triangle_edges(
    adj: Dict[str, Dict[str, float]]
) -> List[Tuple[str, str]]:
    edges: List[Tuple[str, str]] = [
        (u, v) for u in adj for v in adj[u]
    ]
    return [
        (u, v) for u, v in edges
        if not (adj[u].keys() & adj[v].keys())
    ]


def remove_isolated_nodes(adj: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    for node in list(adj.keys()):
        if not adj[node]:
            del adj[node]
    return adj


def remove_edges(
    adj: Dict[str, Dict[str, float]], edges_to_remove: List[Tuple[str, str]]
) -> Dict[str, Dict[str, float]]:
    for u, v in edges_to_remove:
        adj.get(u, {}).pop(v, None)
        adj.get(v, {}).pop(u, None)
    return adj


def remove_nodes(
    adj: Dict[str, Dict[str, float]], nodes_to_remove: List[str]
) -> Dict[str, Dict[str, float]]:
    adj = remove_edges(
        adj, [(u, v) for u in nodes_to_remove for v in adj.get(u, {})]
    )
    for node in nodes_to_remove:
        adj.pop(node, None)
    return adj


def remove_non_triangle_edges(
    raw_adj: Dict[str, Dict[str, float]], cp: bool = True
) -> Dict[str, Dict[str, float]]:
    tri_adj = copy.deepcopy(raw_adj) if cp else raw_adj
    non_tri = non_triangle_edges(tri_adj)
    tri_adj = remove_edges(tri_adj, non_tri)
    tri_adj = remove_isolated_nodes(tri_adj)
    return tri_adj


def get_connected_components(
    adj: Dict[str, Dict[str, float]]
) -> List[List[str]]:
    nodes = set(adj.keys())
    components: List[List[str]] = []
    while nodes:
        start = next(iter(nodes))
        component = {start}
        frontier = {start}
        while frontier:
            frontier = set(itertools.chain.from_iterable(
                adj[n].keys() for n in frontier
            )) - component
            component |= frontier
        components.append(sorted(component))
        nodes -= component
    return components


def non_consecutive_nodes(
    adj: Dict[str, Dict[str, float]], target_node: str
) -> Optional[str]:
    neighbors = list(adj.get(target_node, {}).keys())
    if len(neighbors) < 2:
        return None
    sub_adj = {
        n: {i: adj[n][i] for i in adj[n] if i in neighbors}
        for n in neighbors
    }
    components = get_connected_components(sub_adj)
    return target_node if len(components) > 1 else None


def remove_non_consecutive_nodes(
    tri_adj: Dict[str, Dict[str, float]], p: int, cp: bool = True
) -> Dict[str, Dict[str, float]]:
    con_adj = copy.deepcopy(tri_adj) if cp else tri_adj
    module_genes = set(itertools.chain.from_iterable(modules.values())) if modules else set()
    candidates = {n for n in con_adj if len(con_adj[n]) > 3}
    any_removed = True
    while any_removed:
        q = min(_safe_cpu_count(), p) if len(candidates) > 100 else 1
        with Pool(q) as pool:
            non_con = [
                n for n in pool.starmap(
                    non_consecutive_nodes,
                    [(con_adj, node) for node in candidates],
                )
                if n is not None
            ]
        if not non_con:
            break
        candidates = set(itertools.chain.from_iterable(
            con_adj[n].keys() for n in non_con
        ))
        bmds = set(non_con) & module_genes
        nodes_to_remove = set(non_con) - bmds
        edges_to_remove = []
        if bmds:
            for bmd in bmds:
                bmd_module = gene_attr.get(bmd, {}).get("module")
                if bmd_module:
                    edges_to_remove.extend(
                        (bmd, n) for n in con_adj[bmd]
                        if n not in modules.get(bmd_module, [])
                    )
        con_adj = remove_edges(con_adj, edges_to_remove)
        con_adj = remove_nodes(con_adj, list(nodes_to_remove))
        con_adj = remove_non_triangle_edges(con_adj, cp=False)
        candidates = {n for n in (candidates & con_adj.keys()) if len(con_adj[n]) > 3}
        any_removed = bool(nodes_to_remove) or bool(edges_to_remove)
    return con_adj

# ---------------------------------------------------------------------------
# Module identification & gene assignment
# ---------------------------------------------------------------------------


def identify_modules(
    con_adj: Dict[str, Dict[str, float]],
    tri_adj: Dict[str, Dict[str, float]],
    core_size: int,
    bin_order: int,
    weight: float,
    top_neighbor_fraction: int,
) -> None:
    module_genes = set(itertools.chain.from_iterable(modules.values())) if modules else set()
    new_cores = sorted(
        c for c in get_connected_components(con_adj)
        if len(c) >= core_size and not (set(c) & module_genes)
    )
    for order, members in enumerate(new_cores):
        name = f"{bin_order}_{order}"
        modules[name] = members
        for gene in members:
            set_gene_attr(gene, tri_adj, bin_order, weight, top_neighbor_fraction, core=True)


def belonged_module(geneid: str) -> Optional[str]:
    return next(
        (md for md, genes in modules.items() if geneid in genes), None
    )


def neighbored_module(
    gene_list: List[str], t: int
) -> List[Tuple[str, set]]:
    nmd = {}
    for md, genes in modules.items():
        overlap = set(gene_list) & set(genes)
        if len(overlap) > t:
            nmd[md] = overlap
    return sorted(nmd.items(), key=lambda x: len(x[1]), reverse=True)


def raw_neighbors(
    gene: str, weight: float, top_neighbor_fraction: int, t: int
) -> Tuple[List[str], List[Tuple[str, set]]]:
    temp = mr_rank[gene].sort_values()
    gnbs = temp[temp <= weight].index.tolist()
    if len(gnbs) > top_neighbor_fraction:
        gnbs = gnbs[:top_neighbor_fraction]
    nmd = neighbored_module(gnbs, t)
    return gnbs, nmd


def tri_neighbors(
    gene: str, gnbs: List[str], tri_adj: Dict, t: int
) -> Tuple[set, List[Tuple[str, set]]]:
    tnbs = set(tri_adj.get(gene, {}).keys()) & set(gnbs)
    nmd = neighbored_module(list(tnbs), t)
    return tnbs, nmd


def set_gene_attr(
    gene: str,
    graph: Dict[str, Dict[str, float]],
    bin_order: int,
    weight: float,
    top_neighbor_fraction: int,
    core: bool = False,
) -> None:
    if gene not in gene_attr or gene_attr[gene].get("module"):
        return
    gnbs, gnb_nmd = raw_neighbors(gene, weight, top_neighbor_fraction, 0)
    tnbs, tnb_nmd = tri_neighbors(gene, gnbs, graph, 0)
    bmd = belonged_module(gene)
    gene_attr[gene]["gnb"] = {k[0]: len(k[1]) for k in gnb_nmd} if gnb_nmd else None
    gene_attr[gene]["tnb"] = {k[0]: len(k[1]) for k in tnb_nmd} if tnb_nmd else None
    gene_attr[gene]["bin_order"] = bin_order
    gene_attr[gene]["mr_value"] = weight
    gene_attr[gene]["top_neighbor"] = len(gnbs)
    gene_attr[gene]["core"] = core
    if bmd:
        gene_attr[gene]["md_size"] = len(modules[bmd])
        gene_attr[gene]["module"] = bmd

# ---------------------------------------------------------------------------
# Candidate filtering & assignment
# ---------------------------------------------------------------------------


def raw_md_candidate(
    gnb_nmd: List[Tuple[str, set]], dominance_factor: float
) -> Optional[List[str]]:
    """Return [module_name] if the top raw-connected module dominates."""
    if not gnb_nmd:
        return None
    if len(gnb_nmd) == 1:
        return [gnb_nmd[0][0]]
    if len(gnb_nmd[0][1]) >= dominance_factor * len(gnb_nmd[1][1]):
        return [gnb_nmd[0][0]]
    return None


def tri_md_candidate(
    tnb_nmd: List[Tuple[str, set]],
    top_neighbor_fraction: int,
    triangle_confirmation_fraction: float,
    s: float = 0.5,
) -> List[str]:
    """Return modules whose triangle-connected neighbours meet threshold."""
    candidates = []
    if not tnb_nmd:
        return candidates
    for md, members in tnb_nmd:
        threshold = max(
            3,
            min(
                round(s * top_neighbor_fraction),
                round(triangle_confirmation_fraction * len(modules.get(md, []))),
            ),
        )
        if len(members) >= threshold:
            candidates.append(md)
    return candidates


def if_connected(
    adj: Dict[str, Dict[str, float]], neighbors: set
) -> bool:
    sub_adj = {
        n: {i: adj[n][i] for i in adj[n] if i in neighbors}
        for n in neighbors
    }
    components = get_connected_components(sub_adj)
    return max(len(c) for c in components) > 1


def gene_assignment(
    gene: str,
    tri_adj: Dict[str, Dict[str, float]],
    weight: float,
    dominance_factor: float,
    triangle_confirmation_fraction: float,
    top_neighbor_fraction: int,
) -> Tuple[Tuple[str, str], List[str]]:
    tri_neighbor_genes: List[str] = []
    gnbs, gnb_nmd = raw_neighbors(gene, weight, top_neighbor_fraction, 1)
    tnbs, tnb_nmd = tri_neighbors(gene, gnbs, tri_adj, 1)

    raw_cand = raw_md_candidate(gnb_nmd, dominance_factor)
    tri_cand = tri_md_candidate(tnb_nmd, top_neighbor_fraction, triangle_confirmation_fraction)

    if raw_cand and tri_cand:
        common = list(set(raw_cand) & set(tri_cand))
        if common:
            candidate_md = common[0]
            neighbors = set(tri_adj.get(gene, {}).keys()) & set(modules.get(candidate_md, []))
            if if_connected(tri_adj, neighbors):
                return (candidate_md, gene), list(tnbs)
    return ("not_assigned", gene), tri_neighbor_genes


def gene_assignment_mp(
    tri_adj: Dict[str, Dict[str, float]],
    bin_order: int,
    weight: float,
    top_neighbor_fraction: int,
    dominance_factor: float,
    triangle_confirmation_fraction: float,
    p: int,
) -> None:
    module_genes = set(itertools.chain.from_iterable(modules.values())) if modules else set()
    candidates = set(tri_adj.keys()) - module_genes
    if not candidates:
        return

    any_assigned = True
    while any_assigned:
        any_assigned = False
        q = min(_safe_cpu_count(), p) if len(candidates) > 100 else 1
        with Pool(q) as pool:
            results = pool.starmap(
                gene_assignment,
                [(c, tri_adj, weight, dominance_factor, triangle_confirmation_fraction, top_neighbor_fraction)
                 for c in candidates],
            )
        new_candidates: set = set()
        for (assigned_md, assigned_gene), neighbor_genes in results:
            if assigned_md == "not_assigned":
                set_gene_attr(assigned_gene, tri_adj, bin_order, weight,
                              top_neighbor_fraction, core=False)
            else:
                modules.setdefault(assigned_md, []).append(assigned_gene)
                set_gene_attr(assigned_gene, tri_adj, bin_order, weight,
                              top_neighbor_fraction, core=False)
                new_candidates.update(neighbor_genes)
                any_assigned = True
        module_genes = set(itertools.chain.from_iterable(modules.values())) if modules else set()
        candidates = new_candidates - module_genes
# ---------------------------------------------------------------------------
# Module connectivity
# ---------------------------------------------------------------------------


def module_connectivity(
    adj: Dict[str, Dict[str, float]], bin_order: int
) -> None:
    md_names = list(modules.keys())
    df = pd.DataFrame(index=md_names, columns=md_names, dtype=float)
    for md in md_names:
        neighbors = {
            nb for gene in modules[md] for nb in adj.get(gene, {})
        } - set(modules[md])
        row = {}
        for mdd in md_names:
            overlap = len(neighbors & set(modules[mdd]))
            row[mdd] = (
                overlap / (len(modules[md]) + len(modules[mdd]))
                if overlap > 1 else 0.0
            )
        df[md] = pd.Series(row)
    symmetric = df * df.T
    symmetric[symmetric > 0] = 1
    module_dis[bin_order] = (df * symmetric) + (df * symmetric).T


# ---------------------------------------------------------------------------
# Module distance matrix (AUC-based)
# ---------------------------------------------------------------------------


def _compute_auc_for_module(
    module_name: str, module_dis_dict: Dict[int, pd.DataFrame]
) -> pd.DataFrame:
    """Stack connectivity values for one module across all layers."""
    layers = sorted(module_dis_dict)
    return pd.concat(
        [module_dis_dict[l][module_name].rename(l) for l in layers],
        axis=1,
    )


def compute_auc_dis_mtx(
    module_dis_dict: Dict[int, pd.DataFrame]
) -> pd.DataFrame:
    """Compute the module distance matrix via AUC-based scoring.

    All layers are aligned to the same index/columns (the union from the
    last layer). For each module pair the accumulated connectivity across
    bins is normalised by the number of bins in which they co-appear.
    """
    layers = sorted(module_dis_dict)
    last = module_dis_dict[layers[-1]]

    # Fill missing rows/cols in each layer
    for layer in layers:
        cur = module_dis_dict[layer]
        missing = sorted(set(last.index) - set(cur.index))
        if missing:
            filler = pd.DataFrame(0.0, index=missing, columns=missing)
            cur = pd.concat([cur, filler]).infer_objects().fillna(0.0)
        module_dis_dict[layer] = cur.reindex(index=last.index, columns=last.columns)

    auc_list = []
    background_list = []
    for md in last.index:
        df_auc = _compute_auc_for_module(md, module_dis_dict)
        background = pd.Series(
            {m: max((df_auc.loc[m] > 0).sum(), 1) for m in df_auc.index},
            name=md,
        )
        background_list.append(background)
        auc_list.append(pd.Series(df_auc.sum(axis=1), name=md))

    df_sum = pd.concat(auc_list, axis=1)
    df_bg = pd.concat(background_list, axis=1)
    return 1.0 - df_sum / df_bg


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def write_results(bin_order: str, out_dir: str, f_name: str) -> None:
    md_path = os.path.join(out_dir, f"{f_name}.module.pickle")
    attr_path = os.path.join(out_dir, f"{f_name}.geneattr.pickle")
    with open(md_path, "wb") as fh:
        pickle.dump(modules, fh)
    with open(attr_path, "wb") as fh:
        pickle.dump(gene_attr, fh)


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def run_analysis(
    exp_mtx: str,
    out_dir: str,
    top_neighbor_fraction: float = 0.005,
    core_size: int = 4,
    min_edges_per_bin: int = 300,
    dominance_factor: float = 2.0,
    triangle_confirmation_fraction: float = 0.3,
    min_cell_fraction: float = 0.1,
    cpu_cores: int = 1,
    ercc_correction: bool = True,
) -> None:
    """Run the full CGMFinder pipeline.

    Parameters
    ----------
    exp_mtx:
        Path to the expression matrix (``.feather`` or ``.csv``).
    out_dir:
        Directory for all output files (created if it does not exist).
    top_neighbor_fraction:
        Fraction of top correlated genes retained per gene for the initial network.
        If >1, interpreted as the absolute number of top neighbours.
    core_size:
        Minimum number of genes required for a module to be recognised after
        consecutive network refinement.
    min_edges_per_bin:
        Minimum number of edges accumulated in each bin during sequential
        edge binning.
    dominance_factor:
        Ratio requiring the strongest module connection to be at least this
        many times the second-strongest before assignment.
    triangle_confirmation_fraction:
        Minimum fraction of connections a gene must maintain to members of a
        candidate module in the triangle network to confirm assignment.
    min_cell_fraction:
        Minimum fraction of cells in which a gene must be detected to be
        retained.  Directly controls gene inclusion.
    cpu_cores:
        Number of CPU cores for parallel processing.
    ercc_correction:
        If *True*, perform ERCC spike-in correction when ERCC genes are found.
        Set to *False* to skip (e.g. for toy data without ERCC spike-ins).
        If set, filter modules with fewer than this many genes from the
        final distance matrix (not yet implemented).
    """
    global mr_rank, modules, module_dis, gene_attr

    ts = time.time()
    logger.info("Run started at %s", time.ctime())
    logger.info("CGMFinder call module started at %s", time.ctime())

    exp_mtx = os.path.abspath(exp_mtx)
    out_dir = os.path.abspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    f_base_name = os.path.basename(exp_mtx).rsplit(".", 1)[0]

    # Build a descriptive file name from key parameters
    nt_str = str(int(top_neighbor_fraction)) if top_neighbor_fraction >= 1 else str(top_neighbor_fraction)
    f_name = "_".join([
        f_base_name,
        nt_str,
        str(core_size),
        str(min_edges_per_bin),
        str(dominance_factor),
        str(triangle_confirmation_fraction),
        str(min_cell_fraction),
    ])

    # Write run parameters
    frame = inspect.currentframe()
    args, _, _, values = inspect.getargvalues(frame)
    pars_path = os.path.join(out_dir, f"{f_name}.run.pars")
    with open(pars_path, "w") as fh:
        fh.write("Run parameters:\n")
        for arg in args:
            fh.write(f"{arg}: {values[arg]}\n")
            logger.info("    %s: %s", arg, values[arg])
    logger.info("-" * 70)

    # Reset global state
    modules.clear()
    module_dis.clear()
    gene_attr.clear()
    mr_rank = None  # type: ignore[assignment]

    edata = pre_processing(exp_mtx, min_cell_fraction, flag=ercc_correction)
    mr_rank = compute_corr_mr(edata, min_cell_fraction,
                              out_dir, f_base_name)

    mr_bins, top_neighbor_fraction = get_network_bins(mr_rank, top_neighbor_fraction, min_edges_per_bin)

    # Initialise gene_attr from normalised expression data
    for gene in mr_rank.columns:
        gene_attr[gene] = {
            "n_cell": int(edata.var.loc[gene, "n_cells"]),
            "mean_exp_unnor": float(edata.var.loc[gene, "mean_exp_unnor"]),
            "mean_exp_nor": float(edata.var.loc[gene, "mean_exp_nor"]),
            "mean_exp_log": float(edata.var.loc[gene, "mean_exp_log"]),
            "module": None,
            "md_size": None,
            "gnb": None,
            "tnb": None,
            "bin_order": None,
            "mr_value": None,
            "top_neighbor": None,
            "core": False,
        }

    logger.info("CGM calling started at %s", time.ctime())
    raw_adj: Dict[str, Dict[str, float]] = defaultdict(dict)

    for bin_order, ebin in enumerate(mr_bins):
        t1 = time.time()
        weight = ebin[-1][-1]
        raw_adj = add_weighted_edges_undirected(raw_adj, ebin)
        tri_adj = remove_non_triangle_edges(raw_adj, cp=True)
        con_adj = remove_non_consecutive_nodes(tri_adj, cpu_cores, cp=True)
        identify_modules(con_adj, tri_adj, core_size, bin_order, weight, top_neighbor_fraction)
        gene_assignment_mp(tri_adj, bin_order, weight, top_neighbor_fraction,
                           dominance_factor, triangle_confirmation_fraction, cpu_cores)
        module_connectivity(tri_adj, bin_order)
        elapsed = time.time() - t1
        total_assigned = len(list(itertools.chain.from_iterable(modules.values()))) if modules else 0
        logger.info(
            "Bin:%d finished, raw_nodes:%d, tri_nodes:%d, con_nodes:%d, "
            "module:%d, assigned_genes:%d, time:%.2f seconds",
            bin_order, len(raw_adj), len(tri_adj), len(con_adj),
            len(modules), total_assigned, elapsed,
        )

    write_results("final", out_dir, f_name)
    logger.info("CGM calling finished at %s", time.ctime())
    logger.info("*" * 70)

    # Compute module distance matrix
    logger.info("CGM connection analysis started at %s", time.ctime())
    distance_matrix = compute_auc_dis_mtx(module_dis)
    dis_name = os.path.join(out_dir, f"{f_name}.dismtx.feather")
    distance_matrix.to_feather(dis_name)
    logger.info("CGM connection analysis finished at %s", time.ctime())

    te = time.time()
    logger.info("Run finished at %s", time.ctime())
    logger.info("Total runtime: %.2f seconds", te - ts)
    logger.info("*" * 70)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    parser = argparse.ArgumentParser(
        description="Correlated gene module identification from scRNA-seq data"
    )
    parser.add_argument("exp_mtx", type=str, help="Path to the expression matrix file")
    parser.add_argument("out_dir", type=str, help="Output directory path")
    parser.add_argument(
        "--top_neighbor_fraction", type=float, default=0.005,
        help=(
            "Fraction of top correlated genes retained per gene for network "
            "construction; if >1 treated as absolute neighbour count "
            "(default: 0.005)"
        ),
    )
    parser.add_argument(
        "--core_size", type=int, default=4,
        help="Minimum genes for a module to be recognised (default: 4)",
    )
    parser.add_argument(
        "--min_edges_per_bin", type=int, default=300,
        help="Minimum edges per bin during sequential edge binning (default: 300)",
    )
    parser.add_argument(
        "--dominance_factor", type=float, default=2.0,
        help="Strength ratio for dominant module connection (default: 2.0)",
    )
    parser.add_argument(
        "--triangle_confirmation_fraction", type=float, default=0.3,
        help="Min fraction of triangle connections to confirm assignment "
        "(default: 0.3)",
    )
    parser.add_argument(
        "--min_cell_fraction", type=float, default=0.1,
        help=(
            "Minimum fraction of cells a gene must be detected in "
            "(default: 0.1)"
        ),
    )
    parser.add_argument(
        "--cpu_cores", type=int, default=1,
        help="Number of CPU cores to use (default: 1)",
    )
    parser.add_argument(
        "--ercc_correction", type=bool, default=True,
        help="Perform ERCC spike-in correction when ERCC genes are found (default: True)",
    )

    args = parser.parse_args()
    run_analysis(
        args.exp_mtx,
        args.out_dir,
        top_neighbor_fraction=args.top_neighbor_fraction,
        core_size=args.core_size,
        min_edges_per_bin=args.min_edges_per_bin,
        dominance_factor=args.dominance_factor,
        triangle_confirmation_fraction=args.triangle_confirmation_fraction,
        min_cell_fraction=args.min_cell_fraction,
        cpu_cores=args.cpu_cores,
        ercc_correction=args.ercc_correction,
    )


if __name__ == "__main__":
    main()
