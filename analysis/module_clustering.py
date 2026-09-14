"""
Module clustering analysis using the module distance matrix (*.dismtx.feather).

This script reads the module distance matrix produced by CGMFinder and performs
hierarchical clustering to identify module clades/clusters.
"""

import os
import pickle
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

# ============================================================
# Configuration - adjust these paths
# ============================================================
md_dir = "/path/to/output_dir/"          # CGMFinder output directory
f_name = "input_name_0.005_4_300_2.0_0.3_0.1"  # output filename prefix

# ============================================================
# 1. Load module distance matrix
# ============================================================
dis = pd.read_feather(os.path.join(md_dir, f"{f_name}.dismtx.feather"))
print(f"Module distance matrix: {dis.shape[0]} modules")

# ============================================================
# 2. Hierarchical clustering of modules
# ============================================================
# Extract upper-triangular distance values
dis_np = dis.to_numpy()
dis_1d = dis_np[np.triu_indices(dis_np.shape[0], k=1)]

# Single-linkage clustering
Z = scipy.cluster.hierarchy.linkage(dis_1d, method="single")

# ============================================================
# 3. Plot dendrogram (horizontal)
# ============================================================
sns.set(context="paper", style="ticks")
fig, ax = plt.subplots(figsize=(3, 15))

dendro = scipy.cluster.hierarchy.dendrogram(
    Z,
    leaf_font_size=12,
    labels=dis.columns,
    color_threshold=0.3,
    distance_sort="ascending",
    count_sort="ascending",
    orientation="left",
    above_threshold_color="#000000",
    ax=ax,
)
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8])
ax.tick_params(axis="x", labelsize=12)
plt.tight_layout()
plt.savefig("module_dendrogram.pdf", format="pdf", dpi=300, bbox_inches="tight")
plt.show()

# ============================================================
# 4. Get module order from clustering
# ============================================================
module_order = dendro["ivl"][::-1]  # top-to-bottom order
print(f"Module order: {module_order}")

# ============================================================
# 5. Identify clades / cut the dendrogram
# ============================================================
# Clusters are determined by color_threshold=0.3 in dendrogram
from scipy.cluster.hierarchy import fcluster
clusters = fcluster(Z, t=0.3, criterion="distance")
cluster_df = pd.DataFrame({"module": dis.columns, "cluster": clusters})
print(cluster_df)

# ============================================================
# 6. (Optional) Label modules with clade numbers
# ============================================================
leaves_color = dendro["leaves_color_list"]
print("Module -> cluster colour mapping:")
for md, col in zip(dendro["ivl"], leaves_color):
    print(f"  {md}: {col}")
