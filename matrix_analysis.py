#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns

path = "/Users/cmdb/gordus-rotation/"

# import database of gene locations
loc_db = pd.read_csv(path + "out/gene_loc_db4.csv")

# import segregation matrices
diff_mat = pd.read_csv(path + "out/diff_chrom_mat.csv", index_col=0)
same_mat = pd.read_csv(path + "out/same_chrom_mat.csv", index_col=0)


# sub db for genes where all 4 added species had hits on both chroms
both_chrom = loc_db[ (loc_db["D. plan chrom"]=="both") & (loc_db["L. ele chrom"]=="both") & (loc_db["M. bourn chrom"]=="both") & (loc_db["T. clav chrom"]=="both")]

# make new spreadsheet
# both_chrom.to_csv(path+'out/both_chrom.csv')

print("Count of genes found on both chroms for all 4 added species: " + str(len(both_chrom)))

sixes = diff_mat[diff_mat == 6].count()

pair_count = int(sixes.sum() / 2)
genes_in_pairs = np.count_nonzero(sixes.values)

print("Fully segregating pairs: " + str(pair_count))
print("Count of genes involved in fully segregating pairs: " + str(genes_in_pairs))

# Make heatmap of the differently segregating gene pairs

# fig, ax = plt.subplots(figsize=(8,7))
# sns.heatmap(diff_mat, ax=ax, square=True, cmap="magma_r")
# plt.suptitle("Number of species for which the gene pair are found on different sex chromosomes")
# ax.set_xlabel("Gene 2")
# ax.set_ylabel("Gene 1")

# plt.tight_layout()
# plt.savefig(path + "out/diff_chrom_heatmap.png", dpi=200)

# Make heatmap of the same-segregating gene pairs

fig, ax = plt.subplots(figsize=(8,7))
sns.heatmap(same_mat, ax=ax, square=True, cmap="magma_r")
plt.suptitle("Number of species for which the gene pair are found on the same sex chromosomes")
ax.set_xlabel("Gene 2")
ax.set_ylabel("Gene 1")

plt.tight_layout()
plt.savefig(path + "out/same_chrom_heatmap.png", dpi=200)

# Make histogram of counts of 6s

# fig, ax = plt.subplots()
# sns.histplot(sixes, bins = 30, binrange = (0,30))
# ax.set_xlabel("sixes")
# plt.savefig(path + "out/sixes_hist.png", dpi=200)
