#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import scipy.cluster.hierarchy as sch

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


# compute clusters

# retrieve clusters using fcluster 
d = sch.distance.pdist(diff_mat)
L = sch.linkage(d, method='average')
# 0.2 can be modified to retrieve more stringent or relaxed clusters
clusters = sch.fcluster(L, 0.1*d.max(), 'distance')

# clusters indicices correspond to incides of original df
cluster_list = [[diff_mat.index[i], cluster] for i,cluster in enumerate(clusters)]

clusters_df = pd.DataFrame(cluster_list, columns = ['gene', 'cluster'])
# print(clusters_df)

sort the clusters 
sorted_df = clusters_df.sort_values(by=['cluster'])
print(sorted_df)

# Locating the relevant clusters using genes indicated on the heatmap

print("Group A, block 1") # cluster 13
print(clusters_df[clusters_df['gene']=='mRNA16344'])
# print(clusters_df[clusters_df['gene']=='mRNA19239'])
# print(clusters_df[clusters_df['gene']=='mRNA16056'])

print("Group B, block 1") # cluster 16
print(clusters_df[clusters_df['gene']=='mRNA7981'])
# print(clusters_df[clusters_df['gene']=='mRNA7245'])
# print(clusters_df[clusters_df['gene']=='mRNA6555'])
# print(clusters_df[clusters_df['gene']=='mRNA6391'])

print("Group A, block 2") # cluster 12
print(clusters_df[clusters_df['gene']=='mRNA7367'])
# print(clusters_df[clusters_df['gene']=='mRNA7046'])

print("Group B, block 2") # cluster 17
print(clusters_df[clusters_df['gene']=='mRNA20205'])
# print(clusters_df[clusters_df['gene']=='mRNA19863'])
# print(clusters_df[clusters_df['gene']=='mRNA17756'])

print("Group A, block 3") # cluster 9
print(clusters_df[clusters_df['gene']=='mRNA6646'])

print("Group B, block 3") # cluster 18
print(clusters_df[clusters_df['gene']=='mRNA20021'])
# print(clusters_df[clusters_df['gene']=='mRNA18991'])
# print(clusters_df[clusters_df['gene']=='mRNA18800'])

blocks_df = clusters_df.loc[clusters_df['cluster'].isin([13,16,12,17,9,18])]

# make a column indicating which block of genes we're looking at
blocks_df['block'] = 0
blocks_df.loc[blocks_df['cluster'].isin([13,16]), 'block'] = 1
blocks_df.loc[blocks_df['cluster'].isin([12,17]), 'block'] = 2
blocks_df.loc[blocks_df['cluster'].isin([9,18]), 'block'] = 3

blocks_df['group'] = 'C'
blocks_df.loc[blocks_df['cluster'].isin([13,12,9]), 'group'] = 'A'
blocks_df.loc[blocks_df['cluster'].isin([16,17,18]), 'group'] = 'B'

# print(blocks_df)
sorted_df = blocks_df.sort_values(by=['block', 'group', 'gene'])
sorted_df.to_csv(path+'out/blocks.csv')


# Make clustermap of the differently segregating gene pairs

fig, ax = plt.subplots()

sns.clustermap(diff_mat, row_linkage = L, col_linkage = L)
plt.suptitle("Number of species for which the gene pair are found on different sex chromosomes")
ax.set_xlabel("Gene 2")
ax.set_ylabel("Gene 1")
plt.tight_layout()
plt.savefig(path + "out/diff_chrom_clustermap3.png", dpi=200)

# Make clustermap of the differently segregating gene pairs (OLD WAY)

fig, ax = plt.subplots()
sns.clustermap(diff_mat)
plt.suptitle("Number of species for which the gene pair are found on different sex chromosomes")
ax.set_xlabel("Gene 2")
ax.set_ylabel("Gene 1")

plt.tight_layout()
plt.savefig(path + "out/diff_chrom_clustermap2.png", dpi=200)


# Make heatmap of the same-segregating gene pairs

fig, ax = plt.subplots(figsize=(8,7))
sns.heatmap(same_mat, ax=ax, square=True)
plt.suptitle("Number of species for which the gene pair are found on the same sex chromosomes")
ax.set_xlabel("Gene 2")
ax.set_ylabel("Gene 1")

plt.tight_layout()
plt.savefig(path + "out/same_chrom_heatmap.png", dpi=200)

# Make histogram of counts of 6s

fig, ax = plt.subplots()
sns.histplot(sixes, bins = 30, binrange = (0,30))
ax.set_xlabel("sixes")
plt.savefig(path + "out/sixes_hist.png", dpi=200)

# Finding genes always on a different chrom from mRNA6646
print(diff_mat[diff_mat["mRNA6646"]==6]["mRNA6646"])
