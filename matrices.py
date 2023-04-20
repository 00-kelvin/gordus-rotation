#!/usr/bin/env python3

import numpy as np
import pandas as pd

path = "/Users/cmdb/gordus-rotation/"

# import database of gene locations
db = pd.read_csv(path + "out/gene_loc_db4.csv")

# extract U. div gene IDs
div_ids = db["Uloborus diversus"].tolist()

# initialize matrices
diff_chrom_mat = pd.DataFrame(0, columns=div_ids, index=div_ids)
same_chrom_mat = pd.DataFrame(0, columns=div_ids, index=div_ids)

# list of species columns in database
species_list = ["U. div chrom", 
				"A. breun chrom", 
				"D. plan chrom", 
				"L. ele chrom", 
				"M. bourn chrom",
				"T. clav chrom"]

# check each pair of genes
for id1 in range(len(db)):
	for id2 in range(len(db)):
		for species in species_list:

			# for now, ignore pairs in which either gene is on both chroms 
			if db[species][id1] == "both":
				break
			elif db[species][id2] == "both":
				break
			elif db[species][id1] == db[species][id2]:
				same_chrom_mat.loc[div_ids[id1],div_ids[id2]] = same_chrom_mat[div_ids[id1]][div_ids[id2]] + 1
			else:
				diff_chrom_mat.loc[div_ids[id1],div_ids[id2]] = diff_chrom_mat[div_ids[id1]][div_ids[id2]] + 1

# export matrices
diff_chrom_mat.to_csv(path+'out/diff_chrom_mat.csv')
same_chrom_mat.to_csv(path+'out/same_chrom_mat.csv')