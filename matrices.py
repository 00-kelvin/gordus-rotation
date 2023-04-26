#!/usr/bin/env python3

import numpy as np
import pandas as pd

path = "/Users/cmdb/gordus-rotation/"

# import database of gene locations
db_pre = pd.read_csv(path + "out/gene_loc_db4.csv")

# remove genes with multiple copies in any species
db = db_pre[ (db_pre["D. plan chrom"]!="both") & (db_pre["L. ele chrom"]!="both") & (db_pre["M. bourn chrom"]!="both") & (db_pre["T. clav chrom"]!="both")]

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
for id1 in div_ids:
	for id2 in div_ids:
		for species in species_list:

			if db[db["Uloborus diversus"] == id1][species].values[0] == db[db["Uloborus diversus"] == id2][species].values[0]:
				same_chrom_mat.loc[id1,id2] += 1
			else:
				diff_chrom_mat.loc[id1,id2] += 1

# export matrices
diff_chrom_mat.to_csv(path+'out/diff_chrom_mat.csv')
same_chrom_mat.to_csv(path+'out/same_chrom_mat.csv')

