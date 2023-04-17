#!/usr/bin/env python3

import numpy as np
import re
import pandas as pd

path = "/Users/cmdb/gordus-rotation/"
gene_list = pd.read_csv(path + "data/sex_genes_USADLMT_ids.txt", sep='\t')

div_ids = gene_list["Uloborus diversus"].tolist()
mim_ids = gene_list["Stegodyphus mimosarum"].tolist()
bruen_ids = gene_list["Argiope bruennichi"].tolist()

gene_list["U. div chrom"] = None
gene_list["A. breun chrom"] = None
gene_list["D. plan chrom"] = None
gene_list["L. ele chrom"] = None
gene_list["M. bourn chrom"] = None
gene_list["T. clav chrom"] = None

Udiv_chrom3_genes = open(path + "data/Udiv_scaff_3_sex_genes.txt", 'r').read()
Udiv_chrom10_genes = open(path + "data/Udiv_scaff_10_sex_genes.txt", 'r').read()

Abruen_chrom9_genes = open(path + "data/Abruen_chrom_9_sex_genes.txt", 'r').read()
Abruen_chrom10_genes = open(path + "data/Abruen_chrom_10_sex_genes.txt", 'r').read()

Dplan_chrom5_genes = open(path + "data/sex_genes_Dplan_chrom_5.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Dplan_chrom12_genes = open(path + "data/sex_genes_Dplan_chrom_12.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

Lele_chrom1_genes = open(path + "data/sex_genes_Lele_chrom_1.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Lele_chrom9_genes = open(path + "data/sex_genes_Lele_chrom_9.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

Mbourn_chromX1_genes = open(path + "data/sex_genes_Mbourn_chrom_X1.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Mbourn_chromX2_genes = open(path + "data/sex_genes_Mbourn_chrom_X2.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

Tclav_chrom12_genes = open(path + "data/sex_genes_Tclavata_chrom_12.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Tclav_chrom13_genes = open(path + "data/sex_genes_Tclavata_chrom_13.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

def search_genes(id_species, out_species, chromA_file, chromB_file, chromA_name, chromB_name):
	for i in range(len(gene_list)): 
		chromA = False
		chromB = False

		if gene_list[id_species][i] in chromA_file:
			chromA = True
		if gene_list[id_species][i] in chromB_file:
			chromB = True

		if chromA:
			if chromB: 
				gene_list[out_species][i] = "both"
			else: 
				gene_list[out_species][i] = chromA_name

		elif chromB:
			gene_list[out_species][i] = chromB_name

		else: 
			gene_list[out_species][i] = "NEITHER"


diff_chrom_mat = pd.DataFrame(0, columns=div_ids, index=div_ids)

if __name__ == "__main__":
	search_genes(
		"Uloborus diversus", 
		"U. div chrom", 
		Udiv_chrom3_genes, 
		Udiv_chrom10_genes, 
		"chrom_3", 
		"chrom_10"
		)

	search_genes(
		"Argiope bruennichi", 
		"A. breun chrom", 
		Abruen_chrom9_genes, 
		Abruen_chrom10_genes, 
		"chrom_9", 
		"chrom_10"
		)

	search_genes(
		"Stegodyphus mimosarum", 
		"D. plan chrom", 
		Dplan_chrom5_genes, 
		Dplan_chrom12_genes, 
		"chrom_5", 
		"chrom_12"
		)

	search_genes(
		"Stegodyphus mimosarum", 
		"L. ele chrom", 
		Lele_chrom1_genes, 
		Lele_chrom9_genes, 
		"chrom_1", 
		"chrom_9"
		)

	search_genes(
		"Stegodyphus mimosarum", 
		"M. bourn chrom", 
		Mbourn_chromX1_genes, 
		Mbourn_chromX2_genes, 
		"chrom_X1", 
		"chrom_X2"
		)

	search_genes(
		"Stegodyphus mimosarum", 
		"T. clav chrom", 
		Tclav_chrom12_genes, 
		Tclav_chrom13_genes, 
		"chrom_12", 
		"chrom_13"
		)

	gene_list.to_csv(path+'out/gene_loc_db3.csv', index=False)