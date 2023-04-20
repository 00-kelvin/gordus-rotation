#!/usr/bin/env python3

import numpy as np
import pandas as pd

path = "/Users/cmdb/gordus-rotation/"

# import the list of genes found in all 7 species
gene_list = pd.read_csv(path + "data/sex_genes_USADLMT_ids.txt", sep='\t')

# get gene names from different species
div_ids = gene_list["Uloborus diversus"].tolist()
mim_ids = gene_list["Stegodyphus mimosarum"].tolist()
bruen_ids = gene_list["Argiope bruennichi"].tolist()

# make columns for which chromosome that gene is on in each of the 6 species
# (does not include S. mimosarum, since we do not know which chromosomes the sex
# genes came from)
gene_list["U. div chrom"] = None
gene_list["A. breun chrom"] = None
gene_list["D. plan chrom"] = None
gene_list["L. ele chrom"] = None
gene_list["M. bourn chrom"] = None
gene_list["T. clav chrom"] = None

# import files listing the genes on each chromosome for each species

## U. diversus
Udiv_chrom3_genes = open(path + "data/Udiv.scaffold_3.cds.fa", 'r').read()
Udiv_chrom10_genes = open(path + "data/Udiv.scaffold_10.cds.fa", 'r').read()

## A. bruennechi
Abruen_chrom9_genes = open(path + "data/Abruen.chrom_9.gff", 'r').read()
Abruen_chrom10_genes = open(path + "data/Abruen.chrom_10.gff", 'r').read()

## D. plantarius
Dplan_chrom5_genes = open(path + "data/sex_genes_Dplan_chrom_5.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Dplan_chrom12_genes = open(path + "data/sex_genes_Dplan_chrom_12.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

## L. elegans
Lele_chrom1_genes = open(path + "data/sex_genes_Lele_chrom_1.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Lele_chrom9_genes = open(path + "data/sex_genes_Lele_chrom_9.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

## M. bourneti
Mbourn_chromX1_genes = open(path + "data/sex_genes_Mbourn_chrom_X1.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Mbourn_chromX2_genes = open(path + "data/sex_genes_Mbourn_chrom_X2.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

## T. clavata
Tclav_chrom12_genes = open(path + "data/sex_genes_Tclavata_chrom_12.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()
Tclav_chrom13_genes = open(path + "data/sex_genes_Tclavata_chrom_13.blastn.outfmt6_max-hits-1_max-hsps-1", 'r').read()

# create function to find which chromosome each gene falls on for a given species
def search_genes(id_species, out_species, chromA_file, chromB_file, chromA_name, chromB_name):
	
	# for each gene in the list
	for i in range(len(gene_list)): 
		
		# initialize bools for whether the gene is on each chromosome
		chromA = False
		chromB = False

		# determine whether the gene is on each chromosome
		if gene_list[id_species][i] in chromA_file:
			chromA = True
		if gene_list[id_species][i] in chromB_file:
			chromB = True

		if chromA:
			if chromB:
				# some of the genes turned up on both chromosomes 
				gene_list[out_species][i] = "both"
			else: 
				# if only on 1st but not 2nd chromosome
				gene_list[out_species][i] = chromA_name

		elif chromB:
			# if only on 2nd but not 1st chromosome
			gene_list[out_species][i] = chromB_name

		else: 
			# this should never happen
			gene_list[out_species][i] = "NEITHER"

# when running the program, search the gene list for all 6 species
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

	# export the new database of gene locations
	gene_list.to_csv(path+'out/gene_loc_db4.csv', index=False)