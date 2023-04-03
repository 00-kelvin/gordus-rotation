#!/usr/bin/env python3

import numpy as np
import re
import pandas as pd
import csv

path = '/Users/cmdb/gordus-rotation/'

# Open files
filename = path+'/Udiv-Lelegans.anchors.no-interanchors.no-local.txt'
with open(filename) as f:
    reader = csv.reader(f, delimiter="\t")
    Udiv_Lelegans_data = list(reader)

filename = path+'/Udiv-Dplan.anchors.unfiltered.txt'
with open(filename) as f:
    reader = csv.reader(f, delimiter="\t")
    Udiv_Dplan_data = list(reader)

filename = path+'/Udiv-Mbourn.anchors.unfiltered.txt'
with open(filename) as f:
    reader = csv.reader(f, delimiter="\t")
    Udiv_Mbourn_data = list(reader)

filename = path+'Udiv-Abru.anchors.unfiltered.txt'
with open(filename) as f:
    reader = csv.reader(f, delimiter="\t")
    Udiv_Abru_data = list(reader)

filename = path+'/sex_genes_narrow.txt'
with open(filename) as f:
    reader = csv.reader(f, delimiter="\t")
    sex_list = list(reader)

# Sex scaffold IDs

#this is an issue
Lelegans_sex = 'chr_1','chr_9'
Dplan_sex = 'chrom-5','chrom-12'
Mbourn_sex = 'X1', 'X2'
#what should this be?
Abru_sex = 'CM', 'CM'

def make_gene_list(gene_start, gene_end):

  genes = list(np.arange(gene_start,gene_end+1))

  for count, x in enumerate(genes):
    genes[count] = 'mRNA'+str(genes[count])

  return genes


def search_scaffold_file(sex_gene, scaffold_file, sex_scaffolds):

  line = 0
  scaffold_number = 'not found'
  Udiv_chrom = 'not found'

  while (line < len(scaffold_file)):

    if bool(re.fullmatch('#block begin', scaffold_file[line][0])):
      if bool(re.fullmatch('scaffold_3', scaffold_file[line+1][0])) or bool(re.fullmatch('scaffold_10', scaffold_file[line+1][0])):
        if bool(re.search(sex_scaffolds[0], scaffold_file[line+1][3])) or bool(re.search(sex_scaffolds[1], scaffold_file[line+1][3])):

          gene_start = scaffold_file[line+1][7]
          gene_start = int(gene_start[4:])
      
    if bool(re.fullmatch('#block end', scaffold_file[line][0])):
      if bool(re.fullmatch('scaffold_3', scaffold_file[line-1][0])) or bool(re.fullmatch('scaffold_10', scaffold_file[line-1][0])):
        if bool(re.search(sex_scaffolds[0], scaffold_file[line-1][3])) or bool(re.search(sex_scaffolds[1], scaffold_file[line-1][3])):

          gene_end = scaffold_file[line-1][7]
          gene_end = int(gene_end[4:])

          gene_list = make_gene_list(gene_start, gene_end)

          #this might not be working right. but also i think there may be some overlap
          if sex_gene in gene_list:
            scaffold_number = scaffold_file[line-1][3]
            Udiv_chrom = scaffold_file[line-1][0]

    line+=1

  return scaffold_number, Udiv_chrom

# def make_gene_dict(sex_gene):
#   gene_dict = dict.fromkeys(['gene', 'Udiv_chrom', 'Lele_chrom', 'Dplan_chrom', 'Mbour_chrom', 'Abrue_chrom'])
#   gene_dict['gene'] = sex_gene
#   gene_dict['Udiv_chrom'] = search_scaffold_file(sex_gene, Udiv_Dplan_data)[1]
#   gene_dict['Lele_chrom'] = search_scaffold_file(sex_gene, Udiv_Lelegans_data)[0]
#   gene_dict['Dplan_chrom'] = search_scaffold_file(sex_gene, Udiv_Dplan_data)[0]
#   gene_dict['Mbour_chrom'] = search_scaffold_file(sex_gene, Udiv_Mbourn_data)[0]
#   gene_dict['Abrue_chrom'] = search_scaffold_file(sex_gene, Udiv_Abru_data)[0]

#   return(gene_dict)

def make_gene_row(sex_gene):
  gene_row = []
  gene_row.append(sex_gene)

  #list of just u div chr?
  gene_row.append(search_scaffold_file(sex_gene, Udiv_Dplan_data, Dplan_sex)[1])
  gene_row.append(search_scaffold_file(sex_gene, Udiv_Lelegans_data, Lelegans_sex)[0])
  gene_row.append(search_scaffold_file(sex_gene, Udiv_Dplan_data, Dplan_sex)[0])
  gene_row.append(search_scaffold_file(sex_gene, Udiv_Mbourn_data, Mbourn_sex)[0])
  gene_row.append(search_scaffold_file(sex_gene, Udiv_Abru_data, Abru_sex)[0])

  return(gene_row)

if __name__ == "__main__":

  rows = []

  for i, gene in enumerate(sex_list):
    row = make_gene_row(gene[0])
    rows.append(row)

  df = pd.DataFrame(rows, columns=['gene', 'Udiv_chrom', 'Lele_chrom', 'Dplan_chrom', 'Mbour_chrom', 'Abrue_chrom'])
  df.to_csv('gene_loc_db3.csv', sep='\t', index=False)