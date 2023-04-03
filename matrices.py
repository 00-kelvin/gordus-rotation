#!/usr/bin/env python3

import numpy as np
import re
import pandas as pd

gene_list = open("/Users/cmdb/gordus-rotation/sex_genes_narrow.txt", 'r').read().split("\n")

diff_chrom_mat = pd.DataFrame(0, columns=gene_list, index=gene_list)

print(diff_chrom_mat)