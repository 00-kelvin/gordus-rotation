#!/bin/bash
# run from /Kelvin/autosome_blast/Dsilv directory

for i in {1,2,3,4,5,6,U1}
do

	tblastn -subject ../../../sex_scaffs/Dsilv/Dsilv_chrom_${i}.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Dsilv_chrom_${i}_8sp.blastn.outfmt6
	cut -f1 Dsilv_chrom_${i}_8sp.blastn.outfmt6 | sort | uniq > Dsilv_chrom_${i}_8sp_hits.txt

done

cat *_8sp_hits.txt | sort | uniq > Dsilv_8sp_hits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Dsilv_8sp_hits.txt | uniq -u > Dsilv_8sp_nohits.txt

