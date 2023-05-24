#!/bin/bash
# run from /Kelvin/autosome_blast/Lele directory

for i in {2,3,4,5,6,7,8,10,11,12,13,14}
do

	tblastn -subject ../../../sex_scaffs/Lelegans/Lele_chrom_${i}_cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out Lele_chrom_${i}_7sp.blastn.outfmt6
	tblastn -subject ../../../sex_scaffs/Lelegans/Lele_chrom_${i}_cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Lele_chrom_${i}_8sp.blastn.outfmt6

	cut -f1 Lele_chrom_${i}_7sp.blastn.outfmt6 | sort | uniq > Lele_chrom_${i}_7sp_hits.txt
	cut -f1 Lele_chrom_${i}_8sp.blastn.outfmt6 | sort | uniq > Lele_chrom_${i}_8sp_hits.txt

done

cat *_7sp_hits.txt | sort | uniq > Lele_7sp_hits.txt
cat *_8sp_hits.txt | sort | uniq > Lele_8sp_hits.txt

sort ../../intersect_lists_sex_genes/sex_genes_USADLMT_mim_IDs.txt Lele_7sp_hits.txt | uniq -u > Lele_7sp_nohits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Lele_8sp_hits.txt | uniq -u > Lele_8sp_nohits.txt

