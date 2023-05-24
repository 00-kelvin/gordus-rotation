#!/bin/bash
# run from /Kelvin/autosome_blast/Mbourn directory

for i in {1..11}
do

	tblastn -subject ../../../sex_scaffs/Mbourn/Mbourn_chrom_${i}.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out Mbourn_chrom_${i}_7sp.blastn.outfmt6
	tblastn -subject ../../../sex_scaffs/Mbourn/Mbourn_chrom_${i}.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Mbourn_chrom_${i}_8sp.blastn.outfmt6

	cut -f1 Mbourn_chrom_${i}_7sp.blastn.outfmt6 | sort | uniq > Mbourn_chrom_${i}_7sp_hits.txt
	cut -f1 Mbourn_chrom_${i}_8sp.blastn.outfmt6 | sort | uniq > Mbourn_chrom_${i}_8sp_hits.txt

done

cat *_7sp_hits.txt | sort | uniq > Mbourn_7sp_hits.txt
cat *_8sp_hits.txt | sort | uniq > Mbourn_8sp_hits.txt

sort ../../intersect_lists_sex_genes/sex_genes_USADLMT_mim_IDs.txt Mbourn_7sp_hits.txt | uniq -u > Mbourn_7sp_nohits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Mbourn_8sp_hits.txt | uniq -u > Mbourn_8sp_nohits.txt

