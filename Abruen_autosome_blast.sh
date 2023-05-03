#!/bin/bash
# run from /Kelvin/autosome_blast/Abruen directory

for i in {1,2,3,4,5,6,7,8,11,12,13}
do

	tblastn -subject ../../../sex_scaffs/Abruen/Abruen.chrom_${i}.cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out Abruen_chrom_${i}_7sp.blastn.outfmt6
	tblastn -subject ../../../sex_scaffs/Abruen/Abruen.chrom_${i}.cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Abruen_chrom_${i}_8sp.blastn.outfmt6

	cut -f1 Abruen_chrom_${i}_7sp.blastn.outfmt6 | sort | uniq > Abruen_chrom_${i}_7sp_hits.txt
	cut -f1 Abruen_chrom_${i}_8sp.blastn.outfmt6 | sort | uniq > Abruen_chrom_${i}_8sp_hits.txt

done

cat *_7sp_hits.txt | sort | uniq > Abruen_7sp_hits.txt
cat *_8sp_hits.txt | sort | uniq > Abruen_8sp_hits.txt

sort ../../intersect_lists_sex_genes/sex_genes_USADLMT_mim_IDs.txt Abruen_7sp_hits.txt | uniq -u > Abruen_7sp_nohits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Abruen_8sp_hits.txt | uniq -u > Abruen_8sp_nohits.txt

