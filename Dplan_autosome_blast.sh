#!/bin/bash
# run from /Kelvin/autosome_blast/Dplan directory

for i in {1,2,3,4,6,7,8,9,10,11,13,14}
do

	tblastn -subject ../../../sex_scaffs/Dplan/Dplan_chrom_${i}.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out Dplan_chrom_${i}_7sp.blastn.outfmt6
	tblastn -subject ../../../sex_scaffs/Dplan/Dplan_chrom_${i}.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Dplan_chrom_${i}_8sp.blastn.outfmt6

	cut -f1 Dplan_chrom_${i}_7sp.blastn.outfmt6 | sort | uniq > Dplan_chrom_${i}_7sp_hits.txt
	cut -f1 Dplan_chrom_${i}_8sp.blastn.outfmt6 | sort | uniq > Dplan_chrom_${i}_8sp_hits.txt

done

cat *_7sp_hits.txt | sort | uniq > Dplan_7sp_hits.txt
cat *_8sp_hits.txt | sort | uniq > Dplan_8sp_hits.txt

sort ../../intersect_lists_sex_genes/sex_genes_USADLMT_mim_IDs.txt Dplan_7sp_hits.txt | uniq -u > Dplan_7sp_nohits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Dplan_8sp_hits.txt | uniq -u > Dplan_8sp_nohits.txt

