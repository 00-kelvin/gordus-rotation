#!/bin/bash
# run from /Kelvin/autosome_blast/Udiv directory

for i in {1,2,4,5,6,7,8,9}
do

	tblastn -subject ../../../sex_scaffs/Udiv/Udiv.scaffold_${i}.cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out Udiv_scaff_${i}_7sp.blastn.outfmt6
	tblastn -subject ../../../sex_scaffs/Udiv/Udiv.scaffold_${i}.cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Udiv_scaff_${i}_8sp.blastn.outfmt6

	cut -f1 Udiv_scaff_${i}_7sp.blastn.outfmt6 | sort | uniq > Udiv_scaff_${i}_7sp_hits.txt
	cut -f1 Udiv_scaff_${i}_8sp.blastn.outfmt6 | sort | uniq > Udiv_scaff_${i}_8sp_hits.txt

done

cat *_7sp_hits.txt | sort | uniq > Udiv_7sp_hits.txt
cat *_8sp_hits.txt | sort | uniq > Udiv_8sp_hits.txt

sort ../../intersect_lists_sex_genes/sex_genes_USADLMT_mim_IDs.txt Udiv_7sp_hits.txt | uniq -u > Udiv_7sp_nohits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Udiv_8sp_hits.txt | uniq -u > Udiv_8sp_nohits.txt

