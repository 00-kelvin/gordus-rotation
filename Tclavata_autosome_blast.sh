#!/bin/bash
# run from /Kelvin/autosome_blast/Tclavata directory

for i in {1..11}
do

	tblastn -subject ../../../sex_scaffs/Tclavata/Tclavata_chrom_${i}_cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out Tclavata_chrom_${i}_7sp.blastn.outfmt6
	tblastn -subject ../../../sex_scaffs/Tclavata/Tclavata_chrom_${i}_cds.fa -query ../../intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out Tclavata_chrom_${i}_8sp.blastn.outfmt6

	cut -f1 Tclavata_chrom_${i}_7sp.blastn.outfmt6 | sort | uniq > Tclavata_chrom_${i}_7sp_hits.txt
	cut -f1 Tclavata_chrom_${i}_8sp.blastn.outfmt6 | sort | uniq > Tclavata_chrom_${i}_8sp_hits.txt

done

cat *_7sp_hits.txt | sort | uniq > Tclavata_7sp_hits.txt
cat *_8sp_hits.txt | sort | uniq > Tclavata_8sp_hits.txt

sort ../../intersect_lists_sex_genes/sex_genes_USADLMT_mim_IDs.txt Tclavata_7sp_hits.txt | uniq -u > Tclavata_7sp_nohits.txt
sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Tclavata_8sp_hits.txt | uniq -u > Tclavata_8sp_nohits.txt

