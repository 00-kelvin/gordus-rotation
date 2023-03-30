# gordus-rotation

## 30.03.23

Getting started on the matrix of genes on same vs. different chromosomes in 5 species:

1. U. diversus
2. A. bruennichi
3. D. plantarius
4. L. elegans
5. M. bourneti

	echo -e "#gene\tUdiv_chrom\tLele_chrom\tDplan_chrom\tMbour_chrom\tAbrue_chrom">sex_genes_location_db.txt
	cat sex_genes_narrow.txt >> sex_genes_location_db.txt

Ok concerning... when I search for the sex genes in each of the scaffold files, nowhere near all 526 show up:

	grep -f sex_genes_narrow.txt Udiv-Abru.anchors.unfiltered.txt | wc -l
      41
	grep -f sex_genes_narrow.txt Udiv-Lelegans.anchors.no-interanchors.no-local.txt | wc -l
      57
	grep -f sex_genes_narrow.txt Udiv-Dplan.anchors.unfiltered.txt | wc -l
      64
	grep -f sex_genes_narrow.txt Udiv-Mbourn.anchors.unfiltered.txt | wc -l
      61

