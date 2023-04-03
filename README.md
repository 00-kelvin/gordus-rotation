# gordus-rotation

## 30.03.23

Getting started on the matrix of genes on same vs. different chromosomes in 5 species:

1. U. diversus
2. A. bruennichi
3. D. plantarius
4. L. elegans
5. M. bourneti

```
echo -e "#gene\tUdiv_chrom\tLele_chrom\tDplan_chrom\tMbour_chrom\tAbrue_chrom">sex_genes_location_db.txt
cat sex_genes_narrow.txt >> sex_genes_location_db.txt
```
Ok concerning... when I search for the sex genes in each of the scaffold files, nowhere near all 526 show up:

```
grep -f sex_genes_narrow.txt Udiv-Abru.anchors.unfiltered.txt | wc -l
  41
grep -f sex_genes_narrow.txt Udiv-Lelegans.anchors.no-interanchors.no-local.txt | wc -l
  57
grep -f sex_genes_narrow.txt Udiv-Dplan.anchors.unfiltered.txt | wc -l
  64
grep -f sex_genes_narrow.txt Udiv-Mbourn.anchors.unfiltered.txt | wc -l
  61
```

I could list these in my new file, but i feel like it makes more sense to just wait until i get a response from Andrew in case he just accidentally sent the wrong files

Update from andrew: "In the scaffold files, there are notation lines that say either “#block begin” or “#block end”. The data between these lines are anchor loci that were used to identify syntenic blocks, but they are not exhaustive of all the genes in the block. To identify genes in the block, use the first and last genes as indices."

Made scripts ```gene_loc_db.py``` and ```matrices.py```

Result from database script still showed lots of genese not found in all 5 species. 

3 problems with Andrews code: 

* Searching for "chr_1"
* Genes that fall into multiple blocks
* Also: "I intended to calculate the sex genes that overlapped in all 3 genomes (gene is in genome A AND B AND C), BUT my code actually calculated A OR B OR C. This explains why Calvin can’t find certain genes."

## 03.04.23

Changed ```re.search``` to ```re.fullmatch``` in Andrew's data to see if it helps. It definitely does not -- now nothing found anywhere because the sex scaffold IDs are not the full entry in that column. Seems like it would be an easy fix to just use the full IDs and then use fullmatch

Meeting with Jeremiah to figure out what to do about the data issue. 

* We can't use orthofinder yet because we don't have transcriptomes for M bourneti, D plantarius or L elegans
* Use liftover coordinates to make bed files of gene coordinates, then use ```bedtools getfasta``` to get coding sequences (cds) for those genes
	* Chat GPT script to extract chrStart and Stop for mRNA lines
	* convert to tab delimited
	* bedtools getfasta
	* nope, apparently this won't work either
* Jeremiah has another plan and will get back to me later this evening

To do:

* add a git ignore and update github repo
* email barbara about UTL 389 access

