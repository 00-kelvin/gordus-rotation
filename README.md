# gordus-rotation

## 06.04.23

Searched for transcriptomes and chromosome-level genome assemblies

#### Transcripts
_Metellina segmentata_: not found
_Dysdera silvatica_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313901/
_Dolomedes plantarius_: not found
_Meta bourneti_: not found

#### Chromosome-level assemblies
_Trichonephila clavipes_: not found
_Trichonephila antipodiana_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA627506
_Trichonephila clavata_: https://spider.bioinfotoolkits.net/main/species-info/-1

## 05.04.23

### Selection project: 

Looking for any genes under positive selection in orb weavers or relaxed selection in non-orb weavers

21 species prepped and currently being run in orthofinder to generate multiple sequence alignments

* should not produce an alignment if there is >50% gap in column (of nucleotides)  

Then --> concatenate, clean and trim (need to have no stop codons, be multiple of 3) using Mega

Then --> HyPhy

will/mimic/hyphy/SequenceTreeAlignment.fa

Try running Hyphy on the 6 species test files; see what errors it throws, then try to use Mega to edit the files according to the errors

### Sex-linked genes project: 

Need transcriptome data to plug them into Anchorwave. Don't have those for _D. plantarius_ or _M. bourneti_. BUT might exist for _H. graminicola_ and _M. segmentata_. Already done for _U. div-A. bruen_. Also have _S. min_ scaffold-level assembly and list of sex-linked genes -- won't be able to tell if these gene pairs are on same or different chromosomes, but can still use them to come up with the list of genes.

We do have transcriptome data for _H. graminicola_! 
Search around to see if there is transcriptome data for _M. segmentata_ and _D. silvatica_, can also check again on _D. plantarius_ and _M. bourneti_ around the internet. 

Also look for genome assemblies for _Trichonephila clavipes_ and _Tetragnatha_ spp (chromosome-scale) and any other spider chromosome-scale assemblies

_D. silvatica_ has one big X chromosome? 

Even though _D. plantarius_ and _M. bourneti_ don't have transcripts, can do a BLAST search with genes we come up with from the other 4-6

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

