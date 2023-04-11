# gordus-rotation 

## 11.04.23

Smim.sex.prot.anns -> sort to only the 534 genes in sex_genes_div_mim_brue.txt

* Made a copy of the txt file of IDs in /media/will/mimic/Kelvin directory
* made a new file of only the Smim IDs:

```cut -f 2 -d$'\t' sex_genes_div_mim_brue.txt > sex_genes_mim_IDs.txt```

* copied the Smim.sex.prot.anns file to Kelvin directory
* filtered the protein sequences from the list of 534 using ```seqkit```:

```seqkit grep -f sex_genes_mim_IDs.txt Smim.sex.prot.anns -o sex_genes_div_mim_brue.fa```

* checked that all 534 are there using ```seqkit stats```

* Downloaded L.elegans transcripts and annotation file to ```sex_scaffs/Lelegans```, from http://gigadb.org/dataset/view/id/102210/File_page/2: will need to pull out genes that are annotated as on sex chromosomes and BLAST against those

BLAST command:
```
tblastn -subject [TRANSCRIPTOME/GENOME] -query [GENE LIST] -evalue 1e-5 -num_threads 36 [FOR TRANSCRIPTOMES:] -max_target_seqs 5 -out [FILENAME.blastn.out]
```

1. L elegans: chrom 1 and 9 (cds)
	* have txpts and ann files, need to filter to the ones on sex chroms
2. D. plantarius (chrom 5 and 12)
3. M. bourneti (chrom X1 and X2)
4. T. clavata chrom 12 and 13 (cds)
	* have the genome downloaded and split into chroms
	* have the cds and annotation file downloaded but cant get the files to unzip
5. D. silvatica (pending updated data -- email them?)

Used ```seqkit grep -p 'chr_1\.' -r Latrodectus_elegans_EVM.out.gff3.cds -o ../../Kelvin/Lele_chrom_1_cds.fa``` to get genes from L.elegans chrom 1 and 9, use -r flag to allow partial match, verified that chroms 10-14 were NOT included 

Going to go ahead and try BLASTing against the L.elegans cds, wish me luck. Ran these commands from /Kelvin directory:

```tblastn -subject Lele_chrom_1_cds.fa -query sex_genes_div_mim_brue.fa -evalue 1e-5 -num_threads 36 -max_target_seqs 5 -out ./sex_genes_blast/sex_genes_Lele_chrom_1.blastn.out```

```tblastn -subject Lele_chrom_9_cds.fa -query sex_genes_div_mim_brue.fa -evalue 1e-5 -num_threads 36 -max_target_seqs 5 -out ./sex_genes_blast/sex_genes_Lele_chrom_9.blastn.out```

```tblastn -subject ../sex_scaffs/Dplan/Dplan_chrom_5.fa -query sex_genes_div_mim_brue.fa -evalue 1e-5 -num_threads 36 -out ./sex_genes_blast/sex_genes_Dplan_chrom_5.blastn.out``` and same for chrom 12

```tblastn -subject ../sex_scaffs/Mbourn/Mbourn_chrom_X1.fa -query sex_genes_div_mim_brue.fa -evalue 1e-5 -num_threads 36 -out ./sex_genes_blast/sex_genes_Mbourn_chrom_X1.blastn.out``` and same for X2

## 10.04.23

Instructions from Jeremiah: 

> We have genome, transcriptome, and annotation files for Uloborus diversus and Argiope bruennichi.
We know that these data exist for Latrodectus elegans, so we will want to check that we have them or get whatever we don't already have. This one is in GigaDB, and not in NCBI, so that's where we'll need to go to get this data (and my experience with this is that it can be incredibly slow).
We know that the genome and transcriptome of Hylyphantes graminicola are available. We need to check on the annotation file, or see if we can generate one by mapping the transcriptome to the genome.
There are no transcriptomes that we've been able to find for Meta bourneti or Dolomedes plantarius, so we'll have to count those out for Anchorwave, but once we have a set of candidate genes, we can still map those to these species to see which of the proposed X chromosomes they fall onto and whether the pattern of segregation is the same.
We need to find out if we can get transcriptomes and annotation files for Dysdera silvatica and Metellina segmentata... especially D. silvatica, since that is the one with the single X chromosome.
I think you mentioned that there were two Tetragnatha species and Trichonephila clavipes genomes available, at least one of which was chromosome-scale assembly, so we'd want to look for the associated files for those as well.

>To get these files in their proper place, navigate to:
/media/will/mimic/sex_scaffs/
and create a folder (e.g, Hgram) for the species. Then, you can just use wget to get the files in that folder.
After we have the multi-fasta with the assembly in it, we'll want to break it down into chromosomes or scaffolds.
You can filter out the smaller scaffolds first, using a one-liner like the following:
bioawk -c fastx '{ if(length($seq) > 750000) { print ">"$name; print $seq }}' species_assembly.fa > species_assembly.chroms.fa
Then, you can break that multi-fasta into single fasta files, each with a single chromosome in it.

>If I recall correctly, you can use seqkit split for that task. Then, you'll want to rename the files using the naming convention that we've used for the other species, which looks like [species_abbreviation]dot[chrom_X].fa
Or it might be [species_abbreviation]_[chrom_X].fa
Something like that just to make it easier for us to keep track of

* Downloaded _Hylyphantes graminicola_ genome assembly

```
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_023701765.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT&filename=GCA_023701765.1.zip" -H "Accept: application/zip"

unzip GCA_023701765.1.zip 
```

* Filtered out smaller scaffolds with 

```
bioawk -c fastx '{ if(length($seq) > 750000) { print ">"$name; print $seq }}' GCA_023701765.1_ASM2370176v1_genomic.fna > ../../../Hgram_assembly.chroms.fa
```

* Split into chromosomes with 

``` 
seqkit split -i Hgram_assembly.chroms.fa
```

* Renamed according to naming convention [species_abbreviation]_[chrom_X].fa
* Found annotation file at https://www.scidb.cn/en/detail?dataSetId=b2a377bee8184231bbf8bee49bd98928

``` 
wget https://download.scidb.cn/download?fileId=61aec71689f14b48842fbae7&dataSetType=personal&fileName=IOZCAS_Hgram_genomeAssembly_1.0.gff
```

* Renamed Hgram_genomeAssembly_1.0.gff 
* Moved transcriptome file (?) /media/will/mimic/tsas/cd-hit_complete/Hylyphantes_graminicola.fsa_nt to sex_scaffs/Hgram folder

### Meeting with Andrew and Jeremiah

* average number of transcripts in orthofinder call per species is ~15,000
* Verify that BUSTED and RELAX pipeline work using example data to make sure we can reproduce their results (Jeremiah)
* Me: BLAST the sex-linked genes in sex chromosomes
	* Use annotation file to pull out transcripts that are on the sex chromosomes
	* _L. elegans_ has a transcriptome, should try to get that downloaded
	* BLAST against the _D. silvatica_ big sex chromosome, but just the transcripts that come from that chromosome
	* _M. bourneti_ and _D. plantarius_: do a BLAST search against each of the sex chromosomes
	* restrict outputs to e value (second to last column) of 1E-5 or less
	* Last column is the bit score: bigger = better

```	
blastn -subject [TRANSCRIPTOME/GENOME] -query [GENE LIST] -evalue 1e-5 -num_threads 36 [FOR TRANSCRIPTOMES:] -max_target_seqs 5 -out [FILENAME.blastn.out]
```
 
* Looked for D. sylvatica genome and annotation files -- there's a pre-print about chromosome-level assembly, but all the links are to the previous scaffold-level version: https://www.ncbi.nlm.nih.gov/assembly/GCA_006491805.2
* Downloaded T. clavata genome; ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9932165/ indicates sex chromosomes are 12 and 13

```
bioawk -c fastx '{ if(length($seq) > 750000) { print ">"$name; print $seq }}' dna.all.fa > Tclavata_assembly.chroms.fa

seqkit split -i Tclavata_assembly.chroms.fa

```

Waiting for someone to send me that list of 534 genes. In the meantime, sending H.gram through the following pipeline to prep it for orthofinder:

1. TransDecoder.LongOrfs

```TransDecoder.LongOrfs -t Hylyphantes_graminicola.cd-hit-est -O ../../transdecoder/clustered/Hylyphantes_graminicola```

2. blastp

run from the blastp folder:

```blastp -subject ../TD-LO_complete/Hylyphantes_graminicola.cd-hit-est -query /media/will/demogorgon/refdbs/uniprot/uniprot_sprot.fasta -evalue 1e-5 -num_threads 36 -max_target_seqs 1 -outfmt 6 -out Hylyphantes_graminicola.cd-hit-est.blastp.out```

3. pfam

run from pfam directory:

```
hmmscan --cpu 36 --domtblout Hylyphantes_graminicola.cd-hit-est.pfam.domtblout /media/will/demogorgon/refdbs/Pfam/Pfam-A.hmm /media/will/mimic/transdecoder/clustered/TD-LO_complete/Hylyphantes_graminicola/Hylyphantes_graminicola.cd-hit-est.transdecoder_dir/longestest_orfs.pep
```

4.  TransDecoder.Predict

## 07.04.23

Finished reading Jeremiah's paper. Read up a bit on HyPhy

## 06.04.23

Searched for transcriptomes and chromosome-level genome assemblies

#### Transcripts

* _Metellina segmentata_: not found
* _Dysdera silvatica_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA313901/
* _Dolomedes plantarius_: not found
* _Meta bourneti_: not found

#### Chromosome-level assemblies

* _Trichonephila clavipes_: not found
* _Trichonephila antipodiana_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA627506
* _Trichonephila clavata_: https://spider.bioinfotoolkits.net/main/species-info/-1

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

