# gordus-rotation 


## 12.05.23

### Are the hox genes in the "linkage blocks"? 

mRNA6585
mRNA7140
mRNA7638
mRNA7683
mRNA7687
mRNA7926
mRNA16425
mRNA18127
mRNA18588
mRNA19316
mRNA19326
mRNA20134
mRNA20147


Linkage block 1A.1: mRNA16056 -> mRNA16514

* mRNA16425

Linkage block 1A.2: mRNA20258 -> mRNA20263

Linkage block 1B.1: mRNA5846 -> mRNA6591

* mRNA6585

Linkage block 1B.2: mRNA7857 -> mRNA8047

* mRNA7926

Only 3 of the conserved HOX genes are found in the 3 biggest "linkage blocks" and none of them in the same ones

### Are the sex-chrom-only hits in the linkage blocks? 

21 of the 59 genes from the 7species intersection that have no autosome hits are found in the 3 biggest "linkage blocks"

Lots of them are "consensus disorder prediction" or just "coil"

## 08.05.23

Drosophila gene sequences in ```/demogorgon/refseqs/neuro-targets/Dmel``` somewhere

Any receptors that are not already in that folder of those against Udiv annotated transcripts, download from FlyBase


Results of blast against Udiv CDS already done for the ones Jeremiah already has downloaded: ```/mimic/blast/neuro-targets/dovetail-002/braker-out-001```

but we want to get Udiv protein sequences, so take those fly fasta files and blastp them against the Udiv genome --> get protein sequences
--> check those sequences for the top hits, see how they cluster? 

Next step is to concatenate the hits from Udiv and BLAST those against each other spider species

```xclip -sel c``` (command to copy from a file)

### Orthofinder meeting

* 1st run: 6 species, Transdecoder > CD-HIT, only ~300 orthogroups in all species (should be more like ~1500 BUSCOs)
* 2nd run: 21 species just using BUSCO transcripts, still only giving ~300 orthogroups
* 3rd run: using protein species since those should be more conserved, recovered almost 1000
* 2 affects of saturated sites: alignment will be hard; also could be 2 mutations A>T>A
* Insertions and reversions are much easier to capture in protein data
* Orthologs are breaking apart into multiples because the program can easily be thrown off by large insertions
* Recommend: find the transcripts corresponding to orthogroups according to protein sequences and do alignment
* Typically alignment will look like a very conserved regions and then a bunch of chaotic regions -- could identify these regions and use them to test hypotheses
* Have to think about what claim you're making: overall positive selection or positive selection for this domain?
* Would it make the most sense to use genes that are actually related to web weaving function?
* Possible to do reciprocal BLAST of the protein sequence for a small set of genes that are relevant -- GCPRs, ion channels relevant to neurons, dopamine receptors? Use drosophila as an outgroup
* Followup: where are these genes expressed in orb weavers vs others? knock out those genes and see how it affects orb weaving behavior

## 07.05.23

* Stuff to potentially do later:
	* Add the new chrom-level assemblies to analysis
	* More on dsx TF binding sites

### Hox genes on sex chromosomes

List of genes from 8-species intersect that are coming up as hox genes (or at least having homeobox domain): 

mRNA7638
mRNA7687
mRNA18127
mRNA19316
mRNA19326
mRNA20134
mRNA20147
mRNA7683
mRNA7926
mRNA16425
mRNA18588
mRNA6585
mRNA7140



### Doublesex on sex chromosomes

* Udiv: 2, 7 and **10**
* Abruen: 8, **9, 10**, 11
* Lele: **1** and 7
* Tclav: copies on NCBI protein database didn't list chromosomes
	* Did BLAST of one of their sequences against my Tclav cds, had hits on chroms 10, 11, **12 and 13**
* Dplan: nothing comes up on NCBI (no transcriptome)
	* Blasted the Tclav dsx sequence against the sex chrom genome seqences and there was a hit on both **5 and 12**
* Mbourn: nothing comes up (no transcriptome)
	* Did same thing as Dplan and had a hit on chrom **X1** but not X2
* Dsilv: nothing comes up (transcriptome exists but predates chrom-level assembly)
	* Did same thing as Dplan/Mbourn and had a hit on the big **X**
* nothing comes up for Stegodyphus mimosarum in NCBI, and no hits when I blast the sex scaffolds with the Tclav sequence (explains why dsx didn't show up in the intersect list..?)

###TF binding sites around dsx orthologs:

MotifMap (https://motifmap.ics.uci.edu) is a tool that does this but only for specific genomes

JASPAR (https://jaspar.genereg.net) can take in a TF sequence and give potential consensus sequence (didn't give anything for mRNA6646

This site has a lot of potentially useful databases: https://footprintdb.eead.csic.es/?databases 

### Syntenic blocks

I think these aren't exactly syntenic blocks because the list already has a lot of gaps in it from intersecting it across different species but w/e. Need to figure out what Andrew actually wants with this info

**List 1A.1**
mRNA16056
mRNA16058
mRNA16063
mRNA16065
mRNA16088
mRNA16092
[mRNA16101]
mRNA16187
[mRNA16248]
mRNA16250
mRNA16294
mRNA16297
mRNA16309
mRNA16344
[mRNA16347]
mRNA16475
[mRNA16482]
mRNA16489
mRNA16505
mRNA16514

**List 1A.2**
mRNA20258
mRNA20262
mRNA20263

vs

**List 1B.1**
mRNA5846	
mRNA5850	
mRNA5859	
[mRNA5980]	
mRNA6118	
mRNA6120	
mRNA6155	
mRNA6186	
mRNA6199	
[mRNA6240]	
[mRNA6247]	
mRNA6283	
mRNA6319	
mRNA6330	
mRNA6350	
[mRNA6382]	
mRNA6391	
mRNA6426	
mRNA6446	
mRNA6495	
mRNA6547	
mRNA6555	
mRNA6579	
mRNA6585	
mRNA6591

**List 1B.2**
mRNA7857	
mRNA7893	
mRNA7898	
mRNA7926	
mRNA7981	
[mRNA7992]	
[mRNA8005]	
mRNA8047

These are essentially my Blocks 1A and 1B from my clustermap. Messed up my ```diff_chrom_mat.csv``` file, will need to regenerate if I need to use it for something else later

## 04.05.23

With a blastp search, mRNA6646 is coming up as "Segmentation protein cap'n'collar" in a bunch of spider species --> has a bunch of different targets in different species https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7588902/ aging, proteasome, chromatin remodeling...

Potential consensus sequence for mRNA6646: RTGACTCAGCA
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4911266/)

ATGACTCAGCA
GTGACTCAGCA
TGCTGAGTCAC
TGCTGAGTCAT

none of these were in the 10kb region upstream of the dsx ortholog on chromosome 7

```
grep "ATGACTCAGCA" dsx_chr7_upstream.fa 
grep "GTGACTCAGCA" dsx_chr7_upstream.fa 
grep "TGCTGAGTCAC" dsx_chr7_upstream.fa 
grep "TGCTGAGTCAT" dsx_chr7_upstream.fa 
```

also checked the chromosome 2 ortholog, nothing there

There is also a doublesex on chr10 which has 2 isoforms

Could do a tree of bZIP TF proteins and see where ours lands

Find dsx in spider genome, go 10kb up and down, dump into txp factor data base to find binding sites/partners

The 7species intersect list of sex-only genes has 59 genes, 36 of them are from our blocks of seperately segregating genes (seems like a pretty high ratio out of the 363 genes from the original intersection)

* Look upstream of dsx orthologs for TF binding sides
* look for dsx on sex chroms of other spiders
* syntenic blocks of little subgene list segregating apart always from my original heatmap?
* in rotation talk mention that there are hox genes on the sex chroms of all species

## 03.05.23

Done: 

* U div, cds (all but 3 and 10)
* A bruen, cds (all but 9 and 10)

Left: 

* D plan, genome (all but 5 and 12)
* L ele, cds (all but 1 and 9)
* M bourn, genome (all but X1 and X2)
* T clavata, cds (all but 12 and 13)
* D silv, genome (all but X)

Had to make separate CDS files for the chromosomes for L.elegans and T.clavata:

```seqkit grep -p "Tc01G" -r cds.fa -o Tclavata_chrom_1_cds.fa```

or

```seqkit grep -p 'chr_10\.' -r Latrodectus_elegans_EVM.out.gff3.cds -o Lele_chrom_10_cds.fa```

**Tallying no-hit counts:**

8-species intersection

* U.div: 47
* A.bruen: 48
* D.plan: 69
* L.ele: 48
* M.bourn: 65
* T.clavata: 49
* D.silv: 73

7-species intersection

* U.div: 74
* A.bruen: 73
* D.plan: 108
* L.ele: 77
* M.bourn: 103
* T.clavata: 82

Keep in mind: there may be genes on other non-chromosomal scaffolds...

It was def a waste of time to do the 8- and 7-species intersections as separate BLAST runs instead of just doing the 7species and then filtering to the 8species list. but oh well, it's almost done

```cat *7sp_nohits.txt | sort-uniq-count-rank > no_hits_7sp_scur.txt```

```grep "6	" no_hits_7sp_sucr.txt | cut -f2 > no_hits_any_sp_7sp_mim_IDs.txt```

then used the ```grep -f``` against the ```sex_genes_div_mim_brue.txt``` to get IDs for div/bruen

6646 is one of them....... So is 19375

## 02.05.23

Created fasta files of sequences for the 7- and 8- species intersect lists.

```seqkit grep -f sex_genes_USADLMT_mim_IDs.txt ../../sex_scaffs/Smim/S_mimosarum_proteins.txt -o sex_genes_USADLMT.fa```

Test BLAST of sex genes against Udiv autosome 1

```tblastn -subject ../sex_scaffs/Udiv/Udiv.scaffold_1.cds.fa -query ./intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out ./autosome_blast/Udiv_scaff_1_7sp.blastn.outfmt6```

Also did with no evalue limit

```tblastn -subject ../sex_scaffs/Udiv/Udiv.scaffold_1.cds.fa -query ./intersect_lists_sex_genes/sex_genes_USADLMT.fa -outfmt 6 -out ./autosome_blast/Udiv_scaff_1_7sp_no-elim.blastn.outfmt6```

Now for rest of Udiv scaffolds, using the 7-species intersection, with and without e-value limits...

```for i in {2,4,5,6,7,8,9}; do tblastn -subject ../sex_scaffs/Udiv/Udiv.scaffold_${i}.cds.fa -query ./intersect_lists_sex_genes/sex_genes_USADLMT.fa -outfmt 6 -out ./autosome_blast/Udiv_scaff_${i}_7sp_no-elim.blastn.outfmt; tblastn -subject ../sex_scaffs/Udiv/Udiv.scaffold_${i}.cds.fa -query ./intersect_lists_sex_genes/sex_genes_USADLMT.fa -evalue 1e-5 -outfmt 6 -out ./autosome_blast/Udiv_scaff_${i}_7sp.blastn.outfmt; done```

... and the 8-species intersection

```for i in {1,2,4,5,6,7,8,9}; do tblastn -subject ../sex_scaffs/Udiv/Udiv.scaffold_${i}.cds.fa -query ./intersect_lists_sex_genes/sex_genes_USADLMTD.fa -outfmt 6 -out ./autosome_blast/Udiv_scaff_${i}_8sp_no-elim.blastn.outfmt; tblastn -subject ../sex_scaffs/Udiv/Udiv.scaffold_${i}.cds.fa -query ./intersect_lists_sex_genes/sex_genes_USADLMTD.fa -evalue 1e-5 -outfmt 6 -out ./autosome_blast/Udiv_scaff_${i}_8sp.blastn.outfmt; done```

Made a script to compress the BLAST hits into a unique list for each condition ```list_hits_Udiv.sh```

Created a list of genes that had hits on any of the autosomees in the 4 conditions: 
```cat *_7sp_hits.txt | sort | uniq > Udiv_7sp_hits.txt```

With no e-value limit, there are no genes in either the 7-species or 8-species list that had no hits on any of the Udiv autosomes. But with e-value limit of 1e-5, there are some genes with no hits. 

Made lists of genes that had no hits on the autosomes (with the e-value limit) for the 7- and 8-species intersections:

```sort ../../intersect_lists_sex_genes/sex_genes_USADLMTD_mim_IDs.txt Udiv_8sp_hits.txt | uniq -u > Udiv_8sp_nohits.txt```

Made a new script ```Abruen_autosome_blast.sh``` that does all the above steps for Abruen.

The counts of genes that had no autosome hits were almost identical for the U.div and A.bruen: 

8-species intersection

* U.div: 47
* A.bruen: 48

7-species intersection

* U.div: 74
* A.bruen: 73

Maybe this is a good sign? Or maybe it's due to some error of mine. We shall see!


## 27.04.23

Collecting gene annotations for the genes in the blocks that segregate separately.

In block 1A mRNA20493 and mRNA16063 don't have annotations

1 missing from block 1B: mRNA6350

None missing from block 2A or 2B

None missing from block 3A or 3B

BLASTed mRNA18921 sequence: only 2 results are uncharacterized genes in U.div and *Parasteatoda tepidariorum*, but maybe there are structural homologs in other spider species

Interesting candidate genes: 

* 6646: bZIP Mef-like TF (regulates doublesex in crustacean *Daphnia magna* https://doi.org/10.1371/journal.pgen.1006953)
* Segregating partners of 6646:
	* 18181: suppressor of RAS, has pleckstrin domain, downstream targets of SRY
	* 18895: related to chromatin organization
	* 18921: uncharacterized
	* 18991: ribonuclease inhibitor, interesting?
	* 19375: MA3 domain might be related to sex determination
	* 19979: TF CoE-1
	* 20009: LIM-domain-BP, repro dev https://doi.org/10.1159/000518323
	* 20021: SH3-domain-BP

Also got list of annotations for the genes that were always on both chromosomes

### Meeting with Andrew

* look for BLAST hits of these same genes on autosomes -- if there are any that are only on sex chromosomes those would be interesting
* for any TFs in our candidate list, if anything known about consensus sequences, look for genes on other chromosomes that are their targets
* look at known dsx orthologs and build a tree with those + our dsx gene family

## 26.04.23

Count of genes found on both chroms for all 4 added species: 70
Fully segregating pairs: 1061
Count of genes involved in fully segregating pairs: 150

Made heat map and histogram

Clustered the genes to find blocks of genes that are always segregating apart from other blocks of genes, saved lists in blocks.csv

## 24.04.23

Fixed the ```matrices.py``` script, now I have my 2 matrices, woo! 

* Find out how many of the genes are on both chromosomes for all species
* How many pairs are invovled and how many genes are in pairs (rows / 2)
* Check out go annotations for these genes
* Make a heat map of the matrix
* Histogram: for each row, number of columns that have 6s

## 20.04.23

Have to redo the intersections because I didn't make the lists for each species unique before narrowing down, so the ```uniq -d``` command also included ones that were just on both chromosomes for one species

**NEW COUNTS:** 

* Start: 534
* add Dplan: 494
* add Lele: 386
* add Mbourn: 380
* add Tclavata: 363
* add Dsilv: 275

```
 2027  grep -f ./intersect_lists_sex_genes/sex_genes_USADL.txt sex_genes_div_mim_brue.txt >> sex_genes_USADL_ids.txt 
 2033  grep -f ./intersect_lists_sex_genes/sex_genes_USADLM.txt sex_genes_div_mim_brue.txt >> sex_genes_USADLM_ids.txt 
 2053  grep "KAF" Abruen.chrom_9.CDS
 2054  grep "KAF" Abruen.chrom_9.longestCDS.ids 
 2076  grep -f ./intersect_lists_sex_genes/sex_genes_USADLMT.txt sex_genes_div_mim_brue.txt >> sex_genes_USADLMT_ids.txt 
 2105  grep -f sex_genes_USADLMT_div_IDs.txt Udiv.scaffold_3.cds.fa > Udiv_scaff_3_sex_genes.txt
 2108  grep -f sex_genes_USADLMT_div_IDs.txt Udiv.scaffold_10.cds.fa > Udiv_scaff_10_sex_genes.txt
 2110  grep -f sex_genes_USADLMT_bruen_IDs.txt Abruen.chrom_9.gff > Abruen_chrom_9_sex_genes.txt
 2113  grep -f sex_genes_USADLMT_bruen_IDs.txt Abruen.chrom_10.gff > Abruen_chrom_10_sex_genes.txt
```

^ used that to make the lists of IDs for the different intersections

This solved the problem of getting genes that were on neither chromosome for certain species once I fed it back in to my ```gene_loc_db2.py``` script.

Have my ```matrices.py``` up and running but the output looks wrong, need to figure out how to increment entries correctly and also what to do about genes that are on both chromosomes.

## 17.04.23

Sent Andrew the new narrowed-down gene lists. 

Now, need to sort them by chromosome

A. bruennechi chroms 9 and 10, U. div chromosomes 3 and 10

To sort U.div genes by which chromosome they're from: 

* made a file of just Udiv IDs for the intersection between the 7 species with XX/XXXX

```cut -f 1 sex_genes_USADLMT_ids.txt | tail -n +2 > sex_genes_USADLMT_div_IDs.txt```

* search in Udiv scaffold 3 and 10 cds fasta files for those genes

```grep -f sex_genes_USADLMT_div_IDs.txt Udiv.scaffold_3.cds.fa > Udiv_scaff_3_sex_genes.txt```

* 215 were on scaffold 3, 161 were on scaffold 3!
* Same for A. bruen, but have to use Abruen.chrom_X.gff files, which have a bunch of other stuff in them...
* For the other four species, use the blastn.outfmt6 files

Exported all those files to local computer and wrote a script in ```matrices.py``` to collect which chromosome each one is on. Some of them are showing up as on "neither" chromosome, which is bad, need to figure out why that is.

## 13.04.23

Problem: the Smim.sex.prot.anns file only had the first 60 aa of each sequence. Recreated list using S_mimosarum_proteins.txt file which has full length sequences

```seqkit grep -f sex_genes_mim_IDs.txt ../sex_scaffs/Smim/S_mimosarum_proteins.txt -o sex_genes_div_mim_brue.fa```

Rerun BLAST using full sequences. This time, limiting to max target sequence of 1 (for CDS) and max "high scoring pairs" (HSPs) of 1 as well, so that each query only has 1 alignment per chromosome. just create a table of outputs. 

Examples:

```tblastn -subject Lele_chrom_9_cds.fa -query sex_genes_div_mim_brue.fa -evalue 1e-5 -max_target_seqs 1 -max_hsps 1 -outfmt 6 -out ./sex_genes_blast/sex_genes_Lele_chrom_9.blastn.outfmt6_max-hits-1_max-hsps-1```

```tblastn -subject ../sex_scaffs/Dplan/Dplan_chrom_5.fa -query sex_genes_div_mim_brue.fa -evalue 1e-5 -max_hsps 1 -outfmt 6 -out ./sex_genes_blast/sex_genes_Dplan_chrom_5.blastn.outfmt6_max-hits-1_max-hsps-1```

**Count of hits:** 

* L.ele chrom 1: 290
* L.ele chrom 9: 257
	* L.ele total unique: 408
* D.plan chrom 5: 336
* D.plan chrom 12: 272
	* D.plan total unique: 494
* M.bourn chrom X1: 340
* M.bourn chrom X2: 295
	* M.bourn total unique: 517
* T.clavata chrom 12: 321
* T.clavata chrom 13: 333
	* T.clavata total unique: 499
* D.silvatica chrom X: 358

**To add:** T. clavata, D. silvatica

Got the T.clav files "unzipped" (they weren't actually tar zipped)

Finally found the chromosome-level D.silvatica genome: https://www.ncbi.nlm.nih.gov/Traces/wgs/QLNU02?display=contigs

Downloaded the assembly as 2 separate files, then used bioawk to filter to scaffolds >20Mbp (since those were the pseudochromosomes 1-6, U1 and X listed in the paper https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13471)

Split using seqkit, then sorted by size and renamed each scaffold to their chromosome names according to lengths indicated in the paper. The largest one at >300Mbp is the X chromosome.

Also found and downloaded transcripts and annotation file for D.silvatica (ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5381604/), but these are unlikely to be annotated according to which chromosome they are from, since the transcriptome data predates the chromosome-level assembly. Might just BLAST against the X chromosome genomic data.

TO DO: 

* BLAST against T. clavata and D. silvatica
* Filter list of 534 to genes that had hits on all 5 other species
* Make a spreadsheet for keeping track of links that I've been collecting for different species

Pattern for T.clavata cds is TcXX where XX = chrom #, use that to make separate fasta files for Chrom 12 and 13 cds.

Ran BLAST on T.clavata cds and D.silv big X chromosome.

Created text files for each species containing all gene IDs with hits for that species (cat the 2 files for those with 2 x chromosomes, pipe to uniq) however, even when piping to uniq the combined files still had >534 genes. Don't know what that's about but extras should be eliminated once these lists are intersected with the original list

[this was just because i didn't sort before using uniq]

Ok, intersected sequentially using ```sort sex_genes_mim_IDs.txt Dplan_sex_genes.txt | uniq -d > sex_genes_USAD.txt``` in this order:

* Start: 534 (Udiv + Smim + Abreun)
* add Dplan: 494
* add Lele: 387
* add Mbourn: 386
* add Tclavata: 376
* add Dsilv: 285

final list is called ```sex_genes_USADLMTD.txt```

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

