###########################################################################################################
2023
leanie
omics 
nthi diversity panel
###########################################################################################################
# establishing diversity strain panel
Collaborators collected NTHi samples from healthy and otitis media-prone children during healthy or an occurrence of OM. Culturing followed by genomic sequencing identified patterns of transience, persistence, and recurrence in healthy and sick children.  Knowing that haemophilus influenzae has large genetic diversity across the species, we aimed to identify common biological signatures which are shared among strains however were capable of distinguishing _in vitro_ phenotypes which associate with commensals or pathogens.

A diversity panel of NTHi strains was created from a database of genome sequences from clinical isolates. A gene presence/absence matrix was determined. Computationally, we selected a set of strains which would maximize the amount of genetic diversity (measured by number of genes present within strains found at a frequency greater than ten percent) in the fewest number of strains. This resulted in 9 clinical isolates which captured genetic diversity. In addition, 3 well-studied controls were included in the “diversity panel”.

| strain   | reference_genome | b_num  | bs_num  | hs_num     | patient    | dob        | visit_date | visit_type | months_old | body_site | op      | healthy | has_follow_up | is_follow_up |     |     |
| -------- | ---------------- | ------ | ------- | ---------- | ---------- | ---------- | ---------- | ---------- | ---------- | --------- | ------- | ------- | ------------- | ------------ | --- | --- |
| Rd KW20  | Rd KW20          |        |         |            |            |            |            |            |            |           |         |         |               |              |     |     |
| 86 028Np | 86 028Np         |        |         |            |            |            |            |            |            |           |         |         |               |              |     |     |
| R2866    | R2866            |        |         |            |            |            |            |            |            |           |         |         |               |              |     |     |
| SD_1763  |SD_1763| B5880            | BS2182 | HS52838 | 14-02-604  | 11/12/2014 | 3/1/2016   | AOM6       | 15         | RMEFB      | NOP       | sick    | no      | FALSE         |              |     |     |
| SD_1474  |SD_1474 | B5916            | BS2218 | HS52874 | 12-02-0447 | 7/20/2012  | 7/3/2013   | AOM5       | 11         | RMEF       | OP        | sick    | yes     | FALSE         |              |     |     |
| SD_1430  |SD_1430 | B5930            | BS2232 | HS52888 | 12-02-0433 | 4/22/2012  | 3/14/2013  | AOMTF1     | 10         | RMEFB      | OP        | sick    | no      | TRUE          |              |     |     |
| SD_1200  |SD_1200 |B5978            | BS2322 | HS53037 | 11-02-0368 | 4/14/2011  | 11/11/2013 | V7         | 31         | NW         | NOP       | healthy | no      | FALSE         |              |     |     |
| SD_1117  |SD_1117 |B6013            | BS2357 | HS53072 | 11-02-0355 | 1/15/2011  | 1/10/2012  | AOM1       | 12         | RMEFB      | OP        | sick    | no      | FALSE         |              |     |     |
| SD_993   |SD_993 |B6017            | BS2361 | HS53076 | 10-02-0325 | 9/25/2010  | 9/26/2013  | V7         | 36         | NWB        | OP        | healthy | no      | FALSE         |              |     |     |
| SD_989   |SD_989 |B6021            | BS2365 | HS53080 | 10-02-0325 | 9/25/2010  | 2/21/2012  | AOM6       | 17         | RMEF       | OP        | sick    | no      | FALSE         |              |     |     |
| SD_537   |SD_537 |B6068            | BS2412 | HS53127 | 8/3/1978   | 8/17/2008  | 5/20/2009  | V2         | 9          | NW         | OP        | healthy | no      | FALSE         |              |     |     |
| SD_1566  |SD_1566 |B6243            | BS2587 | HS53134 | 13-02-512  | 2/8/2013   | 2/21/2014  | V3         | 12         | NW         | NOP       | healthy | no      | FALSE         |              |     |     |
|          |                  |        |         |            |            |            |            |            |            |           |         |         |               |              |     |     |


# NTHi biofilm development time course 
RNA was extracted from diversity panel strains in shaking logarithmic phase, shaking stationary phase, and static early (6hr) biofilms or late supernatant (24hr) and biofilms (24hr; t=1:5; figure). Libraries were prepared using NEBNext Ultra II Directional RNA Library prep kit for Ilumina, and a CRISPR-based ribosomal RNA depletion protocol from JumpCode Genomics was used. RNA sequencing was performed on Illumina NextSeq500 at 1x75 nucleotides targeting 1 million reads per transcriptome.

# pantranscriptome pipeline
## pangenome clustering
Panaroo was used for [bacterial pangenome clustering](https://github.com/gtonkinhill/panaroo). Genes were clustered based on a 90% nucleotide identity (blastP threshold). The 90% strict setting was used as this is the most common (panaroo documents), and is useful for removing potential sources of contamination while retaining important genes. Stricter settings can be adjusted to identify rare plamids; this was not useful for the desired downstream analysis.

The [output from panaroo](https://gtonkinhill.github.io/panaroo/#/gettingstarted/output) is a combined DNA CDS fasta file, which can be used to create each strain's respective transcriptome for salmon as well as GFF needed for annotating genomes for BWA. In the diversity panel, this resulted in a total of 3227 unique genes and their resulting locus ID, or where the gene is found in each diversity panel strain.  
## transcriptome comparisons

### bwa
Each strain's genome sequence was used as its respective reference genome with genes defined by panaroo clusters. Reference genomes were acquired from NCBi or were determined by Pacbio long read sequencing at Drexel University’s Center for Advanced Microbial Processing. Reads were aligned to their reference genome and counts were generated using BEDTOOLS coverage. The resulting count table was imported into RStudio for DEG analysis.

A custom snakemake pipeline was developed for sample processing and alignment.  

### salmon 
Each strain’s transcriptome was used for pseudo aligning by Salmon. Transcriptomes were defined by panaroo clusters and subsetted for an individual strain to reduce off-target aligning. Pseudocounts were imported into RStudio using taxi port for DEG analysis. 

# diverse NTHi have complex population dynamics and phenotypic growth profiles
These 12 strains were characterized in a variety of _in vitro_ phenotypic assays establishing quantitative traits such as growth rates, final planktonic optical density, biofilm biomass (crystal violet stain), response to pH changes, serum resistance, among others and **multi-antibiotic tolerance**. There were no obvious correlations between phenotypes and clinical isolation source, there were notably too few strains for an association study. Strains do vary in phenotypes _in vitro_ in all assays performed.

# 