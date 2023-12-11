###########################################################################################################
2023
leanie
omics 
nthi diversity panel
###########################################################################################################
# establishing diversity strain panel --> not omics
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


# NTHi biofilm development time course --> not omics
RNA was extracted from diversity panel strains in shaking logarithmic phase, shaking stationary phase, and static early (6hr) biofilms or late supernatant (24hr) and biofilms (24hr; t=1:5; figure). Libraries were prepared using NEBNext Ultra II Directional RNA Library prep kit for Ilumina, and a CRISPR-based ribosomal RNA depletion protocol from JumpCode Genomics was used. RNA sequencing was performed on Illumina NextSeq500 at 1x75 nucleotides targeting 1 million reads per transcriptome.

![divPanel expDesign](/home/jupyter-purplekestrel/repos/OMICS/projects/purpleKestrel/divPanel_expDesign.png)

# pantranscriptome pipeline --> OMICS
## pangenome clustering
Panaroo was used for [bacterial pangenome clustering](https://github.com/gtonkinhill/panaroo). Genes were clustered based on a 90% nucleotide identity (blastP threshold). The 90% strict setting was used as this is the most common (panaroo documents), and is useful for removing potential sources of contamination while retaining important genes. Stricter settings can be adjusted to identify rare plamids; this was not useful for the desired downstream analysis.

The [output from panaroo](https://gtonkinhill.github.io/panaroo/#/gettingstarted/output) is a combined DNA CDS fasta file, which can be used to create each strain's respective transcriptome for salmon as well as GFF needed for annotating genomes for BWA. In the diversity panel, this resulted in a total of 3227 unique genes and their resulting locus ID, or where the gene is found in each diversity panel strain.  
## transcriptome comparisons

### bwa 
Each strain's genome sequence was used as its respective reference genome with genes defined by panaroo clusters. Reference genomes were acquired from NCBi or were determined by Pacbio long read sequencing at Drexel University’s Center for Advanced Microbial Processing. Reads were aligned to their reference genome and counts were generated using BEDTOOLS coverage. The resulting count table was imported into RStudio for DEG analysis.
**pan genome n= 3227 genes**
**core genome n=1263 genes**

| strain   | BWA genes |
| -------- | --------- |
| RdKW20   | 1621      |
| 86-028Np | 1791      |
| R2866    | 1812      |
| SD537    | 1802      |
| SD989    | 1812      |
| SD993    | 1867      |
| SD1117   | 1908      |
| SD1200   | 1768      |
| SD1430   | 1880      |
| SD1474   | 1840      |
| SD1566   | 1711      |
| SD1763   | 1733      |

A custom snakemake pipeline was developed for sample processing and alignment.  
![BWA snakemake](/home/jupyter-purplekestrel/repos/OMICS/projects/purpleKestrel/snakemake_workflow.png)

#### snakemake pipeline --> OMICS!

On the cluster, found at `/projects/snakemake_hej14/`
BWA snakemake pipeline which aligns a single RNAseq experiment (paired-end reads) to a single genome.
Snakefile:
``` python
import pandas as pd
from os import path

configfile: "config.yaml"
source = "/home/leanie/projects/snakemake_hej14"
SAMPLES = pd.read_csv(config["samples"], sep=',', header=0,index_col=0)
sample_ids = SAMPLES.index.tolist()
exp = config["experiment"]
genome = config["genome"]

rule all:
        input:
                expand("{source}/sorted_reads/{sample}.bam", source=source, sample=sample_ids),
                expand("{source}/results/{exp}.samtools_flagstat.tsv", source=source, exp=exp),
                expand("{source}/results/{sample}.counts", source=source, sample=sample_ids),
                expand("{source}/results/{exp}.counts.csv", source=source, exp=exp),
                expand("{source}/results/{exp}.map_total_ratio_hist.png", source=source, exp=exp)

rule bwa_index:
        input:
                gnome = config["genome"],
        output:
                idx = multiext(str(genome), ".amb", ".ann", ".bwt", ".pac", ".sa"),
        shell:
                "bwa index {input.gnome} > {output.idx}"

rule bwa_map:
        input:
                indexed_genome = rules.bwa_index.output.idx,
                genome = config["genome"],
                fq1 = lambda wildcards: SAMPLES.loc[wildcards.sample, 'Fastq1'],
                fq2 = lambda wildcards: SAMPLES.loc[wildcards.sample, 'Fastq2'],
        output:
                bam = "{source}/sorted_reads/{sample}.bam",
        params:
                rg  = lambda wildcards: "@RG\\tID:"+wildcards.sample+"\\tSM:"+wildcards.sample,
        threads: config["threads"]
        shell:
                "bwa mem -M -t {threads} -R '{params.rg}' {input.genome} {input.fq1} {input.fq2}| "
                "samtools view -hb - | "
                "samtools sort -m 3G - > {output.bam}"
rule samtools_flagstat:
        input:
                alndir = "{source}/sorted_reads"
        output:
                flagstat = "{source}/results/{exp}.samtools_flagstat.tsv",
                hist = "{source}/results/{exp}.map_total_ratio_hist.png",
        shell:
                "python3 scripts/getSamtoolsFlagstat.py --input_dir {input.alndir} --output_file {output.flagstat} --output_hist {output.hist}"

rule col_headers:
        output:
                sampCounts = "{source}/results/{sample}.temp.counts",
        shell:
                "echo -e 'chromosome\tstart\tstop\tgeneName\tscore\tstrand\tcounts\t' > {output.sampCounts}"
rule bedtools_coverage:
        input:
                sampCounts = rules.col_headers.output.sampCounts,
                map = "{source}/sorted_reads/{sample}.bam",
                bed = config["bed"],
        output:
                covCounts = "{source}/results/{sample}.counts",
        shell:
                "/home/leanie/miniconda3/bin/bedtools coverage -s -a {input.map} -b {input.bed} -counts >> {output.covCounts}"

rule make_count_table:
        input:
                resultsdir = "{source}/results"
        output:
                "{source}/results/{exp}.counts.csv",
        shell:
                "python3 scripts/combineCount.py --res_dir {input.resultsdir} --output_file {output}"
```



### salmon 
Each strain’s transcriptome was used for pseudo aligning by Salmon. Transcriptomes were defined by panaroo clusters and subsetted for an individual strain to reduce off-target aligning. Pseudocounts were imported into RStudio using taxi port for DEG analysis. 

![BWA vs Salmon DEGs](/home/jupyter-purplekestrel/repos/OMICS/projects/purpleKestrel/totalDEGS_Picture1.png)

| strain   | salmon genes |
| -------- | --------- |
| RdKW20   | 1633      |
| 86-028Np | 1808      |
| R2866    | 1820      |
| SD537    | 1806      |
| SD989    | 1830      |
| SD993    | 1878      |
| SD1117   | 1924      |
| SD1200   | 1782      |
| SD1430   | 1895      |
| SD1474   | 1844      |
| SD1566   | 1722      |
| SD1763   | 1737      |


# "discussion" and final thoughts

Snakemake workflow managment is definitely a difficult tool to grasp, however it is likely what I will need to perform any additonal pangenome and pantranscriptome analyses. Once I can successfully integrate this workflow with the custom scripts required for multiple strains, I will be able to do many of the analyses required for my aim1. 

Between alignment methods, raw gene counts (n=1) were comparable. 

| sample                    | gene | BWA counts | SAL counts |
| ------------------------- | ---- | ---------- | ---------- |
| biofilm.RdKW20.A.S49      | relA | 970        | 918        |
| biofilm.RdKw20B.S109      | relA | 4282       | 4125       |
| biofilm.RdKw20.C.S169     | relA | 1439       | 1385       |
| early.RdKW20.A.S25        | relA | 720        | 701        |
| early.RdKw20.B.S85        | relA | 1274       | 1203       |
| early.RdKw20.C.S145       | relA | 1673       | 1577       |
| exponential.RdKW20.A.S1   | relA | 3634       | 3438       |
| exponential.RdKW20.B.S61  | relA | 3703       | 3519       |
| exponential.RdKw20.C.S121 | relA | 1590       | 1502       |
| late.RdKW20.B.S97         | relA | 6334       | 6059       |
| late.RdKW29.C.S157        | relA | 1840       | 1757       |
| stationary.RdKW20.A.S14   | relA | 534        | 509        |
| stationary.RdKW20.B.S73   | relA | 1405       | 1345       |
|                           |      |            |            |
**samples which meet raw count threshold > 1 million reads / transcriptome**

![relA geneCounts BWA vs Salmon](/home/jupyter--purplekestrel/repos/OMICS/projects/purpleKestrel/relA_geneCounts_bwaVSsal_compare.png)


