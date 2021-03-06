# SameStr
We developed **SameStr** as a bioinformatic tool for the identification of shared microbial strains in metagenomic shotgun sequencing data. SameStr is related to StrainPhlAn as both use the same taxon-specific MetaPhlAn markers to identify and compare species-specific SNV profiles. Difference lies in the handling of SNVs, as StrainPhlAn processes only the majority variant at each position, whereas SameStr considers all possible variants in the alignments. 

SameStr's shared strains are specific to related but not unrelated sample pairs and can therefore be used to track strains across biological samples. As demonstrated with strain co-occurrence networks, this enables further applications such as for the quality screening of mislabelled data and possible contamination, or personal identification which raises further questions regarding study participant privacy. 

SameStr is currently implemented in python 2.7 (should be updated to python 3.X in the future). While the program does not reconstruct conspecific marker sequences, SameStr's outputs (numpy format) can be used for strain composition modelling with probabilistic algorithms.

Preprint: https://www.medrxiv.org/content/10.1101/2020.09.29.20203638v1

# Acknowledgements
We wrapped tools or borrowed code from these repositories:
- [MetaPhlan2](https://bitbucket.org/biobakery/metaphlan2 "MetaPhlAn2 repository") (v2.6.0): Duy Tin Truong, Eric A Franzosa, Timothy L Tickle, Matthias Scholz, George Weingart, Edoardo Pasolli, Adrian Tett, Curtis Huttenhower & Nicola Segata. Nature Methods 12, 902-903 (2015)
- [StrainPhlAn](https://bitbucket.org/biobakery/metaphlan2 "MetaPhlAn2 repository"): Duy Tin Truong, Adrian Tett, Edoardo Pasolli, Curtis Huttenhower, & Nicola Segata. Genome Research 27:626-638 (2017)
- [Strain Finder](https://github.com/cssmillie/StrainFinder "StrainFinder repository"): Strain Tracking Reveals the Determinants of Bacterial Engraftment in the Human Gut Following Fecal Microbiota Transplantation (https://doi.org/10.1016/j.chom.2018.01.003)
- [kpileup](https://github.com/cssmillie/StrainFinder "StrainFinder repository"): (Katherine Huang, Broad Institute)
- [Kneaddata](https://bitbucket.org/biobakery/kneaddata "Kneaddata repository") (v0.6.1)
- [Fastq-Stats](https://expressionanalysis.github.io/ea-utils/ "ea-utils repository") (v1.01): ea-utils

# Overview
**SameStr** identifies shared strains between pairs of metagenomic samples based on the similarity of SNV profiles.
Here, we present a pipeline to process data starting from raw single or paired-end metagenomic shotgun sequencing files to output tables summarizing shared strain calls for individual species and overall strain co-occurrence between sample pairs.

[Installation](#installation)
- [Requirements](#requirements)

**SameStr** must be used from the command line and encompasses multiple modules which can be called by using the following syntax: **`samestr <command>`**. Help for specific command-line usage is available by using the `--help` option with the `samestr` command (`samestr --help`) or any of its modules (`samestr convert --help`).

[Module Description](#description)
- [align](#align): preprocess and align `fastq` files to MetaPhlAn2 markers `sam` 
- [convert](#convert): convert MetaPhlAn alignments `sam` to SNV profiles `npy`
- [extract](#extract): extract SNV profiles `npy` from reference genomes `fasta`
- [merge](#merge): merge SNV profiles `npy` + `npy` from multiple sources
- [filter](#filter): filter SNV profiles `npy`
- [compare](#compare): compare SNV profiles `npy` to determine their similarity and overlapping coverage `tsv`
- [summarize](#summarize): summarize shared strains and strain co-occurrence `tsv`
- [stats](#stats): report coverage stats `tsv` for SNV profiles `npy`
- [db](#db): regenerate species marker db from MetaPhlAn markers `mpa-pkl`

# Installation
Currently, installation is possible with conda by following these 4 steps:

1. Clone this repository, recreate the environment with conda and install SameStr with pip:
```
git clone https://github.com/danielpodlesny/samestr.git
cd samestr
conda env create
conda activate samestr
pip install .
```

2. Get the directory to where pip installed SameStr during the previous command:
```
SAMESTR_INSTALL_DIR=$HOME/.conda/envs/samestr/lib/python2.7/site-packages/samestr
# OR: pip list | grep 'samestr' 
```

3. Add the SameStr installation directory to your path:
```
export PATH=$PATH:${SAMESTR_INSTALL_DIR}/
export PATH=$PATH:${SAMESTR_INSTALL_DIR}/convert/
```

4. Set the following file permissions:
```
chmod +x ${SAMESTR_INSTALL_DIR}/convert/*py
```

## Requirements
SameStr has been tested with the following tool versions:
- [MetaPhlan2](https://bitbucket.org/biobakery/metaphlan2 "MetaPhlAn2 repository") (>=v2.6)
- [Kneaddata](https://bitbucket.org/biobakery/kneaddata "Kneaddata repository") (v0.6.1)
- [Fastq-Stats](https://expressionanalysis.github.io/ea-utils/ "ea-utils repository") (v1.01)
- Samtoools (v0.1.19)
- MUSCLE (v3.8.31)
- blastn (2.8.1+)

# Description
## align
OPTIONAL. This command conveniently wraps `kneaddata`, `fastq-stats`, and `MetaPhlAn2`. Use it to first perform basic quality control measures such as read length trimming and host-genome mapping, gather qc and fastq statistics, and finally map cleaned sequence reads against MetaPhlAn2's species-specific marker sequences.

Although it is not required to follow this exact pre-processing scheme for the **SameStr** analysis, we highly recommend quality processing of your raw sequencing data. When using alternative qc protocols, make sure to align sequences with MetaPhlAn2, further specifying the `-s` option to save the alignments in `sam` format. MetaPhlAn2 `sam` files are required for downstream processing.


### Usage example
The input to **`samestr align`** are `fastq` files with the file extension `fastq` or `fastq.gz`. Use `--input-sequence-type` to specify `single-end` or `paired-end` sequence format. By default, **`samestr align`** expects `paired-end` files and therefore searches for file pairs by considering common read-pair file endings such as `R1.fastq, R2.fastq` or `_1.fastq, _2.fastq`. For parallel processing, specify `--nprocs`.
```
samestr align \
--input-files RAW/*fastq.gz \
--input-sequence-type paired \
--host-bowtie2db Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
--metaphlan2-exe /opt/metaphlan2/metaphlan2.py \
--mpa /opt/metaphlan2/db_v20/mpa_v20_m200 \
--mpa-pkl /opt/metaphlan2/db_v20/mpa_v20_m200.pkl \
--nprocs 30 \
--output-dir out_align/
```

### Output Files
| Tool | Output File Extension | Description |
| :-: | :-: | :- |
| `Kneaddata` | .fastq | Quality processed fastq read file. Excluding reads mapped to the specified host genome. In case of paired-end data, including both orphan and paired read data |
| `Kneaddata` | .log | Protocol of quality processing |
| `Kneaddata` | .qc_stats.txt | Tabulated statistics of quality processing |
| `fastq-stats` | .fastq_summary | Tabulated fastq read statistics |
| `MetaPhlan2` | .bowtie2out | Intermediate bowtie alignment |
| `MetaPhlan2` | .sam.bz2 | Marker sequence alignments |
| `MetaPhlan2` | .profile.txt | Taxonomic assignment & relative abundance table

## samestr convert
Convert MetaPhlAn2 marker alignments to nucleotide variant profiles. 

### Usage example
The input to **`samestr convert`** are MetaPhlAn2 marker alignments with the file extension `.sam` or `.sam.bz2`. For parallel processing, specify `--nprocs`. Note: This step requires a MetaPhlAn database regenerated with [samestr db](#db).
```
samestr convert \
--input-files out_align/metaphlan/*sam.bz2 \
--marker-dir marker_db/ \
--min-vcov 5 \
--nprocs 30 \
--output-dir out_convert/
```
### Output Format
Per Species:
| Nucleotide | Marker 1 (Pos 1) | Marker 1 (Pos n) | Marker 2 (Pos 1) | Marker 2 (Pos n)| Marker n (Pos n) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| A | 0 | 6 | 0 | 0 | .. |
| C | 12 | 4 | 0 | 0 | .. |
| G | 0 | 0 | 0 | 0 | .. |
| T | 0 | 0 | 14 | 0 | .. |

## extract
The module **`samestr extract`** obtains MetaPhlAn2 marker sequences from reference genomes by using BLASTN as described in the StrainPhlAn paper. 

### Usage example
The input to **`samestr extract`** are genomic sequences in `fasta` format. For parallel processing, specify `--nprocs`. Note: This step requires a MetaPhlAn database regenerated with [samestr db](#db).
```
samestr extract \
--input-files reference_genomes/*.fasta \
--marker-dir marker_db/ \
--nprocs 30 \
--output-dir out_extract/
```

## merge
Use **`samestr merge`** to merge nucleotide variant profiles obtained from metagenomic samples (`convert`) and reference genomes (`extract`) for further filtering and pairwise comparisons.

### Usage example
For parallel processing, specify `--nprocs`.
```
samestr merge \
--input-files out_extract/*.npy out_convert/*.npy \
--nprocs 30 \
--output-dir out_merge/
```

## filter
The module **`samestr merge`** provides extensive options for global, per sample, marker, and position filtering of nucleotide variant profiles based on absolute (including n, mean, sd) and relative alignment horizontal and vertical alignment coverage. 

### Usage example
For parallel processing, specify `--nprocs`.
```
samestr filter \
--input-files out_merge/*.npy \
--input-names out_merge/*.names.txt \
--marker-dir marker_db/ \
--samples-min-n-hcov 5000 \
--species-min-samples 2 \
--marker-trunc-len 20 \
--global-pos-min-n-vcov 2 \
--sample-pos-min-n-vcov 5 \
--sample-var-min-f-vcov 0.1 \
--nprocs 30 \
--output-dir out_filter/
```

## compare
Due to the nature of the consensus rule, StrainPhlAn's accuracy is only given when the dominant strain has above 50% relative abundance in the species population (e.g. 70%/20%/10%). Genotypes reconstructed with StrainPhlAn are confounded when the dominant strain does not hold the majority at all positions, such as in more balanced mixtures (e.g. 40%/20%/20%/20%) or in low-coverage regions where subdominant strains can hold the majority (40% < 20%+20%+20%) over the most abundant strain. 

This can lead to false-negatives when assessing strain co-occurrence between samples in which strain mixtures should be expected (such as FMT). By considering all nucleotide variants, SameStr's shared strain calls remain sensitive to more complex strain mixtures without compromising specificity by requiring at least 5kb coverage overlap between alignment pairs and limiting comparisons to variants with ≥10% variant frequency at a given position. 

The module **`samestr compare`** can be used to perform the pairwise species alignment comparison. We recommend using `filter` to compare alignments only at meaningful sites. 

### Usage example
Similarity can be determined by calculating the maximum variant profile similarity `MVS` or consensus similarity `CVS` for all species and sample pairs. For parallel processing, specify `--nprocs`.
```
samestr compare \
--input-files out_filter/*.npy \
--input-names out_filter/*.names.txt \
--nprocs 30 \
--output-dir out_compare/
```

## summarize
Use **`samestr summarize`** to call and summarize shared strains based on defined criteria. This module also generates taxonomic co-occurrence tables at the kingdom, phylum, class, order, family, genus, species, and strain-level.

### Usage example
Specify the minimum required overlap (`aln-pair-min-overlap`) and similarity (`aln-pair-min-similarity`) to define shared strains between nucleotide variant profiles.
```
samestr summarize \
--input-dir out_compare/ \
--mp-profiles-dir out_align/metaphlan/ \
--output-dir out_summarize/
```

## stats
The module **`samestr stats`** can be used to generate statistics related to coverage and nucleotide diversity in nucleotide variant profiles at any step in the analysis pipeline. 

### Usage example
```
samestr stats \
--input-files out_filter/*.npy \
--input-names out_filter/*.names.txt \
--nprocs 30 \
--output-dir out_stats_/
```

## db
The module **`samestr db`** has to be used after installation of SameStr in order to generate database files from MetaPhlAn2 `mpa-pkl` and `all_markers.fasta`. Database files are required for further processing and can be generated for individual species or all MetaPhlAn2 species that are available.

### Usage example
```
samestr db \
--mpa-pkl metaphlan2/db_v20/mpa_v20_m200.pkl \
--mpa-markers metaphlan2/all_markers/all_markers.fasta \
--output-dir marker_db/
```

Note that SameStr is currently implemented in python 2.7 and MetaPhlAn3 has been 
updated to python 3. Due to updates to python's pickle library MetaPhlAn3 `mpa_pkl` files have to
be converted to the old pickle format before they can be used with SameStr:
```
conda create --name py3.7 py=3.7
conda activate py3.7

# within python 3 environment 
import pickle
import bz2
mpa_pkl = 'mpa_v30_CHOCOPhlAn_201901.pkl'
mpa_pkl = pickle.load(bz2.BZ2File(mpa_pkl))

f = bz2.BZ2File(mpa_pkl.replace('.pkl', '.py2.pkl'), 'wb')
pickle.dump(mpa_pkl, f, protocol = 0)
```

 
