#
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/danielpodlesny/samestr/main/assets/SameStr-Logo.v1-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/danielpodlesny/samestr/main/assets/SameStr-Logo.v1-light.svg">
  <img src="https://raw.githubusercontent.com/danielpodlesny/samestr/main/assets/SameStr-Logo.v1-light.svg" width="200px"> 
</picture>

#
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![PyPI version](https://img.shields.io/pypi/v/samestr)](https://pypi.python.org/pypi/samestr)
[![Bioconda version](https://anaconda.org/bioconda/samestr/badges/version.svg)](http://bioconda.github.io/recipes/samestr/README.html)
[![Bioconda downloads](https://anaconda.org/bioconda/samestr/badges/downloads.svg)](http://bioconda.github.io/recipes/samestr/README.html)


**SameStr** was developed as a bioinformatic tool for the identification of shared microbial strains in metagenomic shotgun sequencing data. SameStr uses clade-specific markers from [MetaPhlan](https://github.com/biobakery/MetaPhlAn) or [mOTUs](https://github.com/motu-tool/mOTUs) to identify and compare SNV profiles. In this, SameStr is not constrained to comparing the majority variants but considers all possible variants in the alignments.

SameStr's shared strains are specific to related but not unrelated sample pairs and can therefore be used to track strains across biological samples. As demonstrated with strain co-occurrence networks, this enables further applications such as for the quality screening of mislabelled data and possible contamination, or personal identification which raises further questions regarding study participant privacy. 

While the program does not reconstruct conspecific marker sequences, SameStr's outputs (numpy format) can be used for strain composition modelling with probabilistic algorithms.

## Recent Updates:
- 2024.8 added precomputed reference databases to Zenodo, added example data, and removed defaults from the anaconda channels
- 2024.2 fixed some bugs introduced in 2024.1
- 2024.1 added support for [mOTUs](https://github.com/motu-tool/mOTUs) markers (v2, 2.5, 3, +)
- 2024.1 added support for MetaPhlAn SGB-level markers (mpa_vJun23_CHOCOPhlAnSGB_202307, mpa_vOct22_CHOCOPhlAnSGB_202212, mpa_vJun23_CHOCOPhlAnSGB_202307, +)
- 2024.1 improved the database footprint and structure and introduced a database integrity check
- 2024.1 is a major rewrite to support custom databases in the future. Command line options have slightly changed in naming.

SameStr databases generated from [mOTUs](https://github.com/motu-tool/mOTUs) and [MetaPhlan](https://github.com/biobakery/MetaPhlAn) markers are now available on Zenodo in the [Zenodo SameStr Community](https://zenodo.org/communities/samestr):
- SameStr databases: mOTUs [DOI: 10.5281/zenodo.10640396](https://zenodo.org/doi/10.5281/zenodo.10640396)
- SameStr databases: MetaPhlAn [DOI: 10.5281/zenodo.10636302](https://zenodo.org/doi/10.5281/zenodo.10636302)

Instead of generating the database files from scratch with `samestr db`, you can download the database version that you used for alignment, extract it, and start processing your data with `samestr convert` or `samestr extract`.

## Citation
Podlesny, D., Arze, C., Dörner, E., Verma, S., Dutta, S., Walter, J., & Fricke, W. F. (2022). Metagenomic strain detection with SameStr: identification of a persisting core gut microbiota transferable by fecal transplantation. Microbiome, 10(1), 53. https://doi.org/10.1186/s40168-022-01251-w

## Other studies using SameStr
Podlesny, D., Durdevic, M., Paramsothy, S., Kaakoush, N. O., Högenauer, C., Gorkiewicz, G., Walter, J., & Fricke, W. F. (2022). Identification of clinical and ecological determinants of strain engraftment after fecal microbiota transplantation using metagenomics. Cell Reports Medicine, 3(8), 100711. https://doi.org/10.1016/j.xcrm.2022.100711

Podlesny, D., & Fricke, W. F. (2021). Strain inheritance and neonatal gut microbiota development: A meta-analysis. International journal of medical microbiology : IJMM, 311(3), 151483. https://doi.org/10.1016/j.ijmm.2021.151483

# Acknowledgements
SameStr was inspired by [StrainPhlAn](https://github.com/biobakery/MetaPhlAn/wiki/StrainPhlAn-4) in initially using the MetaPhlAn markers to identify and compare species-specific SNV profiles. Difference lies in the handling of SNVs, as StrainPhlAn processes only the majority variant at each position, whereas SameStr considers all possible variants in the alignments.

We wrapped tools or borrowed code from these repositories:
- [MetaPhlan](https://github.com/biobakery/MetaPhlAn "MetaPhlAn repository") (v2.6.0): Duy Tin Truong, Eric A Franzosa, Timothy L Tickle, Matthias Scholz, George Weingart, Edoardo Pasolli, Adrian Tett, Curtis Huttenhower & Nicola Segata. Nature Methods 12, 902-903 (2015)
- [StrainPhlAn](https://github.com/biobakery/MetaPhlAn "StrainPhlAn repository"): Duy Tin Truong, Adrian Tett, Edoardo Pasolli, Curtis Huttenhower, & Nicola Segata. Genome Research 27:626-638 (2017)
- [Strain Finder](https://github.com/cssmillie/StrainFinder "StrainFinder repository"): Strain Tracking Reveals the Determinants of Bacterial Engraftment in the Human Gut Following Fecal Microbiota Transplantation (https://doi.org/10.1016/j.chom.2018.01.003)
- [kpileup](https://github.com/cssmillie/StrainFinder "StrainFinder repository"): (Katherine Huang, Broad Institute)

# Overview
**SameStr** identifies shared strains between pairs of metagenomic samples based on the similarity of SNV profiles.
Here, we present a pipeline to process data starting from raw single or paired-end metagenomic shotgun sequencing files to output tables summarizing shared strain calls for individual clades and overall strain co-occurrence between sample pairs.

[Installation](#installation)
- [Requirements](#requirements)

**SameStr** must be used from the command line and encompasses multiple modules which can be called by using the following syntax: **`samestr <command>`**. Help for specific command-line usage is available by using the `--help` option with the `samestr` command (`samestr --help`) or any of its modules (`samestr convert --help`).

[Module Description](#description)
- [db](#db): regenerate species marker db from MetaPhlAn markers `markers-pkl`
- [convert](#convert): convert MetaPhlAn alignments `sam` to SNV profiles `npy`
- [extract](#extract): extract SNV profiles `npy` from reference genomes `fasta`
- [merge](#merge): merge SNV profiles `npy` + `npy` from multiple sources
- [filter](#filter): filter SNV profiles `npy`
- [stats](#stats): report coverage stats `tsv` for SNV profiles `npy`
- [compare](#compare): compare SNV profiles `npy` to determine their similarity and overlapping coverage `tsv`
- [summarize](#summarize): summarize shared strains and strain co-occurrence `tsv`

# Installation
It is highly recommended to install SameStr in a conda environment, which will install the required dependencies automatically.

## Requirements
SameStr requires python>=3.9 and has been tested with the following software versions:
- Samtoools (v0.1.19)
- MUSCLE (v3.8.31)
- Glibc (v2.28)

SameStr is fully compatible with the databases from mOTUs (v2, 2.5, 3, +) and MetaPhlAn (mpa_vJun23_CHOCOPhlAnSGB_202307, mpa_vOct22_CHOCOPhlAnSGB_202212, mpa_vJun23_CHOCOPhlAnSGB_202307, +).

## bioconda
SameStr can be install through conda by using the following command:
```
conda install -c bioconda -c conda-forge samestr
```

## conda
You can also clone this repository, recreate the environment with conda and install SameStr with pip:
```
git clone https://github.com/danielpodlesny/samestr.git
cd samestr
conda env create -f environment.yml
conda activate samestr
pip install .
```

## pypi
SameStr can be installed from pypi with pip by using the following command:
```
pip install samestr
```

For the pip installation to work, make sure to install the following dependencies manually:

```
blast>=2.6.0
glibc>=2.3
muscle==3.8.1551 and/or mafft==7.515
python>=3.9
samtools==0.1.19
```

# Description
You can refer to the command-line help at any of the following commands:
```
samestr --help
usage: samestr [-h] [--version] [--citation [STR]] [--verbosity STR] {db,convert,extract,merge,filter,stats,compare,summarize} ...

Welcome to SameStr! SameStr identifies shared strains between pairs of metagenomic samples based on the similarity of their Single Nucleotide Variant (SNV) profiles.

optional arguments:
  -h, --help            show this help message and exit
  --version             Show version of the program and exit.
  --citation [STR]      Print citation and exit. Options: Text, BibTex, Endnote, RIS, DOI. (default: None)
  --verbosity STR       Set the verbosity of the program. Options: DEBUG, INFO, WARNING, ERROR, CRITICAL. (default: INFO)

commands:
  {db,convert,extract,merge,filter,stats,compare,summarize}
                        Use one of the following commands for different tasks:
    db                  Make database from MetaPhlAn or mOTUs markers.
    convert             Convert sequence alignments to SNV Profiles.
    extract             Extract SNV Profiles from Reference Genomes.
    merge               Merge SNV Profiles from multiple sources.
    filter              Filter SNV Profiles.
    stats               Report alignment statistics.
    compare             Calculate pairwise sequence similarity.
    summarize           Summarize Taxonomic Co-Occurrence.
```

## db
The module **`samestr db`** has to be used after installation of SameStr in order to generate database files from MetaPhlAn `pickle` file (e.g. `mpa_vJun23_CHOCOPhlAnSGB_202307.pkl`) and `marker` file (e.g. `mpa_vJun23_CHOCOPhlAnSGB_202307.fna.bz2`). Database files are required for further processing and can be generated for individual species or all MetaPhlAn species that are available. 

### Usage example (MetaPhlAn)
Obtain the latest MetaPhlAn database from Segata Lab's FTP server and extract it. As of writing, this is mpa_vJun23_CHOCOPhlAnSGB_202307. 
```
# download
wget -P db_MetaPhlAn/ http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJun23_CHOCOPhlAnSGB_202307.tar

# untar
tar -xvf db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.tar
```

You might want to combine the two separate marker files that have been published with the recent versions of MetaPhlAn to generate a database for all species. This can be done with the following commands:
```
# concatenate
cat db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307_SGB.fna.bz2 \
    db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307_VSG.fna.bz2 > \
        db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.fna.bz2
```
Note that you don't need to decompress the `fna.bz2` file before concatenating or using it with **SameStr**.

Generate a SameStr database from the source files.
```
samestr db \
--markers-info db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.pkl \
--markers-fasta db_MetaPhlAn/mpa_vJun23_CHOCOPhlAnSGB_202307.fna.bz2 \
--db-version db_MetaPhlAn/mpa_latest \
--output-dir samestr_db/
```

### Usage example (mOTUs)
Obtain the latest mOTUs database from the mOTUs Zenodo link and extract it. As of writing, this is mpa_vJun23_CHOCOPhlAnSGB_202307. 

```
# download
wget https://zenodo.org/records/7778108/files/db_mOTU_v3.1.0.tar.gz

# untar
tar -xzvf db_mOTU_v3.1.0.tar.gz
```

Generate a SameStr database from the source files.
```
samestr db \
--markers-info db_mOTU/db_mOTU_taxonomy_ref-mOTUs.tsv db_mOTU/db_mOTU_taxonomy_meta-mOTUs.tsv \
--markers-fasta db_mOTU/db_mOTU_DB_CEN.fasta \
--db-version db_mOTU/db_mOTU_versions \
--output-dir samestr_db/
```

## align
Deprecated. The **`samestr align`** command previously wrapped `kneaddata`, `fastq-stats`, and `MetaPhlAn2` and you could use it to first perform basic quality control measures such as read length trimming and host-genome mapping, gather qc and fastq statistics, and finally map cleaned sequence reads against MetaPhlAn2's species-specific marker sequences.

Please refer to the respective tools to conduct QC preprocessing and alignment. Although it is not required to follow this exact pre-processing scheme for the **SameStr** analysis, we highly recommend quality processing of your raw sequencing data. For example with Kneaddata:

Please refer to the [kneaddata documentation](https://github.com/biobakery/kneaddata) for the latest commands.
```
kneaddata \
    -i ${ID}.R1.fastq.gz \
    -i ${ID}.R2.fastq.gz \
    -db /path-to-db/Homo_sapiens_Bowtie2_v0.1/Homo_sapiens \
    -p 2 \
    -t 15 \
    --max-memory ${RAM} \
    --output-prefix ${ID} \
    --remove-intermediate-output \
    -o QC/ 
```

### MetaPhlAn
When using MetaPhlAn as the source database, make sure to align sequences further specifying the `-s` or `--samout` option to save the alignments in `sam` format. MetaPhlAn `sam` files are required for downstream processing.
```
metaphlan QC/${ID}.R1.fastq.gz,QC/${ID}.R2.fastq.gz \
    --bowtie2db db_MetaPhlAn/ \
    --input_type fastq \
    --nproc 30 \
    -t rel_ab \
    --bowtie2out out_align/${ID}.bowtie2out \
    --samout out_align/${ID}.sam.bz2 \
    -o out_align/${ID}.profile.txt
```

### mOTUs
When using mOTUs as the source database, make sure to generate a mOTUs profile with its `mOTUs profile` command, and further align sequences with `mOTUs map_snv` to generate and save alignments in `bam` format. mOTUs `bam` files are required for downstream processing.
```
# generate profile
motus profile \
    -f QC/${ID}.R1.fastq.gz -r QC/${ID}.R2.fastq.gz \
    -g 1 \
    -t 30 \
    -y insert.raw_counts \
    -o out_align/${ID}.profile.txt 

# generate alignment
motus map_snv \
    -f QC/${ID}.R1.fastq.gz -r QC/${ID}.R2.fastq.gz \  
    -t 30 \
    -o out_align/${ID}.bam 
```

## convert
Convert MetaPhlAn or mOTUs marker alignments to nucleotide variant profiles. 

### Usage example
The input to **`samestr convert`** are marker alignments with the file extension `.sam`, `.sam.bz2`, or `.bam`. For parallel processing, specify `--nprocs`. Note: This step requires a SameStr database regenerated with [samestr db](#db).
```
samestr convert \
--input-files out_align/*sam.bz2 \ 
--marker-dir samestr_db/ \
--nprocs 30 \
--min-vcov 5 \
--output-dir out_convert/
```

### Output Format
This step converts the marker alignment to a nucleotide variant profile. The outputs are per-clade numpy files (`.npy`, `.npz`, `.np.gz` are recognized throughout the script) with the following format:
| Nucleotide | Marker 1 (Pos 1) | Marker 1 (Pos n) | Marker 2 (Pos 1) | Marker 2 (Pos n)| Marker n (Pos n) |
| :---: | :---: | :---: | :---: | :---: | :---: |
| A | 0 | 6 | 0 | 0 | .. |
| C | 12 | 4 | 0 | 0 | .. |
| G | 0 | 0 | 0 | 0 | .. |
| T | 0 | 0 | 14 | 0 | .. |

## extract
The module **`samestr extract`** obtains MetaPhlAn or mOTUs marker sequences from reference genomes by using BLASTN. 

### Usage example
The input to **`samestr extract`** are genomic sequences in `fasta` format. For parallel processing, specify `--nprocs`. Note: This step requires a MetaPhlAn database regenerated with [samestr db](#db).
```
samestr extract \
--input-files reference_genomes/*.fasta \
--marker-dir samestr_db/ \
--nprocs 30 \
--output-dir out_extract/
```

## merge
Use **`samestr merge`** to merge nucleotide variant profiles obtained from metagenomic samples (`convert`) and/or reference genomes (`extract`) for further filtering and pairwise comparisons.

### Usage example
For parallel processing, specify `--nprocs`.
```
samestr merge \
--input-files out_extract/*.npz out_convert/*.npz \
--marker-dir samestr_db_metaphlan/ \
--nprocs 30 \
--output-dir out_merge/
```

## filter
The module **`samestr filter`** provides extensive options for global, per sample, marker, and position filtering of nucleotide variant profiles based on absolute (including n, mean, sd) and relative alignment horizontal and vertical alignment coverage. Refer to the command line `--help` for more information.

### Usage example
For parallel processing, specify `--nprocs`.
```
samestr filter \
--input-files out_merge/*.npy \
--input-names out_merge/*.names.txt \
--marker-dir samestr_db/ \
--clade-min-n-hcov 5000 \
--clade-min-samples 2 \
--marker-trunc-len 20 \
--global-pos-min-n-vcov 2 \
--sample-pos-min-n-vcov 5 \
--sample-var-min-f-vcov 0.1 \
--nprocs 30 \
--output-dir out_filter/
```

## stats
The module **`samestr stats`** can be used to generate statistics related to coverage and nucleotide diversity in nucleotide variant profiles at any step in the analysis pipeline. 

### Usage example
```
samestr stats \
--input-files out_filter/*.npy \
--input-names out_filter/*.names.txt \
--marker-dir samestr_db/ \
--nprocs 30 \
--output-dir out_stats_/
```

## compare
The module **`samestr compare`** can be used to perform the pairwise clade-specific alignment comparison. We recommend using `filter` to compare alignments only at meaningful sites. 

### Usage example
Similarity can be determined by calculating the maximum variant profile similarity `MVS` or consensus similarity `CVS` for all species and sample pairs. For parallel processing, specify `--nprocs`.
```
samestr compare \
--input-files out_filter/*.npy \
--input-names out_filter/*.names.txt \
--marker-dir samestr_db/ \
--nprocs 30 \
--output-dir out_compare/
```

Difference to StrainPhlAn when using MetaPhlAn as the source database: Due to the nature of the consensus rule, StrainPhlAn's accuracy is only given when the dominant strain has above 50% relative abundance in the species population (e.g. 70%/20%/10%). Genotypes reconstructed with StrainPhlAn are confounded when the dominant strain does not hold the majority at all positions, such as in more balanced mixtures (e.g. 40%/20%/20%/20%) or in low-coverage regions where subdominant strains can hold the majority (40% < 20%+20%+20%) over the most abundant strain. 

This can lead to false-negatives when assessing strain co-occurrence between samples in which strain mixtures should be expected (such as FMT). By considering all nucleotide variants, SameStr's shared strain calls remain sensitive to more complex strain mixtures without compromising specificity by requiring at least 5kb coverage overlap between alignment pairs and limiting comparisons to variants with ≥10% variant frequency at a given position. 

## summarize
Use **`samestr summarize`** to call and summarize shared strains based on defined criteria. This module also generates taxonomic co-occurrence tables at the kingdom, phylum, class, order, family, genus, species, and strain-level.

### Usage example
Specify the minimum required overlap (`aln-pair-min-overlap`) and similarity (`aln-pair-min-similarity`) to define shared strains between nucleotide variant profiles.
```
samestr summarize \
--input-dir out_compare/ \
--tax-profiles-dir out_align/ \
--marker-dir samestr_db/
--output-dir out_summarize/
```
