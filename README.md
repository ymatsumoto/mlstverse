# mlstverse
mlstverse is a general purpose identification software.
Currently database for Mycobacterium species (consists of MTB complex and nontuberculous mycobacteria) is available.

# Installation
## Get source codes 
```
> git clone https://github.com/ymatsumoto/mlstverse
> git clone https://github.com/ymatsumoto/mlstverse.Mycobacterium.db
```
## Install software and dependencies
```
> R
>>> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
>>> BiocManager::install("Rsamtools", version = "3.8")
>>> install.packages("devtools")
>>> devtools::install("mlstverse")
>>> devtools::install("mlstverse.Mycobacterium.db")
```
<!---
## Install mlstverse and database
```
> git clone https://github.com/ymatsumoto/mlstverse
> git clone https://github.com/ymatsumoto/mlstverse.Mycobacterium.db
> R
>>> install.packages("devtools")
>>> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
>>> BiocManager::install("Rsamtools", version = "3.8")
>>> devtools::install("mlstverse")
>>> devtools::install("mlstverse.Mycobacterium.db")
```
--->

## Prepare Loci sequences
Downloads [Loci.fasta](https://github.com/ymatsumoto/mlstverse/raw/master/data/Loci.fasta)
to your local strage. Loci.fasta is also included in mlstverse repository.

# Usage
## Preparing bam Files
Before identifying isolates using R scripts, prepare bam files.
Mapping shotgun sequence reads should be mapped to Loci.fasta using any mapping software and indexed.

```
> bwa mem Loci.fasta input.fastq | samtools sort - -o out.bam
> samtools index out.bam
```

## Identify example (R script)
```R
library(mlstverse)
library(mlstverse.Mycobacterium.db)
filenames <- c("out.bam")
result <- mlstverse(filenames, threads=16)
print(result$score)
```

## Options for mlstverse()
| option | default | description |
|---|---|---|
| filenames | | (required) Character, locations of input bam files
| mlstdb | mlstverse.Mycobacterium.db | (optional) MLST database |
| min_depth | 0 | (optional) Numeric, only use genes larger than the minimum reads depth
| min_ratio | 0.1 | (optional) Numeric, only use genes larger than the ratio to maximum reads depth
| th.pvalue | 0.05 | (optional) Numeric, filter by threshold value in Kolmogorov–Smirnov test
| th.score | 0.1 | (optional) Numeric, filter by threshold value in MLST score
| threads | 1 | (optional) Numeric, number of threads
| normalize | TRUE | (optional) Boolean, if TRUE, normalize coverage
| samfile | NULL | (optional) BamFile object in RSamtools package, input bam file. Excluive with *filenames* option.
| method | "default" | (optional) Character, scoring methods. "default" and "sensitive" are available.

# Other
MLST database used in the mlstverse is compatible with other MLST methodologies. We are currently preparing to be also available
[pubmlst](https://pubmlst.org) database from mlstverse.

## Citation
Not published yet.

## Licence
Released under MIT + file LICENCE.

## Contact
If you find some problem, please contact matsumoto@gen-info.osaka-u.ac.jp

