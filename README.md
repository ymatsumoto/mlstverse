# mlstverse
mlstverse is a general purpose identification software.
Currently database for Mycobacterium species (consists of MTB complex and nontuberculous mycobacteria) is available.

# Installation
## Install R package

```
> git clone ymatsumoto/mlstverse
> R
>>> devtools::install("mlstverse")
```
## Prepare Loci sequences
Downloads [Loci.fasta](https://) to your local strage.
Loci.fasta is also included in mlstverse repository.
https://github.com/ymatsumoto/mlstverse/raw/master/mlstverse_0.1.0_R_x86_64-pc-linux-gnu.tar.gz
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
MLST database used in mlstverse is compatible with other MLST methodology. We are currently preparing to be also available
[pubmlst](https://pubmlst.org) database from mlstverse.

## Citation
Not published yet.

## Licence
Released under MIT + file LICENSE.

## Contact
If you find some problem, please contact matsumoto@gen-info.osaka-u.ac.jp

