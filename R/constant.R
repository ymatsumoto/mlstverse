#loci <- c("BACT000001", "BACT000002", "BACT000003", "BACT000004", "BACT000005",
#          "BACT000006", "BACT000007", "BACT000008", "BACT000009", "BACT000010",
#         "BACT000011", "BACT000012", "BACT000013", "BACT000014", "BACT000015",
#          "BACT000016", "BACT000017", "BACT000018", "BACT000019", "BACT000020",
#          "BACT000021", "BACT000030", "BACT000031", "BACT000032", "BACT000033",
#          "BACT000034", "BACT000035", "BACT000036", "BACT000038", "BACT000039",
#          "BACT000040", "BACT000042", "BACT000043", "BACT000044", "BACT000045",
#          "BACT000046", "BACT000047", "BACT000048", "BACT000049", "BACT000050",
#          "BACT000051", "BACT000052", "BACT000053", "BACT000056", "BACT000057",
#          "BACT000058", "BACT000059", "BACT000060", "BACT000061", "BACT000062",
#          "BACT000063", "BACT000064", "BACT000065")
#len.allele <- c(1526.5, 788.1, 726.1, 615.2, 544.2, 355.3, 471.7, 397.6, 423.3, 311.0,
#                395.8, 386.9, 369.3, 257.5, 270.3, 341.5, 269.3, 250.2, 279.5, 266.1,
#                204.8, 699.3, 831.4, 649.1, 628.7, 550.8, 538.5, 376.0, 467.6, 518.5,
#                430.5, 441.8, 369.3, 449.5, 421.6, 434.0, 364.2, 366.9, 365.0, 336.8,
#                372.7, 303.4, 317.8, 269.6, 229.3, 212.7, 186.2, 237.4, 183.9, 164.0,
#                139.5, 198.2, 118.2)
#names(len.allele) <- loci

#table.loci <- read.csv("data/mlst_loci.csv", header=T, fill=T, stringsAsFactors=F)

#assembly_summary <- asm.sub3

##len.allele -> len.loci


# table.allele <- read.csv("data/mlst_loci.csv", stringsAsFactors = F)
# table.allele$index <- 2000000
# target <- read.csv("data/annotation.csv", stringsAsFactors=F)

# speciesNames <- read.csv("data/MycoStock_BRC_db_20171201.csv", header=F, stringsAsFactors=F)[,1]
# speciesNames <- speciesNames[speciesNames!="Mycobacterium tuberculosis"]


# library(tidyverse)
# table.allele <- read.csv("data/mlst_loci.csv", stringsAsFactors = F)
# table.allele$index <- 2000000
# target <- read.csv("data/annotation.csv", stringsAsFactors=F)

# bigsdb <- rbind(bigsdb, mlst_profiles)
# save(list="bigsdb", file="bigsdb_tmp.RData")
# for (allele in table.allele$mlst) {
#   write.fasta(sequences=mlst_fasta[[allele]],
#               names=paste(allele, names(mlst_fasta[[allele]]), sep="_"),
#               file.out=paste("data/allele_seqs_new/", allele, ".fas", sep=""))
# }


