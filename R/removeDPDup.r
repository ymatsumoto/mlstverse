#library(dplyr)
#library(seqinr)
#library(snowfall)

#genes <- gsub(".dup", "", dir("~/mlstverse/data/locus_seqs_merged", pattern=".dup"))

filter.j <- function(x) {
  g.dup <- x[!is.na(x)][-(1:2)]
  which(gids %in% g.dup)
}
updateGIDs <- function(g, gids) {
  tmp <- read.table(paste("~/mlstverse/data/locus_seqs_2/", g, ".dup", sep=""), header=F, fill=T, stringsAsFactors=F)
  tmp <- apply(tmp, 2, function(x) {gsub("Locus_[0-9]+_([0-9]+),*", "\\1", x)})
  tmp <- apply(tmp, c(1,2), as.numeric)

  j <- apply(tmp, 1, filter.j)
  for (i in 1:nrow(tmp)) {
    gids[j[[i]]] <- tmp[i,2]
  }
  gids
}

updateGID <- function(x, locus) {
  tmp <- c()
  for (i in 1:length(x)) {
    k <- which(unlist(lapply(duplicated, function(d) {any(d%in%x[i])})))
    if (length(k) > 0) {
      tmp <- c(tmp, duplicated[[k]][1])
    } else {
      tmp <- c(tmp, x[i])
    }
  }
  return(sort(tmp))
}

#mlstdb <- mlst_profiles
#mlstdb.my <- mlst.my$profiles
#sfInit(parallel=TRUE, cpus=16)
#for (locus in loci) {
#  if (locus=="Locus_0021") {
#    next
#  }
#  dupfile <- paste("~/mlstverse/data/locus_seqs_merged/", locus, ".dup", sep="")
#  if(!file.exists(dupfile)) {
#    next
#  }
#  tmp <- scan(dupfile, what="character", sep="\t")
#  tmp <- strsplit(tmp[seq(length(tmp))%%2==0], ", ")
#  duplicated <- lapply(tmp, function(x) {as.numeric(gsub(paste(locus,"_",sep=""), "", x))})
#  sfExport("duplicated")
#  mlstdb.my[[locus]] <- sfLapply(mlstdb.my[[locus]], updateGID, locus)
#}
#mlstdb <- bind_rows(mlstdb, mlstdb.my)


#for (i in seq(mlstdb[[locus]])) {
#  for (j in seq(mlstdb[[locus]][[i]])) {
#    tmp <- c()
#    k <- which(unlist(lapply(duplicated, function(d) {any(d%in%j)})))
#
#    #x <- unlist(lapply(duplicated, function(d) {if(j%in%d) {return(d[1])}}))
#    if (length(k) > 0) {
#      tmp <- c(tmp, duplicated[[k]][1])
#    } else {
#      tmp <- c(tmp, mlstdb[[locus]][[i]][j])
#    }
#  }
#  mlstdb[[locus]][[i]] <- tmp
#}

#   mlstdb[[locus]][[1]] <- c(1)
#   mlstdb[[locus]]
#   gids <- mlstdb[[locus]]
#   tmp <- read.table(paste("~/mlstverse/data/locus_seqs_2/", locus, ".dup", sep=""), header=F, fill=T, stringsAsFactors=F)
#   for (duplicated in tmp) {
#     lapply()
#   }
#   function(gid) {
#
#   }
#
#   tmp <- gsub(paste(locus,"_",sep=""), "", tmp)
#   tmp[seq(length(tmp))%%2==0][1]
#   tmp <- apply(tmp, 2, function(x) {gsub("Locus_[0-9]+_([0-9]+),*", "\\1", x)})
#   tmp <- apply(tmp, c(1,2), as.numeric)
#
#   cat(paste("Preparing for ", g, ".\n", sep=""))
#   sfExport(list=list("tmp", "g", "gids"))
#
#   j <- sfApply(tmp, 1, filter.j)
#
#   for (i in 1:nrow(tmp)) {
#     gids[j[[i]]] <- tmp[i,2]
#   }
#   bigsdb[,g] <- gids
# }
#
# sfExport(list=list("tmp", "g", "gids"))
# gids <- bigsdb[, g]
# updateGIDs(g, gids)
#
# j <- apply(tmp, 1, filter.j)
# for (i in 1:nrow(tmp)) {
#   gids[j[[i]]] <- tmp[i,2]
# }
#
# for (g in genes) {
#   gids[j[[i]]] <- tmp[i,2]
# }
# bigsdb[,g] <- gids
#
# sfStop()
#
# bigsdb.sub <- distinct(bigsdb, .keep_all=T,
#                        BACT000001,BACT000002,BACT000003,BACT000004,BACT000005,BACT000006,
#                        BACT000007,BACT000008,BACT000009,BACT000010,BACT000011,BACT000012,
#                        BACT000013,BACT000014,BACT000015,BACT000016,BACT000017,BACT000018,
#                        BACT000019,BACT000020,BACT000021,BACT000030,BACT000031,BACT000032,
#                        BACT000033,BACT000034,BACT000035,BACT000036,BACT000038,BACT000039,
#                        BACT000040,BACT000042,BACT000043,BACT000044,BACT000045,BACT000046,
#                        BACT000047,BACT000048,BACT000049,BACT000050,BACT000051,BACT000052,
#                        BACT000053,BACT000056,BACT000057,BACT000058,BACT000059,BACT000060,
#                        BACT000061,BACT000062,BACT000063,BACT000064,BACT000065,species)
#
# i <- which(bigsdb.sub$species%in%c("Mycobacterium indicus", "Mycobacterium indicus pranii "))
# bigsdb.sub$species[i] <- "Mycobacterium indicus pranii"
#
# #length(unique(subset(bigsdb.sub, grepl("Mycobacterium", species))$species))
# species <- sort(unique(subset(mlstdb, grepl("Mycobacterium", species))$species))
# species.MTBC <- c("Mycobacterium tuberculosis",
#                   "Mycobacterium africanum",
#                   "Mycobacterium orygis",
#                   "Mycobacterium bovis",
#                   "Mycobacterium microti",
#                   "Mycobacterium canetti",
#                   "Mycobacterium caprae",
#                   "Mycobacterium pinnipedii",
#                   "Mycobacterium suricattae",
#                   "Mycobacterium mungi")
# species.leprae <- c("Mycobacterium leprae",
#                     "Mycobacterium lepromatosis")
# species.NTM <- species[!species %in% c(species.MTBC, species.leprae)]
#
#
# #length(unique(subset(bigsdb.sub, grepl("Mycobacterium", species))$species))
# species.NTM[!grepl("subsp.", species.NTM)]
# nrow(subset(bigsdb.sub, species %in% species.MTBC))
# nrow(subset(bigsdb.sub, species %in% species.NTM))
#
# save(list=c("bigsdb.sub"), file="~/mlstverse/bigsdb.RData")
# save(list=c("bigsdb"), file="~/mlstverse/bigsdb_verbose.RData")
#
