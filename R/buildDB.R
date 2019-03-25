table.loci <- read.csv("data/mlst_loci.csv", header=T, fill=T, stringsAsFactors=F)

#' Title
#'
#' @return
#' @export
#'
#' @examples
getAssemblySummary <- function() {
  tmp <- tempfile()
  url <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
  download.file(url,tmp, quiet=T)
  return(read.table(tmp, sep="\t", header=T, stringsAsFactors=F, skip=1, fill=T, comment.char="", quote=""))
}

download.fasta <- function(url, fasta_path) {
  while (1) {
    e <- try(download.file(url, fasta_path, quiet=T), silent=F)
    if (class(e) == "try-error") {
      Sys.sleep(30)
    } else {
      break
    }
  }
}

# importGFF
mlstDB.importGFF <- function(gff_path) {
  gff <- read.table(gff_path, fill=T, stringsAsFactors=F, comment.char="", quote="", skip=1, sep="\t")
  colnames(gff) <- c("seqname", "source", "feature",
                     "start", "end", "score",
                     "strand", "frame", "attribute")
  return(gff)
}


#' Title
#'
#' @param ftp_path
#' @param fasta_dir
#'
#' @return
#' @export
#'
#' @examples
getRefSeqFasta <- function(ftp_path, fasta_dir="temporary") {
  tmp <- strsplit(ftp_path, "/")[[1]]
  refseqID <- tmp[length(tmp)]
  filename <- paste(refseqID, "_cds_from_genomic.fna.gz", sep="")
  if (fasta_dir == "temporary") {
    fasta_path <- tempfile()
    url <- paste(ftp_path, "/", filename, sep="")
    download.fasta(url, fasta_path)
  } else {
    fasta_path <- paste(fasta_dir, "/", filename, sep="")
    if (!file.exists(fasta_path)|file.size(fasta_path)==0) {
      url <- paste(ftp_path, "/", filename, sep="")
      download.fasta(url, fasta_path)
    }
  }
  return(seqinr::read.fasta(file=fasta_path))
}
getRefSeqCDS <- getRefSeqFasta
getRefSeqRNA <- function(ftp_path, fasta_dir="temporary") {
  tmp <- strsplit(ftp_path, "/")[[1]]
  refseqID <- tmp[length(tmp)]
  filename <- paste(refseqID, "_rna_from_genomic.fna.gz", sep="")
  if (fasta_dir == "temporary") {
    fasta_path <- tempfile()
    url <- paste(ftp_path, "/", filename, sep="")
    download.fasta(url, fasta_path)
  } else {
    fasta_path <- paste(fasta_dir, "/", filename, sep="")
    if (!file.exists(fasta_path)|file.size(fasta_path)==0) {
      url <- paste(ftp_path, "/", filename, sep="")
      download.fasta(url, fasta_path)
    }
  }
  return(seqinr::read.fasta(file=fasta_path))
}


#' Title
#'
#' @param fasta
#'
#' @return
#' @export
#'
#' @examples
extractRPs.NCBI <- function(fasta) {
  annot <- unlist(lapply(fasta, attr, "Annot"))
  i <- grep("ribosomal protein", lapply(fasta, attr, "Annot"))
  pnames <- annot[i]
  pnames <- strsplit(pnames, "\\] \\[")
  pnames <- unlist(lapply(pnames,
                          function(pname) {
                            sub("protein=",
                                "",
                                pname[grep("protein=", pname)])}))
  pnames <- gsub(".*([53]0S ribosomal protein [SL][/L0-9]+).*", "\\1", pnames)

  rp <- fasta[i]
  names(rp) <- pnames
  return(rp)
}

extractRPs.RAST <- function(fasta, gff) {
  gff.rp <- subset(gff, grepl("ribosomal protein", attribute))
  rp.ids <- gsub("ID=(.*);.*", "\\1", gff.rp$attribute)
  list.rp <- lapply(fasta[rp.ids], c2s)

  rp.names <- gsub("p", "", gsub(".*([SL][0-9]+p|L7/L12).*", "\\1", gff.rp$attribute))
  suType <- substr(rp.names, 1, 1)
  suSize <- unlist(lapply(suType, switch, L="50S", S="30S"))
  names(list.rp) <- paste(suSize, "ribosomal protein", rp.names)
  return(list.rp)
}

extractrRNA.RAST <- function(fasta, gff) {
  gff.rrna <- subset(gff, grepl("(Large|Small) Subunit Ribosomal RNA", attribute))
  ids.rrna <- gsub("ID=(.*);.*", "\\1", gff.rrna$attribute)
  list.rrna <- lapply(fasta[ids.rrna], c2s)

  tmp <- gsub(".*([SL]SU rRNA).*", "\\1", gff.rrna$attribute)
  names.rrna <- sapply(tmp, switch,
                       "SSU rRNA"="16S ribosomal RNA",
                       "LSU rRNA"="23S ribosomal RNA")

  names(list.rrna) <- names.rrna
  return(list.rrna)
}

extractOrthologs <- function(fasta,
                             table.loci,
                             blast_path,
                             th_pident=80,
                             th_ratioMaxScore=0.8,
                             isRASToutput=FALSE) {
  blast <- read.table(blast_path, stringsAsFactors=F)
  colnames(blast) <- c("qseqid", "sseqid", "pident", "length", "mismatch",
                       "gapopen", "qstart", "qend", "sstart", "send",
                       "evalue", "bitscore", "qlen", "slen")
  #genelist <- "~/mlstverse/data/H37Rv/H37Rv.txt"
  #genenames <- read.table(genelist, sep="\t", header=T, fill=T, stringsAsFactors=F)
  #seqs <- vector("list", nrow(genenames))
  tmp <- subset(table.loci, method=="blast")
  seqs <- vector("list", nrow(tmp))
  names(seqs) <- tmp$mlst_tag
  for (i in 1:nrow(tmp)) {
    target <- tmp$blast_tag[i]
    maxScore <- max(0, subset(blast, grepl(target,sseqid))$bitscore)
    matched <- subset(blast, grepl(target, sseqid) &
                        pident > th_pident &
                        th_ratioMaxScore*maxScore < bitscore)
    if (isRASToutput) {
      q <- gsub("fig\\|", "", matched$qseqid)
    } else {
      q <- matched$qseqid
    }
    j <- unique(unlist(lapply(q, grep, names(fasta))))
    if (length(j) > 0) {
      seqs[[i]] <- lapply(fasta[j], seqinr::c2s)
#      names(seqs[[i]]) <- paste(names(seqs[i]), sprintf("%06d", 1:length(j)), sep="_")
    }
  }
  seqs
}

#lapply(x, function(i){lapply(fasta[i], c2s)})
#  products <- subset(table.loci, method == "product" & source == "cds")$product
extractProduct <- function(fasta, table.loci, element="protein") {
  annot <- unlist(lapply(fasta, attr, "Annot"))
  products <- paste("\\[", element, "=", table.loci$product, "\\]", sep="")
  matched <- lapply(products, grep, annot)
  seqs <- lapply(matched, function(i){lapply(fasta[i], seqinr::c2s)})
  names(seqs) <- table.loci$mlst_tag
  return(seqs)
}

createProfile <- function(assembly_summary,
                          fasta_dir="~/mlstverse/data/CDS",
                          blast_dir="~/mlstverse/data/blast") {
  ftp_path <- assembly_summary["ftp_path"]
  blast_path <- paste(blast_dir, "/",
                      gsub(".*/", "", ftp_path), ".blast", sep="")
  print(blast_path)
  cds <- getRefSeqCDS(ftp_path, fasta_dir=paste(fasta_dir, "/cds", sep=""))
  rna <- getRefSeqRNA(ftp_path, fasta_dir=paste(fasta_dir, "/rna", sep=""))
  cds_seqs <- extractProduct(cds, subset(table.loci, method == "product" & source == "cds"), element="protein")
  rna_seqs <- extractProduct(rna, subset(table.loci, method == "product" & source == "rna"), element="product")
  ortho_seqs <- extractOrthologs(cds, subset(table.loci, method=="blast"), blast_path)

  mlst_seqs <- c(cds_seqs, rna_seqs, ortho_seqs)
  # Add MLST profiles
  m <- sapply(mlst_seqs, length)
  profile <- tibble::as.tibble(lapply(m, function(i) {if(i==0) {list(-1:-1)} else {list(1:i)}}))
  notes <- paste(assembly_summary["X..assembly_accession"],
                 assembly_summary["organism_name"],
                 assembly_summary["infraspecific_name"],
                 assembly_summary["isolate"],
                 assembly_summary["assembly_level"],
                 assembly_summary["seq_rel_date"], sep=",")
  mlst_profiles <-
    dplyr::bind_cols(
      dplyr::data_frame(
        rST=1,
        genus="",
        species=assembly_summary["species_name"],
        subspecies="",
        lineage="",
        sublineage="",
        other_designation="",
        notes=notes),
      profile)
  return(list(seqs=mlst_seqs, profiles=mlst_profiles))
}

createCustomRASTProfile <-
  function(fileprefix,
           gff_dir  ="/imetgpfs/projects/cw/OKNW/n039_to_n049/analysis/05_RAST/01_gff",
           fasta_dir="/imetgpfs/projects/cw/OKNW/n039_to_n049/analysis/05_RAST/02_fna",
           blast_dir="~/mlstverse/data/blast_newlySequenced",
           genus="",
           species_name="",
           strain_name="",
           date="",
           assembly_level="") {
  gff_path <- paste(gff_dir, "/", fileprefix, ".gff", sep="")
  gff <- mlstDB.importGFF(gff_path)
  fasta_path <- paste(fasta_dir, "/", fileprefix, ".fna", sep="")
  tmp <- seqinr::read.fasta(fasta_path)
  names(tmp) <- gsub("ID=(.*?);.*", "\\1", gff$attribute)
  gff <- dplyr::distinct(gff)
  fasta <- tmp[gsub("ID=(.*?);.*", "\\1", gff$attribute)]

  rp <- extractRPs.RAST(fasta, gff)
  rrna <- extractrRNA.RAST(fasta, gff)
  x <- c(rp, rrna)
  r_seqs <- list()
  for (i in which(table.loci$method=="product")) {
    j <- table.loci$product[i] == names(x)
    if (any(j)) {
      r_seqs[[table.loci$mlst_tag[i]]] <- x[j]
    } else {
      r_seqs[[table.loci$mlst_tag[i]]] <- NULL
    }
  }

  blast_path <- paste(blast_dir, "/", fileprefix, ".blast", sep="")
  ortho_seqs <- extractOrthologs(fasta, subset(table.loci, method=="blast"), blast_path, isRASToutput=TRUE)

  mlst_seqs <- c(r_seqs, ortho_seqs)
  # Add MLST profiles
  m <- sapply(mlst_seqs, length)
  profile <- tibble::as.tibble(lapply(m, function(i) {if(i==0) {list(-1:-1)} else {list(1:i)}}))
  notes <- paste("",
                 species_name,
                 "",
                 strain_name,
                 assembly_level,
                 date, sep=",")
  mlst_profiles <-
    dplyr::bind_cols(
      dplyr::data_frame(
        rST=1,
        genus=genus,
        species=species_name,
        subspecies="",
        lineage="",
        sublineage="",
        other_designation="",
        notes=notes),
      profile)
  return(list(seqs=mlst_seqs, profiles=mlst_profiles))
}

#speciesList <- read.csv("~/mlstverse/data/locus_seqs_newlySequenced/species_list.csv",
#                        header=F, stringsAsFactors=F)
#tmp <- list()
#for (i in 1:nrow(speciesList)) {
#  cat(paste("Proccessing ", speciesList[i,1], "\n"))
#  tmp[[speciesList[i,1]]] <- list()
#  tmp[[speciesList[i,1]]][[speciesList[i,2]]] <- createCustomRASTProfile(speciesList[i,2],
#                                                          species_name=speciesList[i,1])
#}
#mlst.my <- mlstDB.resetIDs(tmp, table.loci, rep(1e6, length(loci)), 2e6)
#for (locus in table.loci$mlst_tag) {
#  if (length(mlst.my$seqs[[locus]]) > 0) {
#    write.fasta(sequences=mlst.my$seqs[[locus]],
#                names=paste(locus, names(mlst.my$seqs[[locus]]), sep="_"),
#                file.out=paste("~/mlstverse/data/locus_seqs_newlySequenced/", locus, ".fasta", sep=""))
#  }
#}

mlstDB.subset <- function(speciesName, assembly_summary) {
  target <- assembly_summary$organism_name
  if (grepl("subsp.", speciesName)) {
    i <- which(grepl(speciesName, target))
  } else {
    i <- which(grepl(speciesName, target) &
                 !grepl("subsp.", target))
  }
  return(i)
}

#' Title
#'
#' @param speciesNames
#' @param assembly_summary
#'
#' @return
#' @export
#'
#' @examples
mlstDB.subsetASM <- function(speciesNames, assembly_summary) {
  i <- unlist(lapply(speciesNames, mlstDB.subset, assembly_summary))
  dup <- which(duplicated(i))
  if (length(dup) > 0) {
    print(paste("ambiguous species names were detected:",
                assembly_summary$organism_name[dup]))
    return()
  }
  return(assembly_summary[i,])
}

#' Title
#'
#' @param i
#' @param strains
#' @param speciesName
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom tidyr %>%
buildProfiles <- function(assembly_summary, species, fasta_dir="~/mlstverse/data/refseq") {
  mlst_seqs <- array(list(list()), nrow(table.loci))
  names(mlst_seqs) <- table.loci$mlst
  ftp_path <- assembly_summary["ftp_path"]
  fasta <- getRefSeqFasta(ftp_path, fasta_dir=fasta_dir)
  rp <- extractRPs.NCBI(fasta)

  # Add MLST fasta
  if(length(rp)==0) {
    return()
  }
  for (j in seq(table.loci$mlst)) {
    locus <- table.loci$mlst[j]
    product <- table.loci$product[j]
    k <- which(table.loci$product[j]==names(rp))
    if (length(k) == 0) {
      next
    }
    mlst_seqs[[locus]] <- lapply(rp[k], as.character)
    names(mlst_seqs[[locus]]) <- seq(mlst_seqs[[locus]])
  }

  # Add MLST profiles
  m <- max(sapply(mlst_seqs, length))
  profile <- data.frame(lapply(lapply(mlst_seqs, names), function(x) {
    if (is.null(x)) {
      return(as.numeric(rep_len(-1, m)))
    } else {
      return(as.numeric(rep_len(x, m)))
    }
  }))
  profile <- as.data.frame(eval(parse(text=paste("tidyr::expand(profile,", paste(table.loci$mlst, collapse=", "), ")"))))
  notes <- paste(assembly_summary["X..assembly_accession"],
                 assembly_summary["organism_name"],
                 assembly_summary["infraspecific_name"],
                 assembly_summary["assembly_level"],
                 assembly_summary["seq_rel_date"], sep=",")
  mlst_profiles <- data.frame(rST=1:nrow(profile),
                              profile,
                              genus="",
                              species=species,
                              subspecies="",
                              lineage="",
                              sublineage="",
                              other_designation="",
                              notes=notes,
                              stringsAsFactors=F)
  return(list(seqs=mlst_seqs, profiles=mlst_profiles))
}

#' Title
#'
#' @param speciesNames
#' @param table.loci
#' @param threads
#'
#' @return
#' @export
#'
#' @examples
#'
mlstDB.extract <- function(assembly_summary, speciesNames, fasta_dir="~/mlstverse/data/refseq", threads=16) {
  continueSF <- TRUE
  if (!sfIsRunning()) {
    snowfall::sfInit(parallel=T, cpus=threads)
    continueSF <- FALSE
  }
  list.mlst <- list()
  sfExport("getRefSeqCDS", "getRefSeqRNA", "extractProduct", "table.loci", "extractOrthologs")
  for (speciesName in speciesNames) {
    asm <- subset(assembly_summary, species_name==speciesName)
    if (nrow(asm) == 0) {
      cat(paste("Warning: ", speciesName, " not found.\n", sep=""))
      next
    }
    cat(paste("Proccessing ", speciesName, "\n"))
    list.mlst[[speciesName]] <- snowfall::sfApply(asm, 1, createProfile, fasta_dir)
    names(list.mlst[[speciesName]]) <- asm$X..assembly_accession
  }
  if (!continueSF) {
    snowfall::sfStop()
  }
  return(list.mlst)
}

#' Title
#'
#' @param list.mlst
#' @param initial.rST
#' @param table.loci
#' @param initial.locusCounts
#'
#' @return
#' @export
#'
#' @examples
mlstDB.resetIDs <- function(list.mlst, table.loci, initial.locusCounts, initial.rST) {
  rST <- 0
  locusCounts <- initial.locusCounts
  for (species in names(list.mlst)) {
    cat(paste(species, "\n"))
    for (strain in names(list.mlst[[species]])) {
      for (locus in table.loci$mlst_tag) {
        l <- length(list.mlst[[species]][[strain]]$seqs[[locus]])
        if (l > 0) {
          names(list.mlst[[species]][[strain]]$seqs[[locus]]) <-
            locusCounts[table.loci$mlst_tag == locus] + 1:l
          list.mlst[[species]][[strain]]$profiles[[locus]][[1]] <-
            list.mlst[[species]][[strain]]$profiles[[locus]][[1]] +
            locusCounts[table.loci$mlst_tag == locus]
          locusCounts[table.loci$mlst_tag == locus] <-
            locusCounts[table.loci$mlst_tag == locus] + l
        }
      }
      list.mlst[[species]][[strain]]$profiles$rST <-
        initial.rST + rST + list.mlst[[species]][[strain]]$profiles$rST
      rST <- rST + length(list.mlst[[species]][[strain]]$profiles$rST)
    }
  }
  # merge locus sequences
  mlst_seqs <- list()
  tmp <- unlist(lapply(list.mlst, function(x) {lapply(x, "[[", 1)}), recursive=F)
  for (locus in table.loci$mlst_tag) {
    mlst_seqs[[locus]] <- unlist(lapply(tmp, function(x) {x[[locus]]}), recursive=F)
    if (length(mlst_seqs[[locus]]) > 0) {
      names(mlst_seqs[[locus]]) <- gsub(".*\\.", "", names(mlst_seqs[[locus]]))
    }
  }
  # merge profiles
  tmp <- lapply(list.mlst, lapply, "[[", "profiles")
  mlst_profiles <- dplyr::bind_rows(lapply(tmp, bind_rows))
  return(list(seqs=mlst_seqs, profiles=mlst_profiles))
}

#' Title
#'
#' @param speciesNames
#' @param table.loci
#' @param initial.locusCounts
#' @param initial.rST
#' @param outdir
#'
#' @return
#' @export
#'
#' @examples
mlstDB.build <- function(speciesNames, table.loci, initial.locusCounts, initial.rST, outdir=".") {
  assembly_summary <- getAssemblySummary()
  asm <- mlstDB.subsetASM(speciesNames, assembly_summary)
  list.mlst <- mlstDB.extract(asm, speciesNames)
  mlst <- mlstDB.resetIDs(list.mlst, table.loci, initial.locusCounts, initial.rST)
  for (locus in table.loci$mlst_tag) {
    write.fasta(sequences=mlst$seqs[[locus]],
                names=paste(locus, names(mlst$seqs[[locus]]), sep="_"),
                file.out=paste("~/mlstverse/data/locus_seqs/", locus, ".fasta", sep=""))
  }
  mlstverse.db <- mlst$profiles
  save(mlstverse.db, file=paste(outdir, "/mlstverse.db.RData", sep=""))

  return(mlst)
}

#'
#' @param gff_path
#' @param fasta_path
#' @param speciesName
#' @param strain
#' @param start
#'
#' @importFrom tidyr %>%

db.buildCustomDB <- function(gff_path, fasta_path, speciesName, strain, table.loci) {
  mlst_seqs <- array(list(list()), nrow(table.loci))
  names(mlst_seqs) <- table.loci$mlst
  gff <- mlstDB.importGFF(gff_path)
  fasta <- seqinr::read.fasta(fasta_path)
  rp <- extractRPs.RAST(fasta, gff)

  # Add MLST fasta
  if(length(rp)==0) {
    return()
  }
  for (j in seq(table.loci$mlst)) {
    locus <- table.loci$mlst[j]
    product <- table.loci$product[j]
    k <- which(table.loci$product[j]==names(rp))
    if (length(k) == 0) {
      next
    }
    mlst_seqs[[locus]] <- lapply(rp[k], as.character)
    names(mlst_seqs[[locus]]) <- seq(mlst_seqs[[locus]])
  }

  # Add MLST profiles
  m <- max(sapply(mlst_seqs, length))
  profile <- data.frame(lapply(lapply(mlst_seqs, names), function(x) {
    if (is.null(x)) {
      return(as.numeric(rep_len(-1, m)))
    } else {
      return(as.numeric(rep_len(x, m)))
    }
  }))
  profile <- as.data.frame(eval(parse(text=paste("tidyr::expand(profile,", paste(table.loci$mlst, collapse=", "), ")"))))
  notes <- paste("",
                 speciesName,
                 paste("strain=", strain, sep=""),
                 "",
                 "", sep=",")
  mlst_profiles <- data.frame(rST=1:nrow(profile),
                              profile,
                              genus="",
                              species=species,
                              subspecies="",
                              lineage="",
                              sublineage="",
                              other_designation="",
                              notes=notes,
                              stringsAsFactors=F)

  return(list(fasta=mlst_seqs, profiles=mlst_profiles))
}

### removeDBDuplicates
filter.j <- function(x) {
  g.dup <- x[!is.na(x)][-(1:2)]
  which(gids %in% g.dup)
}
updateGIDs <- function(g, gids) {
  tmp <- read.table(paste("db/", g, ".dup", sep=""), header=F, fill=T, stringsAsFactors=F)
  tmp <- apply(tmp, 2, function(x) {gsub("BACT[0-9]+_([0-9]+),*", "\\1", x)})
  tmp <- apply(tmp, c(1,2), as.numeric)

  j <- apply(tmp, 1, filter.j)
  for (i in 1:nrow(tmp)) {
    gids[j[[i]]] <- tmp[i,2]
  }
  gids
}
tmp.func2 <- function() {
  sfInit(parallel=TRUE, cpus=16, type="SOCK")
  for (g in genes) {
    gids <- bigsdb[, g]
    tmp <- read.table(paste("db/", g, ".dup", sep=""), header=F, fill=T, stringsAsFactors=F)
    tmp <- apply(tmp, 2, function(x) {gsub("BACT[0-9]+_([0-9]+),*", "\\1", x)})
    tmp <- apply(tmp, c(1,2), as.numeric)

    cat(paste("Preparing for ", g, ".\n", sep=""))
    sfExport(list=list("tmp", "g", "gids"))

    cat(paste("Filtering duplicated genes...\n"))
    j <- sfApply(tmp, 1, filter.j)

    cat(paste("Overwriting db of ", g, "\n", sep=""))
    for (i in 1:nrow(tmp)) {
      gids[j[[i]]] <- tmp[i,2]
    }
    bigsdb[,g] <- gids
  }


  sfExport(list=list("tmp", "g", "gids"))
  gids <- bigsdb[, g]
  updateGIDs(g, gids)

  j <- apply(tmp, 1, filter.j)
  for (i in 1:nrow(tmp)) {
    gids[j[[i]]] <- tmp[i,2]
  }

  for (g in genes) {
    gids[j[[i]]] <- tmp[i,2]
  }
  bigsdb[,g] <- gids


  sfStop()
}
tmp.func3 <- function() {
  bigsdb.sub <- distinct(bigsdb, .keep_all = T,
                         BACT000001,BACT000002,BACT000003,BACT000004,BACT000005,BACT000006,
                         BACT000007,BACT000008,BACT000009,BACT000010,BACT000011,BACT000012,
                         BACT000013,BACT000014,BACT000015,BACT000016,BACT000017,BACT000018,
                         BACT000019,BACT000020,BACT000021,BACT000030,BACT000031,BACT000032,
                         BACT000033,BACT000034,BACT000035,BACT000036,BACT000038,BACT000039,
                         BACT000040,BACT000042,BACT000043,BACT000044,BACT000045,BACT000046,
                         BACT000047,BACT000048,BACT000049,BACT000050,BACT000051,BACT000052,
                         BACT000053,BACT000056,BACT000057,BACT000058,BACT000059,BACT000060,
                         BACT000061,BACT000062,BACT000063,BACT000064,BACT000065,species)

  i <- which(bigsdb.sub$species%in%c("Mycobacterium indicus", "Mycobacterium indicus pranii "))
  bigsdb.sub$species[i] <- "Mycobacterium indicus pranii"

  #length(unique(subset(bigsdb.sub, grepl("Mycobacterium", species))$species))
  species <- sort(unique(subset(bigsdb.sub, grepl("Mycobacterium", species))$species))
  species.MTBC <- c("Mycobacterium tuberculosis",
                    "Mycobacterium africanum",
                    "Mycobacterium orygis",
                    "Mycobacterium bovis",
                    "Mycobacterium microti",
                    "Mycobacterium canettii",
                    "Mycobacterium caprae",
                    "Mycobacterium pinnipedii",
                    "Mycobacterium suricattae",
                    "Mycobacterium mungi")
  #                  "Mycobacterium leprae",
  #                  "Mycobacterium lepromatosis")
  species.NTM <- species[!species %in% species.MTBC]
  length(species.MTBC)
  length(species.NTM)
  species.NTM
}
