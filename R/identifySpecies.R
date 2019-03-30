#' read SAM format file
#'
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
readSAM <- function(filename) {
  sam <- suppressMessages(suppressWarnings(readr::read_delim(filename, delim="\t", col_names=F, quote='', progress=F, comment="@")))
  colnames(sam) <- c("qname", # Query template NAME
                     "flag",  # bitwise FLAG
                     "rname", # References sequence NAME
                     "pos",   # 1 - based leftmost mapping POSition
                     "mapq",  # MAPping Quality
                     "cigar", # CIGAR String
                     "rnext", # Ref. name of the mate/next read
                     "pnext", # Position of the mate/next read
                     "tlen",  # observed Template LENgth
                     rep("", ncol(sam) - 9))
  return(sam[,1:11])
}

#' create an input query for MLST score calculation
#'
#' @param sam
#'
#' @return
#' @export
#'
#' @examples
createQuery <- function(loci, depth, len.loci, threads=1) {
  continueSF <- TRUE
  if (!snowfall::sfIsRunning()) {
    suppressMessages(snowfall::sfInit(parallel=T, cpus=threads, useRscript=T))
    continueSF <- FALSE
  }
  query <- snowfall::sfLapply(loci, function(l, depth, len.loci) {
    q <- data.frame()
    depth.sub <- subset(depth, grepl(l, seqnames))
    if (nrow(depth.sub) > 0) {
      gids <- as.character(unique(depth.sub$seqnames))
      for (g in gids) {
        tmp <- subset(depth, g==seqnames)
        q <- rbind(q,
                   data.frame(locus_tag=gsub(paste(l,"_",sep=""), "", g),
                              depth=sum(tmp$count),
                              coverage=sum(tmp$count)/len.loci[g],
                              coverRatio=as.vector(nrow(tmp)/len.loci[g]),
                              len=len.loci[g],
                              mapped=nrow(tmp)))
      }
    }
    return(q)
  }, depth, len.loci)
  if (!continueSF) {
    snowfall::sfStop()
  }
  names(query) <- loci
  return(query)
}

#' calculate MLST score
#'
#' @param threads
#' @param query
#' @param mlstdb
#'
#' @return
#' @export
#'
#' @examples
calcMLSTScore <- function(query,
                          loci,
                          mlstdb=mlstverse.NTM.db,
                          method="default",
                          normalize=TRUE,
                          threads=1) {
  continueSF <- TRUE
  if (!snowfall::sfIsRunning()) {
    suppressMessages(snowfall::sfInit(parallel=T, cpus=threads, useRscript=T))
    continueSF <- FALSE
  }
  g <- function(db_entry, query) {
    db_entry <- lapply(db_entry, "[[", 1)
    f <- function(l, db_entry, query) {
      if (-1 %in% db_entry[[l]] | nrow(query[[l]]) == 0) {
        return(0)
      }
      #found <- query[[l]]$locus_tag %in% db_entry[[l]]
      found <- db_entry[[l]] %in% query[[l]]$locus_tag
      if (any(found)) {
        if (method=="default") {
          return(mean(subset(query[[l]], locus_tag %in% db_entry[[l]])$coverRatio) / length(db_entry[[l]]))
          #return(max(subset(query[[l]], locus_tag %in% db_entry[[l]])$coverRatio))
          #return(sum(subset(query[[l]], locus_tag %in% db_entry[[l]])$coverRatio) / length(db_entry[[l]]))
        } else if (method=="sensitive") {
          return(nrow(subset(query[[l]], locus_tag %in% db_entry[[l]])) / length(db_entry[[l]]))
        }
      } else {
        return(0)
      }
    }
    tmp <- sapply(names(db_entry), f, db_entry, query)
    scoreLimit <- sum(!sapply(db_entry, function(x) {-1%in%x}))
    if (scoreLimit > 0) {
      if (normalize) {
        return(sum(tmp, na.rm=T) / scoreLimit)
      } else {
        return(sum(tmp, na.rm=T))
      }
    } else {
      return(0)
    }
  }
  scores <- snowfall::sfApply(mlstdb[,loci], 1, g, query)
  if (!continueSF) {
    snowfall::sfStop()
  }
  return(scores)
}

#' get read count
#'
#' @param entry
#' @param query
#' @param loci
#' @param method
#' @param fill
#'
getCounts <- function(entry, query, loci, method="coverage", fill=TRUE) {
  x <- c()
  for (l in loci) {
    isFound <- FALSE
    if (nrow(query[[l]]) > 0) {
      tmp <- subset(query[[l]], locus_tag %in% as.character(entry[[l]]))
      if (nrow(tmp) > 0) {
        isFound <- TRUE
        if (method=="coverage") {
          x <- c(x, mean(tmp$coverage)/nrow(tmp))
        } else if (method=="mapped") {
          x <- c(x, mean(tmp$mapped)/nrow(tmp))
        }
      }
    }
    if (!isFound) {
      x <- c(x, NA)
    }
  }
  if (fill) {
    i <- is.na(x) & !sapply(entry[loci], function(x) {-1%in%x})
    x[i] <- 0
  }
  return(x)
}


obtainConsensus <- function(bam, query, locus) {

}


#' Title
#'
#' @param filenames
#' @param mlstdb database for MLST (default: mlstverse.NTM.db)
#' @param min_depth only use genes larger than the minimum reads depth (default: 0)
#' @param min_ratio only use genes larger than the ratio to maximum reads depth (default: 0.1)
#' @param th.pvalue filter by threshold value in Kolmogorov–Smirnov test (default: 0.05)
#' @param th.score filter by threshold value in MLST score (default: 0.1)
#' @param threads number of threads (default: 1)
#' @param normalize boolean, if TRUE, normalize coverage (default: TRUE)
#' @param samfile return value of readSAM() (default: NULL)
#'
#' @return
#' @export
#'
#' @examples
mlstverse <- function(filenames,
                      mlstdb=mlstverse.Mycobacterium.db::mlstverse.Mycobacterium.db,
                      min_depth=0,
                      min_ratio=0.1,
                      th.pvalue=0.05,
                      th.score=0.1,
                      threads=1,
                      normalize=TRUE,
                      samfile=NULL,
                      method="default") {
  suppressWarnings(suppressMessages(snowfall::sfInit(parallel=T, cpus=threads, useRscript=T)))
  #loci <- colnames(mlstdb)[grep("BACT[0-9]+", colnames(mlstdb))]
  loci <- colnames(mlstdb)[grep("Locus_[0-9]+", colnames(mlstdb))]

  score <- list()
  query <- list()
  for (filename in filenames) {
    cat(paste("Processing", filename, "\n"))
    cat(paste("  Loading bam file...\n"))
    if (is.null(samfile)) {
      pileup.param <- Rsamtools::PileupParam(max_depth=10000,
                              min_base_quality=0,
                              min_mapq=0,
                              min_nucleotide_depth=0,
                              min_minor_allele_depth=0,
                              distinguish_strands=FALSE,
                              distinguish_nucleotides=FALSE,
                              ignore_query_Ns=TRUE,
                              include_deletions=TRUE,
                              include_insertions=TRUE,
                              left_bins=NULL,
                              query_bins=NULL,
                              cycle_bins=NULL)
      depth <- Rsamtools::pileup(filename, pileupParam=pileup.param)
      len.loci <- Rsamtools::scanBamHeader(filename)[[1]]$targets
    } else {
      sam <- samfile
    }
    #sfExport("depth", "len.loci", "loci")
    cat(paste("  Generating query...\n"))
    query[[filename]] <-
      lapply(createQuery(loci, depth, len.loci),
             function(q) {
               if (nrow(q) > 0) {
                 min_d <- max(max(q$depth)*min_ratio, min_depth)
                 subset(q, min_d < depth)
                 } else {
                   return(q)
                   }
               })
    cat(paste("  Calculating MLST score...\n"))
    snowfall::sfExport("normalize", "method")
    results <- calcMLSTScore(query[[filename]], loci, mlstdb=mlstdb, threads=threads, method=method, normalize=normalize)
    #results <- calcMLSTScore(query[[filename]], loci, mlstdb=mlstdb)
    #results <- calcMLSTScore(query[[filename]], loci, mlstdb=mlstdb, method="sensitive", threads=threads)
    #results <- calcMLSTScore(query[[filename]][loci[1:53]], loci[1:53], mlstdb=mlstdb, threads=threads)

    if (normalize) {
      i <- results > th.score
    } else {
      i <- results > max(results)*th.score
    }
    speciesNames <- unique(apply(
      mlstdb[i, c("genus", "species")], 1, paste, collapse="|"))
    if (length(speciesNames) == 0) {
      score[[filename]] <- NULL
      next
    }

    tmp <- list()

    tmp$mean <- c()
    tmp$var <- c()
    tmp$dist <- c()
    tmp$score <- c()
    tmp$pvalue <- c()
    tmp$strains <- c()
    cat("  Starting Kolmogorov–Smirnov test...\n")
    for (speciesName in speciesNames) {
      cat(paste("    testing for", speciesName, "\n"))
      x <- strsplit(speciesName, "\\|")[[1]]
      if (is.na(x[2])) {
        j <- which(i & mlstdb$genus==x[1] & mlstdb$species=="")
      } else {
        j <- which(i & mlstdb$genus==x[1] & mlstdb$species==x[2])
      }
      #q.dist <- apply(mlstdb[j,loci], 1, getCounts, q)
      #sfExport("loci", "len.loci")
      q.dist <- snowfall::sfApply(mlstdb[loci][j,], 1, getCounts, query[[filename]], loci, fill=TRUE)

      q.mean <- apply(q.dist, 2, mean, na.rm=T)
      q.var <- apply(q.dist, 2, var, na.rm=T)
      snowfall::sfExport("q.dist", "q.mean", "q.var")

      tmp$ks.test <- snowfall::sfLapply(
        1:ncol(q.dist),
        function(i) {
          m <- q.mean[i]
          v <- q.var[i]
          fit <- MASS::fitdistr(round(na.omit(q.dist[,i])),
                          densfun="negative binomial",
                          start=list(mu=m, size=m^2/(m+v)))
          ks.test(x=q.dist[,i],
                  y="pnbinom",
                  size=fit$estimate[2],
                  prob=fit$estimate[2]/fit$estimate[1])
          #ks.test(x=q.dist[,i], y="pgamma", shape=m^2/v, scale=v/m)
          ##          ks.test(x=q.dist[,i],
          ##                  y="pnbinom",
          ##                  size=q.mean[i]^2/(q.mean[i]+q.var[i]),
          ##                  prob=q.mean[i]/(q.mean[i]+q.var[i]))
        }
      )
      tmp$ks.pvalue <- sapply(tmp$ks.test, "[[", "p.value")

      table.score <- data.frame(index=j, score=results[j], pvalue=tmp$ks.pvalue)
      table.score.p <- subset(table.score, pvalue > th.pvalue | pvalue == max(table.score$pvalue))
      table.score.p.s <- subset(table.score.p, score==max(score))
      l <- table.score.p.s$index[order(table.score.p.s$pvalue, decreasing=T)[1]]

      tmp$mean <- c(tmp$mean, q.mean[j==l])
      tmp$var <- c(tmp$var, q.var[j==l])
      tmp$dist <- c(tmp$dist, list(q.dist[,j==l]))
      tmp$pvalue <- c(tmp$pvalue, tmp$ks.pvalue[j==l])
      tmp$score <- c(tmp$score, results[l])
      tmp$strains <- c(tmp$strains,
                       gsub("strain=",
                            "",
                            lapply(strsplit(mlstdb$notes[l], ","), "[", 3)))
    }
    if (any(is.na(tmp$strains))) {
      tmp$strains[is.na(tmp$strains)] <- ""
    }
    score[[filename]] <-
      dplyr::data_frame(genus=gsub("\\|.*", "", speciesNames),
                 species=gsub(".*\\|", "", speciesNames),
                 strain=tmp$strains,
                 score=tmp$score,
                 p.value=tmp$pvalue,
                 mean=tmp$mean,
                 var=tmp$var)
  #               dist=tmp$dist)
    if (nrow(score[[filename]]) > 0) {
      score[[filename]] <- score[[filename]][order(score[[filename]]$p.value, decreasing=T),]
      score[[filename]] <- score[[filename]][order(score[[filename]]$score, decreasing=T),]
    }
  }
  snowfall::sfStop()
  return(list(query=query, score=score))
}



