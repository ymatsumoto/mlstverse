getGenomicFasta <- function(ftp_path, fasta_dir="temporary") {
  tmp <- strsplit(ftp_path, "/")[[1]]
  refseqID <- tmp[length(tmp)]
  filename <- paste(refseqID, "_genomic.fna.gz", sep="")
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
}

getGFF <- function(ftp_path, gff_dir="temporary") {
  tmp <- strsplit(ftp_path, "/")[[1]]
  refseqID <- tmp[length(tmp)]
  filename <- paste(refseqID, "_genomic.gff.gz", sep="")
  if (fasta_dir == "temporary") {
    gff_path <- tempfile()
    url <- paste(ftp_path, "/", filename, sep="")
    download.fasta(url, gff_path)
  } else {
    gff_path <- paste(gff_dir, "/", filename, sep="")
    if (!file.exists(gff_path)|file.size(gff_path)==0) {
      url <- paste(ftp_path, "/", filename, sep="")
      download.fasta(url, gff_path)
    }
  }
}


getCDSFasta <- function(ftp_path, fasta_dir="temporary") {
  tmp <- strsplit(ftp_path, "/")[[1]]
  refseqID <- tmp[length(tmp)]
  for (suffix in c("_cds_from_genomic.fna.gz",
                   "_protein.faa.gz",
                   "_rna_from_genomic.fna.gz")) {
    filename <- paste(refseqID, suffix, sep="")
    fdir <- paste(fasta_dir,
                  gsub("_.*|.fa.*", "", gsub("_(.*)", "\\1", suffix)),
                  sep="/")
    if (ftp_path == "na") {
      next
    }
    if (fasta_dir == "temporary") {
      fasta_path <- tempfile()
      url <- paste(ftp_path, "/", filename, sep="")
      download.fasta(url, fasta_path)
    } else {
      fasta_path <- paste(fdir, "/", filename, sep="")
      if (!file.exists(fasta_path)|file.size(fasta_path)==0) {
        url <- paste(ftp_path, "/", filename, sep="")
        download.fasta(url, fasta_path)
      }
    }
  }
}

#asm <- read.table("assembly_summary_refseq_20171201.txt", sep="\t", header=T, stringsAsFactors=F, skip=1, fill=T, comment.char="", quote="")
#asm.sub <- subset(asm, grepl("^Mycobacterium", organism_name))
#asm.sub2 <- subset(asm.sub, !grepl("Mycobacterium phage", organism_name))
#asm.sub3 <- subset(asm.sub2, !grepl("Mycobacterium virus", organism_name))
#cds_dir <- "~/mlstverse/data/CDS"
#for (ftp_path in asm.sub3$ftp_path) {
#  print(which(ftp_path==asm.sub3$ftp_path))
#  getCDSFasta(ftp_path, cds_dir)
#}

#for (ftp_path in asm.sub3$ftp_path) {
#  print(which(ftp_path==asm.sub3$ftp_path))
#  getGenomicFasta(ftp_path, "~/mlstverse/data/genome/")
#}

