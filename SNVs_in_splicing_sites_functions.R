get_SS_from_EEJ <- function(read_from_file=TRUE, filepath=NULL, df=NULL){
  
  if (read_from_file) {
    EEJs <- read.table(file=filepath,
                       sep="\t",
                       header=TRUE,
                       quote="\"",
                       as.is=TRUE)
  } else {
    EEJs <- makeGRangesFromDataFrame(df=df,
                                     keep.extra.columns=TRUE)
  }
  
  # Convert genomic coordinates of exon-exon junctions
  #   in genomic coordinates of splice sites.
  fiveSSs <- EEJs
  start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) - 2
  end(x=fiveSSs[strand(x=fiveSSs) == "+", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "+", ]) + 8
  start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
    end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) - 6
  end(x=fiveSSs[strand(x=fiveSSs) == "-", ]) <-
    start(x=fiveSSs[strand(x=fiveSSs) == "-", ]) + 8
  threeSSs <- EEJs
  start(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
    end(x=threeSSs[strand(x=threeSSs) == "+", ]) - 20
  end(x=threeSSs[strand(x=threeSSs) == "+", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "+", ]) + 22
  start(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) - 2
  end(x=threeSSs[strand(x=threeSSs) == "-", ]) <-
    start(x=threeSSs[strand(x=threeSSs) == "-", ]) + 22
  SSs <- list(fiveSSs=fiveSSs, threeSSs=threeSSs)
  
  return(SSs)
}

find_overlaps_GNC <- function(ss_granges, vcf_granges) {
  
  # This function creates Nested Containment List  for
  # GRanges object with SS data and uses it to find overlaps
  # with GRanges-formatted data from VCF.
  #
  # Args:
  # ss_granges : GRanges object
  #   splicing sites data
  # vcf_granges : GRanges object
  #   SNP data
  
  ss_git <- GNCList(ss_granges)
  ss_vcf_overlap <- findOverlaps(vcf_granges, ss_git)
  return(ss_vcf_overlap)
}

find_overlaps_jointSS <- function(ss_coords, vcf_granges, ss_df) {
  
  ss5_overlap <- find_overlaps_GNC(ss_coords$fiveSSs, vcf_granges)
  ss3_overlap <- find_overlaps_GNC(ss_coords$threeSSs, vcf_granges)
  
  ss5_snp <- vcf_granges[ss5_overlap@from,]
  ss3_snp <- vcf_granges[ss3_overlap@from,]
  
  ss_overlapped_idx <- c(ss5_overlap@to, ss3_overlap@to)
  ss_overlapped <- ss_df[ss_overlapped_idx,]
  
  ss_overlapped$ss <- c(rep("5", length(ss5_snp)), 
                        rep("3", length(ss3_snp)))
  
  #ss_overlapped <- ss_overlapped[!duplicated(ss_overlapped),]
  
  ss5 <- ss_coords$fiveSSs[match(ss_df[ss5_overlap@to,]$eej_id, ss_coords$fiveSSs$eej_id)]
  ss3 <- ss_coords$threeSSs[match(ss_df[ss3_overlap@to,]$eej_id, ss_coords$threeSSs$eej_id)]
  
  ss_overlapped$ss_start <- c(ss5@ranges@start, 
                              ss3@ranges@start)
  ss_overlapped$ss_end <- c(ss5@ranges@start +
                              ss5@ranges@width -
                              1,
                            ss3@ranges@start +
                              ss3@ranges@width -
                              1)
  
  ss_overlapped$snp_id <- c(ss5_snp@ranges@NAMES, 
                            ss3_snp@ranges@NAMES)
  ss_overlapped$snp_pos <- c(start(ss5_snp@ranges),
                             start(ss3_snp@ranges))
  ss_overlapped$snp_ref <- c(as.character(ss5_snp$REF),
                             as.character(ss3_snp$REF))
  
  ss5_snp_alt <- sapply(ss5_snp$ALT, 
                        function(x) paste(as.character(x), collapse = ","))
  ss3_snp_alt <- sapply(ss3_snp$ALT, 
                        function(x) paste(as.character(x), collapse = ","))
  
  ss_overlapped$snp_alt <- c(ss5_snp_alt,
                             ss3_snp_alt)
  
  ss_overlapped$snp_pos_in_ss <- ss_overlapped$snp_pos - ss_overlapped$ss_start

  ss_overlapped <- ss_overlapped[!duplicated(ss_overlapped),]
  
  return(ss_overlapped)
}


add_refseqs <- function(ref, df) {
  ss_gr <- GRanges(IRanges(start = df$ss_start, end = df$ss_end), 
                   seqnames = df$seqnames, 
                   strand = rep("+", length(df$strand)))
  df$refseq <- unname(as.character(getSeq(ref, ss_gr)))
  
  return(df)
}


add_altseqs <- function(df) {
  
  # Add the altseq column
  df$altseq <- apply(df, 1, function(row) {
    # Extract row values
    snp_pos <- as.numeric(row["snp_pos_in_ss"])
    snp_ref <- row["snp_ref"]
    snp_alt <- row["snp_alt"]
    refseq <- row["refseq"]
    
    # Check if snp_pos is within bounds
    if (snp_pos < 0 || snp_pos >= nchar(refseq)) {
      warning(sprintf("snp_pos_in_ss out of bounds for eej_id %s. Setting altseq to NA.", row["eej_id"]))
      return(NA)
    }
    
    # Check if refseq nucleotide matches snp_ref
    if (substr(refseq, snp_pos + 1, snp_pos + 1) != snp_ref) {
      warning(sprintf(
        "Mismatch in refseq and snp_ref at eej_id %s: refseq[%d] = %s, expected %s. Setting altseq to NA.",
        row["eej_id"], snp_pos + 1, substr(refseq, snp_pos + 1, snp_pos + 1), snp_ref
      ))
      return(NA)
    }
    
    # Replace the nucleotide and return the new sequence
    altseq <- refseq
    substr(altseq, snp_pos + 1, snp_pos + 1) <- snp_alt
    return(altseq)
  })
  
  return(df)
}