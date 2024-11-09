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
  
  ss_overlapped$snp_id <- c(ss5_snp@ranges@NAMES, 
                            ss3_snp@ranges@NAMES)
  ss_overlapped$snp_pos <- c(start(ss5_snp@ranges),
                             start(ss3_snp@ranges))
  ss_overlapped$snp_ref <- c(as.character(ss5_snp$REF),
                             as.character(ss3_snp$REF))
  ss_overlapped$snp_alt <- c(unlist(lapply(ss5_snp$ALT, as.character)),
                             unlist(lapply(ss3_snp$ALT, as.character)))
  ss_overlapped$snp_ss <- c(rep("5", length(ss5_snp)), 
                            rep("3", length(ss3_snp)))
  
  ss_overlapped <- ss_overlapped[!duplicated(ss_overlapped),]
  
  return(ss_overlapped)
}
