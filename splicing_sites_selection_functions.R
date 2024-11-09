filter_eej_v <- function(v, 
                         use_raw=FALSE, 
                         raw_threshold=10, 
                         cpm_threshold=-1, 
                         calculate_cpm_threshold=FALSE) {
  
  # This function takes vector, calculates cpm and threshold for filtration
  
  if (use_raw) {
    return(as.integer(v >= raw_threshold))
  } else {
    dens <- density(cpm(v, log=TRUE, prior.count=0.5))
    der <- c(diff(dens$y), tail(diff(dens$y), 1)) > 0
    if (calculate_cpm_threshold) {
      cpm_threshold <- tail(dens$x[c(FALSE, diff(der) != 0)], 2)[1]
    }
    return(as.integer(cpm(v, log=TRUE, prior.count=0.5) >= cpm_threshold))
  }
}



filter_eej_cols <- function(df,
                            use_raw=FALSE, 
                            raw_threshold=10, 
                            cpm_threshold=-1, 
                            calculate_cpm_threshold=FALSE) {
  
  # This function applies filter_eej_v() to all columns in the dataframe given
  
  cols <- list()
  for (col in colnames(df)) {
    col_filter <- filter_eej_v(df[col],
                               use_raw, 
                               raw_threshold, 
                               cpm_threshold, 
                               calculate_cpm_threshold)
    cols[[col]] <- col_filter
  }
  df <- data.frame(cols)
  return(df)
}



filter_eej_df <- function(full_df,
                          use_raw=FALSE, 
                          raw_threshold=10, 
                          cpm_threshold=-1, 
                          calculate_cpm_threshold=FALSE) {
  
  # This function applies filter_eej_cols() to initial file
  
  meta_df <- full_df[, colnames(combined_df) %in% 
                       c("eej_id", "gene_id", "seqnames", 
                       "start", "end", "width", "strand", "intron_length")]
  matrix_df <- full_df[, !(colnames(combined_df) %in% 
                         c("eej_id", "gene_id", "seqnames", 
                        "start", "end", "width", "strand", "intron_length"))]
  df_filter <- filter_eej_cols(matrix_df,
                               use_raw,
                               raw_threshold,
                               cpm_threshold,
                               calculate_cpm_threshold)
  df <- cbind(meta_df, df_filter)
  return(df)
}



add_sum_column <- function(df) {
  
  # This function adds to the dataframe a column with number of samples where EEJ passed the threshold
  
  df$sum_samples <- rowSums(df[, !(colnames(combined_df) %in% 
                                    c("eej_id", "gene_id", "seqnames", 
                                    "start", "end", "width", "strand", 
                                    "intron_length"))], 
                            na.rm = TRUE)
  return(df)
}



get_eej <- function(df,
                    use_raw=FALSE, 
                    raw_threshold=10, 
                    cpm_threshold=-1, 
                    calculate_cpm_threshold=FALSE) {
  
  # This function returns subset of the dataframe with EEJ which passed the threshold at least in one sample
  
  df_filter <- filter_eej_df(df,
                             use_raw, 
                             raw_threshold, 
                             cpm_threshold, 
                             calculate_cpm_threshold)
  df_filter <- add_sum_column(df_filter)
  df_filtered <- df_filter[df_filter$sum_samples > 0,]
  
  print(paste0('All EEJs: ', nrow(df_filter)))
  print(paste0('Filtered EEJs: ', nrow(df_filtered)))
  print(paste0('Filtered out: ', nrow(df_filter) - nrow(df_filtered)))
  
  return(df_filtered)
}