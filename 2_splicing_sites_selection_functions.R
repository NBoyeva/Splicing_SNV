get_sum_vector <- function(df) {
  
  sum_vector <- rowSums(df[, !(colnames(df) %in% 
                                 c("eej_id", "gene_id", "seqnames", 
                                   "start", "end", "width", "strand", 
                                   "intron_length", "sum_samples_cpm", 
                                   "sum_samples_raw", "sum_cpm", 
                                   "sum_raw"))], 
                        na.rm = TRUE)
  return(sum_vector)
}

filter_eej_cols <- function(df,
                            use_raw=FALSE, 
                            raw_threshold=10, 
                            cpm_threshold=-1) {
  
  # This function applies filter_eej_v() to all columns in the dataframe given
  
  if (use_raw) {
    sum_vector <- get_sum_vector(df)
    df_filtered <- as.data.frame(lapply(df, 
                                        function(x) { return(as.integer(x >= raw_threshold)) }))
  } else {
    df_filtered <- as.data.frame(lapply(df, 
                                        cpm, log=TRUE, prior.count=0.5))
    sum_vector <- get_sum_vector(df_filtered)
    df_filtered <- as.data.frame(lapply(df, 
                                        function(x) { return(as.integer(x >= raw_threshold)) }))
  }
  
  return(list(df_filtered, sum_vector))
}



filter_eej_df <- function(full_df,
                          use_raw=FALSE,
                          use_cpm=FALSE,
                          raw_threshold=10, 
                          cpm_threshold=-1, 
                          calculate_cpm_threshold=FALSE,
                          n_samples_threshold=1) {
  
  # This function applies filter_eej_cols() to initial file
  
  meta_df <- full_df[, colnames(full_df) %in% 
                       c("eej_id", "gene_id", "seqnames", 
                       "start", "end", "width", "strand",
                       "intron_length", "sum_samples_cpm", 
                       "sum_samples_raw", "sum_cpm", 
                       "sum_raw")]
  matrix_df <- full_df[, !(colnames(full_df) %in% 
                             c("eej_id", "gene_id", "seqnames", 
                               "start", "end", "width", "strand", 
                               "intron_length", "sum_samples_cpm", 
                               "sum_samples_raw", "sum_cpm", 
                               "sum_raw"))]
  
  if (use_raw) {

    list_filter <- filter_eej_cols(matrix_df,
                                 use_raw=TRUE,
                                 raw_threshold)
    meta_df$sum_raw <- list_filter[[2]]
    meta_df$sum_samples_raw <- get_sum_vector(list_filter[[1]])
    filter_vector <- meta_df$sum_samples_raw >= n_samples_threshold
    meta_df <- meta_df[filter_vector,]
    matrix_df <- matrix_df[filter_vector,]
  }
  
  if (use_cpm) {
    list_filter <- filter_eej_cols(matrix_df,
                                 use_raw=FALSE,
                                 cpm_threshold)
    meta_df$sum_cpm <- list_filter[[2]]
    meta_df$sum_samples_cpm <- get_sum_vector(list_filter[[1]])
    filter_vector <- meta_df$sum_samples_cpm >= n_samples_threshold
    meta_df <- meta_df[filter_vector,]
  }
  
  return(meta_df)
}