####### Read VCF files into a list of dataframes #######
create_full_info_table <- function(filepaths) {
  
  info_tables <- list()
  
  for (i in seq_along(filepaths)) {
    vcf <- readVcf(filepaths[i], "hg38")
    info_table <- info(vcf)
    info_tables[[i]] <- as.data.frame(info_table)
  }
  
  return(info_tables)
}


####### Plot given parameter densities for all VCFs represented as infotables #######
plot_density <- function(info_tables, 
                         param_name, 
                         labels, 
                         descriptions, 
                         x_log10 = FALSE, 
                         y_log10 = FALSE) {
  
  data_list <- list()
  
  for (i in seq_along(info_tables)) {
    
    if (param_name %in% names(info_tables[[i]])) {
      parameter_data <- unlist(info_tables[[i]][[param_name]])  
      
      if (length(parameter_data) > 0) {
        data_list[[i]] <- data.frame(
          Value = parameter_data,
          File = labels[i]
        )
      }
    } else {
      warning(paste("Error: parameter", param_name, "not found in", labels[i]))
    }
  }
  
  if (length(data_list) == 0) {
    stop("Error: no valid data found for the specified parameter")
  }
  
  combined_data <- do.call(rbind, data_list)
  
  x_label <- if (x_log10) paste(param_name, "(log10 scale)") else param_name
  y_label <- if (y_log10) "Density (log10 scale)" else "Density"
  
  p <- ggplot(combined_data, aes(x = Value, color = File)) +
    geom_density() +
    labs(title = descriptions[param_name, ]$Description,
         x = x_label,
         y = y_label) +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  if (x_log10) {
    p <- p + scale_x_log10()
  }
  if (y_log10) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}


####### Apply VCF filtration and return filtration statistics (optionally) #######
# Prior to this, VCF should be processed with FilterVariants from GATK.
hardfilter_vcf <- function(unfiltered_vcf, stats=FALSE) {
  
  if (stats) {
    print(paste0("Filtered by FS:  ", 
                 as.character(length(
                   unfiltered_vcf[unfiltered_vcf$FILTER == "FS"]))))
    print(paste0("Filtered by SOR: ", 
                 as.character(length(
                   unfiltered_vcf[unfiltered_vcf$FILTER == "SOR"]))))
    print(paste0("Filtered by QUAL:  ", 
                 as.character(length(
                   unfiltered_vcf[unfiltered_vcf$FILTER == "QUAL"]))))
    print(paste0("Filtered by MQ:  ", 
                 as.character(length(
                   unfiltered_vcf[unfiltered_vcf$FILTER == "MQ"]))))
    print(paste0("Filtered by ReadPosRankSum:  ", 
                 as.character(length(
                   unfiltered_vcf[unfiltered_vcf$FILTER == "ReadPosRankSum"]))))
    print(paste0("Filtered by MQRankSum: ", 
                 as.character(length(
                   unfiltered_vcf[unfiltered_vcf$FILTER == "MQRankSum"]))))
  }
  return(unfiltered_vcf[unfiltered_vcf$FILTER == "PASS"])
}

load_vcf <- function(path, hardfilter=TRUE, stats=FALSE) {
  vcf_file <- readVcf(path, "hg38")
  gr_vcf <- granges(vcf_file)
  if (hardfilter) {
    gr_vcf <- hardfilter_vcf(gr_vcf, stats=stats)
  }
  return(gr_vcf)
}


####### Extract the INFO tables and SNV positions from VCF files #######
extract_info_and_positions <- function(filepath) {
  
  vcf <- readVcf(filepath, "hg38")
  info_table <- info(vcf)
  snv_positions <- rownames(vcf)
  
  return(list(info_table = info_table, 
              positions = snv_positions))
}


####### Identify unique and common SNVs based on positions #######
identify_snv_sets_positions <- function(snv_positions, 
                                        file_labels) {

  unique_snvs <- 
    lapply(seq_along(snv_positions), 
           function(i) {
             all_other_snvs <- unlist(snv_positions[-i])
             unique_positions <- snv_positions[[i]][!(snv_positions[[i]] 
                                                      %in% all_other_snvs)]
             return(data.frame(Position = unique_positions, 
                               File = paste(file_labels[i], "unique")))
             })
  
  # Identify common SNVs in genomic VCFs 
  gen_files <- c("gen1", "gen2", "gen_merged")
  gen_indices <- which(file_labels %in% gen_files)
  gen_common_positions <- Reduce(intersect, 
                                 lapply(gen_indices, 
                                        function(i) snv_positions[[i]]))
  gen_common_snvs <- if (length(gen_common_positions) > 0) {
    data.frame(Position = gen_common_positions, File = "gen common")
  } else {
    NULL
  }
  
  # Identify common SNVs in nascent VCFs
  nas_files <- c("nas1", "nas2", "nas3", "nas_merged")
  nas_indices <- which(file_labels %in% nas_files)
  nas_common_positions <- Reduce(intersect, 
                                 lapply(nas_indices, 
                                        function(i) snv_positions[[i]]))
  nas_common_snvs <- if (length(nas_common_positions) > 0) {
    data.frame(Position = nas_common_positions, File = "nas common")
  } else {
    NULL
  }
  
  # Identify common SNVs across all files
  all_common_positions <- Reduce(intersect, snv_positions)
  all_common_snvs <- if (length(all_common_positions) > 0) {
    data.frame(Position = all_common_positions, File = "all common")
  } else {
    NULL
  }
  
  # Combine all sets, removing NULL entries
  combined_positions <- bind_rows(
    unique_snvs,
    gen_common_snvs,
    nas_common_snvs,
    all_common_snvs
  ) %>% filter(!is.null(.))
  
  return(combined_positions)
}


####### Helper function to extract common metrics for "gen common", "nas common", and "all common" #######
extract_common_metric <- function(info_tables, param_name, snv_set) {
  positions <- snv_set$Position  # SNV positions from the common SNV set

  common_values <- rep(NA, length(positions))
  
  for (i in seq_along(info_tables)) {
    
    if (param_name %in% names(info_tables[[i]]$info_table)) {
      # Extract all parameter values for this file
      param_values <- unlist(info_tables[[i]]$info_table[[param_name]])
      
      # Get the positions that match the common SNV positions
      matched_indices <- match(positions, info_tables[[i]]$positions)
      
      # Extract the parameter values corresponding to the matched positions
      common_values <- ifelse(!is.na(matched_indices), 
                              param_values[matched_indices], 
                              common_values)
    } else {
      next
    }
  }
  
  if (length(common_values) != length(positions)) {
    stop("Mismatch between common parameter values 
         and positions in the common SNV set.")
  }
  
  return(common_values)
}



####### Extract parameter values for SNV sets using stored infotables #######
extract_metric_from_info_tables <- function(info_tables, 
                                            param_name, 
                                            snv_set) {
  
  param_values_list <- list()
  for (i in seq_along(info_tables)) {
    
    # Get the SNV positions for this file in the SNV set
    positions_in_file <- snv_set$Position[snv_set$File == paste0(labels[i], " unique")]
    
    # Check if the parameter exists in the INFO table
    if (param_name %in% names(info_tables[[i]]$info_table)) {
      
      # Extract all the parameter values from the INFO table
      param_values <- unlist(info_tables[[i]]$info_table[[param_name]])
      
      # Filter the parameter values based on the SNV positions
      matched_indices <- match(positions_in_file, 
                               info_tables[[i]]$positions)
      matched_values <- param_values[matched_indices]
      
      # Ensure that the number of values matches the number of positions in the SNV set
      if (length(matched_values) != length(positions_in_file)) {
        stop("Mismatch between the number of matched 
             parameter values and positions in the SNV set.")
      }
      
      # Add the matched values to the list
      param_values_list[[i]] <- matched_values
      
    } else {
      
      # If the parameter doesn't exist, return NA values for the positions
      param_values_list[[i]] <- rep(NA, length(positions_in_file))
      
    }
  }
  
  unique_values <- unlist(param_values_list)
  
  # Update the SNV set with the extracted metric values for unique SNVs
  snv_set$Value[snv_set$File %in% paste0(labels, " unique")] <- unique_values
  
  # Extract and update metrics for the common SNV sets 
  common_files <- c("gen common", "nas common", "all common")
  for (common_file in common_files) {
    common_indices <- which(snv_set$File == common_file)
    
    if (common_file == "gen common") {
      snv_set$Value[common_indices] <- extract_common_metric(info_tables[1:3], 
                                                             param_name, 
                                                             snv_set[common_indices, ])
    } else if (common_file == "nas common") {
      snv_set$Value[common_indices] <- extract_common_metric(info_tables[4:7], 
                                                             param_name, 
                                                             snv_set[common_indices, ])
    } else if (common_file == "all common") {
      snv_set$Value[common_indices] <- extract_common_metric(info_tables, 
                                                             param_name, 
                                                             snv_set[common_indices, ])
    }
  }
  
  return(snv_set)
}


####### Plot the SNV data with different parameter metrics #######
plot_snv_density <- function(snv_set_with_metric, 
                             param_name, 
                             descriptions, 
                             x_log10 = FALSE, 
                             y_log10 = FALSE) {
  
  x_label <- if (x_log10) paste(param_name, "(log10 scale)") else param_name
  y_label <- if (y_log10) "Density (log10 scale)" else "Density"
  
  snv_set_with_metric$LineThickness <- ifelse(snv_set_with_metric$File %in% 
                                                c("gen common", 
                                                  "nas common", 
                                                  "all common"), 1.5, 0.5)
  
  p <- ggplot(snv_set_with_metric, aes(x = Value, 
                                       color = File, 
                                       size = LineThickness)) +
    geom_density() +
    labs(title = descriptions[[param_name, "Description"]],
         x = x_label,
         y = y_label) +
    theme_minimal() +
    scale_size_identity() +
    theme(legend.title = element_blank())
  
  if (x_log10) {
    p <- p + scale_x_log10()
  }
  if (y_log10) {
    p <- p + scale_y_log10()
  }
  
  return(p)
}
