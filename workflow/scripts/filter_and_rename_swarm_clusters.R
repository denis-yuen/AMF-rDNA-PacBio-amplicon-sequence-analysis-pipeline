#!/usr/bin/env Rscript

# This script takes a fasta file as input and applies
# a filter based on a set of specific conditions related
# to the abundance of sequences. It returns a filtered
# fasta file.

# Filtering conditions:
# if unique abundance values < 7 and all abundances = 1, all sequences are removed
# if unique abundance values < 7, but not all abundance = 1, only keep sequence with highest abundance
# if unique abundance values > 7, run 7-cluster K-means, keep top-three groups if max abundance > 150
# if unique abundance values > 7, run 7-cluster K-means, keep top-two groups if max abundance < 150
 
# set args
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# get args (input file name and output dir)
input_file <- args[1]
output_dir <- args[2]

# extract sample name and set output suffixes and ouput file names
input_file_suffix <- ".tax_filtered.swarm.seeds.fasta"
sample <- gsub(input_file_suffix, "", basename(input_file))
output_file_suffix <- ".abd_filtered.swarm.seeds.fasta"
ouput_report_suffix <- ".swarm.kmeans.tsv"
output_filename <- paste0(output_dir,"/", sample, output_file_suffix)
output_reportname <- paste0(output_dir,"/", sample, ouput_report_suffix)

# read input fasta
# if file is not empty, continue
if (file.size(input_file) != 0L) { 
  # read input_file
  input <- read.csv(input_file, header=F)
  # merge header and sequence on the same row
  merged_input <- sapply(seq(1, nrow(input), by = 2), function(i) {
    if (i + 1 <= nrow(input)) {
      paste(input$V1[i], input$V1[i + 1], sep = "&")
    } else {
      input$Value[i]  # If there's an odd row, keep it as is
    }
  })
  # from merged_input, create a vector of all sizes
  sizes <- as.numeric(gsub(".*size=(\\d+);.*", "\\1", merged_input))
  # Handle cases where there are fewer than 7 unique abundance values
  if (length(unique(sizes)) < 7) {
    if (all(sizes == 1)) {
      # Skip entirely if all sequence abundances are 1
      cat("Skipping", basename(input_file), "- All sequences have abundance = 1.\n")
      # Write empty files 
      system(paste0("touch ", output_filename))
      system(paste0("touch ", output_reportname))
    } else {
      # Keep only the highest abundance sequence
      significant_numbers <- max(sizes)
      cutoff_kmeans <- significant_numbers  # Set cutoff as highest abundance
      cat("Processing", basename(input_file), "- Kept only highest abundance sequence:", significant_numbers, "\n")
      # Store results in list
      state <- list(
        Filename = basename(input_file),
        Total_Sequences = length(sizes),
        Cutoff_KMeans = cutoff_kmeans,
        Significant_Sequences = 1,  # Only 1 sequence kept
        Significant_Numbers = significant_numbers
      )
      # Keep only highest size in merged_input
      merged_input <- merged_input[grep(paste0(";size=", significant_numbers), merged_input)]
      # Split at the "&" and unlist into a single column
      merged_input <- unlist(strsplit(merged_input, "&"))
      # Convert to a dataframe
      merged_input <- data.frame(Value = merged_input, stringsAsFactors = FALSE)
      # Save new fasta
      write.table(merged_input, file=output_filename, quote=F, row.names = F, col.names = F)
      # Save state to tsv
      write.table(as.data.frame(state), file=output_reportname, quote=F, row.names = F, col.names = T, sep="\t")
      
    }
  } else {
    # --- K-MEANS CLUSTERING METHOD (Using 7 Clusters) ---
    set.seed(42)  # For reproducibility
    kmeans_result <- kmeans(sizes, centers = 7)  # Use 7 clusters
    # Assign values to their respective clusters
    cluster_labels <- kmeans_result$cluster
    # Compute mean abundance for each cluster
    cluster_means <- tapply(sizes, cluster_labels, mean)
    # Rank clusters from lowest to highest mean
    ranked_clusters <- order(cluster_means)
    # Condition for datasets with low-abundance swarm clusters
    if (max(sizes) <= 150) {
      # Consider only the top 2 clusters
      top_clusters <- ranked_clusters[(length(ranked_clusters) - 1):length(ranked_clusters)]
    } else {
      # Consider the top 3 clusters
      top_clusters <- ranked_clusters[(length(ranked_clusters) - 2):length(ranked_clusters)]
    }
    # Select numbers in these top clusters
    significant_numbers <- sizes[cluster_labels %in% top_clusters]
    # Define a **less strict** cutoff (minimum of the selected cluster)
    cutoff_kmeans <- ifelse(length(significant_numbers) > 0, min(significant_numbers), NA)
    # store results in list
    state <- list(
      Filename = basename(input_file),
      Total_Sequences = length(sizes),
      Cutoff_KMeans = cutoff_kmeans,
      Significant_Sequences = length(significant_numbers),
      Significant_Numbers = ifelse(length(significant_numbers) > 0, paste(significant_numbers, collapse = ", "), "None")
    )
    # Print summary in console
    cat("\nProcessed:", basename(input_file))
    cat("\nCutoff based on K-Means:", cutoff_kmeans)
    cat("\nTotal sequences:", length(sizes))
    cat("\nSignificant sequences:", length(significant_numbers), "\n")
    # Select sizes that are significantly higher than the cutoff
    significant_sizes <- sizes[which(sizes >= cutoff_kmeans)]
    # Keep only significant sizes in merged_input
    pattern = paste0(";size=", significant_sizes)
    merged_input <- merged_input[grepl(pattern=paste(pattern, collapse="|"), merged_input)]
    # Split at the "&" and unlist into a single column
    merged_input <- unlist(strsplit(merged_input, "&"))
    # Convert to a dataframe
    merged_input <- data.frame(Value = merged_input, stringsAsFactors = FALSE)
    # Save new fasta
    write.table(merged_input, file=output_filename, quote=F, row.names = F, col.names = F)
    # Save state to tsv
    write.table(as.data.frame(state), file=output_reportname, quote=F, row.names = F, col.names = T, sep="\t")
  }
} else {
  #write empty files if input is empty
  if (file.size(input_file) == 0) {
    cat("Skipping", basename(input_file), "- No valid abundance numbers found.\n")
    system(paste0("touch ", output_filename))
    system(paste0("touch ", output_reportname))
  }
}
