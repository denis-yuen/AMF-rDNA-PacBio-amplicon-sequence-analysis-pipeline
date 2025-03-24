#!/usr/bin/env Rscript

# This script takes three inputs:
#   - a fasta file (clusters from swarm)
#   - a blast result file (*.pSSU_ITS_pLSU.blasted.txt)
#   - an output folder name
# It returns a fasta file from which the non-mucoromycota clusters
# have been removed.

# get arguments
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# set args values (input file name, blast file and output dir)
input_file <- args[1]
blast_file <- args[2]
output_dir <- args[3]

# extract sample name and set output suffix and file name.
input_file_suffix <- ".swarm.seeds.fasta"
sample <- gsub(input_file_suffix, "", basename(input_file))
output_file_suffix <- ".tax_filtered.swarm.seeds.fasta"
output_filename <- paste0(output_dir,"/", sample, output_file_suffix)

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
  # read blast File
  blast <- read.csv(blast_file, sep="\t")  
  # keep blast results of current sample only
  blast <- blast[grepl(sample, blast$sample_file),] 
  # keep non-mucoromycota
  blast <- blast[blast$phylum != "Mucoromycota",]
  # if all top clusters are non-mucoromycota, return empty file
  if (nrow(blast) == 20) {
    system(paste0("touch ", output_filename)) 
  # if all top clusters are mucoromycota, keep all sequences in merged_input
  } else if (nrow(blast) == 0) {
    # Split at the "&" and unlist into a single column
    merged_input <- unlist(strsplit(merged_input, "&"))  
    # Convert to a dataframe
    merged_input <- data.frame(Value = merged_input, stringsAsFactors = FALSE)  
    # Save new fasta
    write.table(merged_input, file=output_filename, quote=F, row.names = F, col.names = F)
  } else {
    # remove non-mucoromycota from merged_input
    merged_input <- merged_input[!grepl(pattern=paste(blast$query, collapse="|"), merged_input)]  
    # Split at the "&" and unlist into a single column
    merged_input <- unlist(strsplit(merged_input, "&"))
    # Convert to a dataframe
    merged_input <- data.frame(Value = merged_input, stringsAsFactors = FALSE)
    # Save new fasta
    write.table(merged_input, file=output_filename, quote=F, row.names = F, col.names = F)
  }
} else {
  # if input_file is empty, write an empty output file
  system(paste0("touch ", output_filename))
}
