#!/usr/bin/env Rscript

# This script takes a fasta file as input and returns a 
# fasta file with the first 20 sequences of the input 

# set arguments 
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# get args value (input file name and output dir)
input_file <- args[1]
output_dir <- args[2]

# extract sample name and set output filename and suffix
input_file_suffix <- ".swarm.seeds.fasta"
sample <- gsub(input_file_suffix, "", basename(input_file))
output_file_suffix <- ".top20.swarm.seeds.fasta"
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
  # Keep 20 first sequences in merged_input
  merged_input <- merged_input[1:20]
  # Revome NAs in case there are less than 20 seqs
  merged_input <- merged_input[!is.na(merged_input)]
  # Split at the "&" and unlist into a single column
  merged_input <- unlist(strsplit(merged_input, "&"))
  # Convert to a dataframe
  merged_input <- data.frame(Value = merged_input, stringsAsFactors = FALSE)
  # Save new fasta
  write.table(merged_input, file=output_filename, quote=F, row.names = F, col.names = F)
} else {
  # if input file was empty, create an empty output file
  system(paste0("touch ", output_filename))
}


