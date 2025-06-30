#!/usr/bin/env Rscript

# Function to install and load packages
install_and_load <- function(packages) {
  for(package in packages) {
    if(!require(package, character.only = TRUE, quietly = TRUE)) {
      print(paste0("Installing: ", package))
      install.packages(package, quiet = TRUE, repos = "https://cran.rstudio.com/")
      library(package, character.only = TRUE)
    }
  }
}

# Use the function
required_packages <- c("dplyr", "purrr", "readr", "stringr", "writexl", "tidyr")
install_and_load(required_packages)

# Function to convert GTF to BED12
gtf_to_bed12 <- function(gtf_file) {
  
  # Read GTF file
  gtf <- read.table(gtf_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("seqname", "source", "feature", "start", "end", 
                                  "score", "strand", "frame", "attribute"))
  
  extract_attr <- function(attribute, key) {
    str_extract(attribute, str_glue("(?<={key}\\s)[^;]+"))
  }
  
  gtf <- gtf %>%
    mutate(
      gene_id = extract_attr(attribute, "gene_id"),
      transcript_id = extract_attr(attribute, "transcript_id"),
      tpm = as.numeric(extract_attr(attribute, "TPM")) 
    )
  
  # Filter for only records that are part of a transcript
  gtf_transcripts <- gtf %>% filter(!is.na(transcript_id))
  
  # Calculate total TPM from transcript features only
  total_tpm <- gtf_transcripts %>%
    filter(feature == "transcript") %>%
    summarise(total = sum(tpm, na.rm = TRUE)) %>%
    pull(total)
  
  # First, get transcript-level information (name, TPM, etc.)
  transcript_info <- gtf_transcripts %>%
    filter(feature == "transcript") %>%
    select(seqname, strand, transcript_id, gene_id, tpm)
  
  # Second, process exons to get block info for each transcript
  exon_info <- gtf_transcripts %>%
    filter(feature == "exon") %>%
    arrange(transcript_id, start) %>% # Sort exons by position within each transcript
    group_by(transcript_id) %>%
    summarise(
      # BED is 0-based, GTF is 1-based
      chromStart = min(start) - 1, 
      chromEnd = max(end),
      blockCount = n(),
      # Calculate block sizes and starts relative to the transcript start
      blockSizes = paste(end - start + 1, collapse = ","),
      blockStarts = paste(start - (min(start)), collapse = ","),
      # Also calculate ExonCollapsed length for the name field
      ef = sum(end - start + 1)
    ) %>%
    ungroup()
  
  bed12_prelim <- inner_join(transcript_info, exon_info, by = "transcript_id") %>%
    # Handle cases where TPM might be NA
    mutate(
      rel = ifelse(is.na(tpm), 0, tpm / total_tpm),
      # Safely extract numbers for the name field
      gene_num = str_extract(gene_id, "\\d+"),
      transcript_num = str_extract(transcript_id, "\\d+\\.\\d+|\\d+"), # handles STRG.1.1 and STRG.1
      name = sprintf("(%s)G%s|TU%s(%.1f)(%d)(%.2f)", strand, gene_num, transcript_num, tpm, ef, rel),
      score = 1000,
      # itemRgb = case_when(rel == max(rel) ~ "255,0,0",
      #                     .default = "0"),
      itemRgb = case_when(strand == "+" ~ "0,0,255",
                          .default = "255,0,0"),
      thickStart = chromStart,
      thickEnd = chromStart
    )
  
  bed12 <- bed12_prelim %>% 
    select(
      chrom = seqname,
      chromStart,
      chromEnd,
      name,
      score,
      strand,
      thickStart,
      thickEnd,
      itemRgb,
      blockCount,
      blockSizes,
      blockStarts
    ) %>% arrange(chrom, chromStart)
  
  bed12df <- bed12_prelim %>% 
    select(
      tpm,
      relative_fraction = rel,
      effective_length = ef,
      chrom = seqname,
      chromStart,
      chromEnd,
      strand,
      name
    ) %>% arrange(chrom, chromStart)
  
  
  
  list(bed12 = bed12, bed12df = bed12df)

}

main <- function(){
  args = commandArgs(trailingOnly=TRUE)
  input_gtf <- args[1]
  output_bed12 <- args[2]
  output_xlsx <- paste0(tools::file_path_sans_ext(output_bed12),".StringTie.tsv")
  
  bed12_result <- gtf_to_bed12(input_gtf)
  
  # write.table(bed12_result$bed12, output_bed12, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write_tsv(bed12_result$bed12, path = output_bed12, col_names = F)
  write_tsv(bed12_result$bed12df, path = output_xlsx, col_names = T)
}

main()
