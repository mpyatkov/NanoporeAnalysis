#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(stringr)
library(readr)

args = commandArgs(trailingOnly=TRUE)
fgb_featurecounts_tsv <- args[1]

## If gene mentioned more then one time, take full span as min(start) - max(end), counts will be the same
## If two regions has the same start but different ends aggregate them as max(end), counts will be a sum
read_tsv(fgb_featurecounts_tsv, comment = "#", show_col_types = F) %>% 
  rename(counts = 7) %>% 
  filter(counts >= 1) %>% 
  tidyr::separate_longer_delim(c(Chr,Start,End,Strand), delim = ";") %>% 
  summarize(Chr = first(Chr), Start = min(Start), End = max(End), counts = first(counts), .by = "Geneid") %>% 
  summarize(End = max(End), counts = sum(counts), .by = c("Chr", "Start")) %>% 
  mutate(coords=str_glue("{Chr}:{Start}-{End}")) %>% 
  select(coords,counts) %>% 
  write_csv("coords_counts.csv", col_names = F)
