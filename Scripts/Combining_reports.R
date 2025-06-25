library(stringr)
library(writexl)
library(readr)
library(purrr)
library(dplyr)


tsv_files <- list.files(pattern = "tsv", recursive = T) %>% 
  keep(~str_detect(.x, "Filtered|StringTie.tsv")) %>% 
  tibble(filename = .) %>% 
  mutate(shortcut = case_when(str_detect(filename,"unspliced.StringTie.tsv") ~ "StringTie_unspliced",
                              str_detect(filename,"spliced.StringTie.tsv") ~ "StringTie_spliced",
                              str_detect(filename, "ExonCollapsed.Spliced") ~ "FeatCounts.ExonColl.Spliced",
                              str_detect(filename, "ExonCollapsed.Unspliced") ~ "FeatCounts.ExonColl.Unspliced",
                              str_detect(filename, "FullGeneBody.Spliced") ~ "FeatCounts.FullGB.Spliced",
                              .default = "FeatCounts.FullGB.Unspliced"))


summary_files <- list.files(pattern = "summary", recursive = T) %>% 
  tibble(filename = .) %>% 
  mutate(shortcut = case_when(str_detect(filename,"ExonCollapsed.Spliced") ~ "Sum.ExonColl.Spliced",
                              str_detect(filename,"ExonCollapsed.Unspliced") ~ "Sum.ExonColl.Unspliced",
                              str_detect(filename,"FullGeneBody.Spliced") ~ "Sum.FullGB.Spliced",
                              .default = "Sum.FullGB.Unspliced"))

files <- bind_rows(tsv_files, summary_files)

map2(files$filename, files$shortcut, \(filename, shortcut) {
  if (str_detect(shortcut, "StringTie|Sum")) {
    read_tsv(filename, col_names = T, show_col_types = FALSE)
  } else {
    read_tsv(filename, col_names = T, show_col_types = FALSE) %>% 
      select(GeneID = 1, Counts = 2)
  }
}) %>% 
  set_names(files$shortcut) %>% 
  writexl::write_xlsx("Combined_Report.xlsx")
