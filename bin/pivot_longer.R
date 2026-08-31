#! /usr/bin/env Rscript

library(tidyverse)
library(magrittr)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("Usage: Rscript pivot_longer.R <sample_variants.tsv>")
}

variants_file <- args[1]

data <- read_delim(variants_file, delim="\t", na = character()) %>%
  mutate(
    IsPseudo = case_when(IsPseudo==1 ~ TRUE, IsPseudo==0 ~ FALSE, TRUE ~ NA),
    IsGenic = case_when(IsGenic==1 ~ TRUE, IsGenic==0 ~ FALSE, TRUE ~ NA),
    IsTransition = case_when(IsTransition==1 ~ TRUE, IsTransition==0 ~ FALSE, TRUE ~ NA),
    IsSynonymous = case_when(IsSynonymous==1 ~ TRUE, IsSynonymous==0 ~ FALSE, TRUE ~ NA)
  )

id_col <- names(data)[1:19]
val_col <- names(data)[20:ncol(data)]

cdata <- pivot_longer(data, cols=val_col, names_to="Sample", values_to = "Percentage") %>%
  relocate(c("Sample","Percentage")) %>%
  subset(Percentage != ".")

writexl::write_xlsx(cdata, path="samples_variant_results_longer.xlsx")
