#!/usr/bin/env Rscript

# libraries
library(magrittr)
library(flowCore)

# load fcs data
fcs <-
  list.files("data/flow_cytometry", pattern = ".fcs", full.names = TRUE) %>%
  {setNames(., basename(.))} %>%
  purrr::map(read.FCS)

# check area
fcs %>%
  purrr::map(function(fcs_i) {
    fcs_a <- exprs(fcs_i)[, "FSC04-A"]
    hist(fcs_a, breaks = 100, main = "Cell Size Distribution", xlab = "FSC-A")
  })
