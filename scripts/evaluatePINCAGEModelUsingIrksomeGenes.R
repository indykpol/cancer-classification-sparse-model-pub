#!/usr/bin/env Rscript
library(optparse)
# Command Line Arguments
option_list <- list(
  make_option(c("-d","--data"), type='character', help = 'cached data to use', default = "data/RData/data_LUSC_all_lazyReady.R_CACHE"),
  make_option(c("-r","--results"), type='character', help = 'model results to use', default = "data/latestResults/LUSC_Ts_vs_ANs_4foldCV_withFixedPrior.RData"),
  make_option(c("-f","--folds"), type='integer', help = 'number of folds to split the dataset', default = 4),
  make_option(c("-n","--name"), type='character', help = 'name of the analysis performed', default = "MLs_Irksome_LUSC_TvsAN_4foldCV"),
  make_option("--n_top", type='integer', help = 'number of top ranks/genes to train for and combine', default = 100),
  make_option(c("-p","--parallel"), type='integer', help = 'if desired, provide a number of cores to parallelize the execution of the analysis', default = 1)
)
opts <- opt <- parse_args(OptionParser(option_list = option_list))
#opt <- list(data = "data/RData/data_LUSC_all_lazyReady.R_CACHE", results="data/latestResults/LUSC_Ts_vs_ANs_4foldCV_withFixedPrior.RData",  name = "MLs_Irksome_LUSC_TvsAN_4foldCV", folds = 4, n_top = 100, parallel=16)

library(caret)
library(pROC)
library(dplyr)
library(parallel)
source("R/utilities_new.R")