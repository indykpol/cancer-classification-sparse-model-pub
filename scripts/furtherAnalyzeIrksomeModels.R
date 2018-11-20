#!/usr/bin/env Rscript
library(optparse)
option_list <- list(
  make_option(c("-d","--path"), type='character', help = 'optional directory to look into, defaults to current working directory from where the script is called from', default = NULL),
  make_option(c("-g","--gene_names"), type='character', help = 'pointer to the RData file with geneNames', default = '/project/iCancerGenomics/irksome-fibula/data/latestResults/geneNames_BRCA_progression_filtered.RData'),
  make_option(c("-p","--parallel"), type='integer', help = 'if desired, provide a number of cores to parallelize the execution of the analysis', default = 64),
  make_option(c("-t", "--top"), type='integer', default = 100, help = 'number of top genes/ranks to output in the table'),
  make_option(c("-n","--name"), type='character', default = NULL, help = 'desired name of the output table'),
  make_option("--g1_length", type='integer', default = NULL, help = 'the number of group 1 samples in the input table'),
  make_option("--g2_length", type='integer', default = NULL, help = 'the number of group 2 samples in the input table'),
  make_option(c("-r","--ranks"), type='character', help = "pointer to the RData file with the 'ranks' object with ranking across all folds", default = NA)
)
opt <- parse_args(OptionParser(option_list = option_list))
# opt <- list(parallel=64, top=100, name="BRCA_progressing_IrksomeGenes", g1_length=38, g2_length=43, gene_names="/project/iCancerGenomics/irksome-fibula/data/latestResults/geneNames_BRCA_progression_filtered.RData", name="BRCA_progressing", ranks="/project/iCancerGenomics/irksome-fibula/data/latestResults/ranks_BRCA_progressing_10foldCV.RData")
library(pROC)
library(parallel)
library(dplyr)
if(is.null(opt$path)) opt$path <- getwd()
load(opt$gene_names)

sum_asLogs <- function (x) {
	if (x[1] >= x[2]) x[1] + log(1+ exp(x[2]-x[1])) else x[2] + log(1+ exp(x[1]-x[2]))
}

results_list <- dir(path=opt$path, pattern=".result", full.names=TRUE)
cat(paste("Found", length(results_list), "models in", opt$path, "\nBeginning reading using", opt$parallel,"thread(s)\n"))

ptm <- proc.time()[3]
out <- mclapply(1:length(results_list), mc.cores = opt$parallel, FUN = function(i) {
	Zs <- read.table(results_list[i], nrows=1)
	Likelihoods <- read.table(results_list[i], skip=1, header=TRUE)
	list(Zs, Likelihoods)
})

cat(paste("Done reading models in ", sprintf("%.2f", (proc.time()[3]-ptm)/60)," minutes\n", sep=""))

errors <- grep(pattern="Error", out)
while (length(errors) > 0) {
	cat(paste("Encountered",  length(grep(pattern="Error", out)), "reading error(s). Attempting serial reading.\n"))
	ptm <- proc.time()[3]
	for (i in 1:length(errors)) {
		out[[errors[i]]] <- tryCatch(list(read.table(results_list[i], nrows=1), read.table(results_list[i], skip=1, header=TRUE)), error = function(e) return("Error"))
	}
	cat(paste("Done reading errored models in", sprintf("%.2f", (proc.time()[3]-ptm)/60), "minutes\n"))
	errors <- grep(pattern="Error", out)
}


# Make Z scores matrix
Zs <- matrix(nrow=length(results_list), ncol=length(out[[1]][[1]]))
for (i in 1:length(results_list)) Zs[i,] <- unlist(out[[i]][[1]])
if (is.na(opt$ranks)) ranks <- apply(-Zs, 2, rank) else load(opt$ranks)
# Use split to recreate the CV-sets
set.seed(1) # Use the same sets accross genes
G1_cv_sets <- 1:opt$g1_length %>% sample %>% split(1:ncol(Zs))
G2_cv_sets <- (opt$g1_length+1):(opt$g1_length+opt$g2_length) %>% sample %>% split(1:ncol(Zs))

# setup out table and likelihood combination table
auc_top <- NULL
auc_top_combined <- NULL
Likelihoods_combined <- matrix(nrow=opt$g1_length+opt$g2_length, ncol=2, data=0)
for (i in 1:opt$top) {
	Likelihoods_single <- matrix(nrow=opt$g1_length+opt$g2_length, ncol=2, data=0)
	for (j in 1:ncol(Zs)) {
		Likelihoods_single[c(unlist(G1_cv_sets[[j]]), unlist(G2_cv_sets[[j]])),] <- as.matrix(out[[which(ranks[,j]==i)]][[2]][c(unlist(G1_cv_sets[[j]]), unlist(G2_cv_sets[[j]])), 2:3])
	}
	Likelihoods_combined <- Likelihoods_combined+Likelihoods_single
	
	response <-  c(rep(1, opt$g1_length), rep(0, opt$g2_length))
	auc_top[[i]] <- tryCatch(auc(response = response, predictor = Likelihoods_single[,1]/apply(Likelihoods_single, 1, sum_asLogs))[1], error = function(e) return(NA))
  auc_top_combined[[i]] <- tryCatch(auc(response = response, predictor = Likelihoods_combined[,1]/apply(Likelihoods_combined, 1, sum_asLogs))[1], error = function(e) return(NA))
}

IDs <- as.numeric(unlist(lapply(1:length(results_list), FUN = function(x) strsplit(results_list[x],  ".result")[1])))
top_genes <- ranks[(order(apply(ranks, 1, mean))),]
top_genes <- head(data.frame(geneID = geneNames[IDs[order(apply(ranks, 1, mean))]], top_genes), n=opt$top)

out_table <- data.frame(top= 1:length(auc_top), auc_top, auc_top_combined, matrix(nrow = length(auc_top), ncol = ncol(Zs), data = unlist(lapply(1:opt$top, FUN=function(x) unlist(lapply(1:ncol(Zs), FUN=function(y) geneNames[IDs[which(ranks[,y]==x)]])))), byrow=TRUE))
write.table(out_table, file=paste("aucTable_PINCAGE_", opt$name, ".tab", sep=""), quote=FALSE, sep = "\t", row.names = FALSE)
write.table(top_genes, file=paste("aucTable_PINCAGE_", opt$name, ".tab", sep=""), quote=FALSE, sep = "\t", row.names = FALSE, append=TRUE)
