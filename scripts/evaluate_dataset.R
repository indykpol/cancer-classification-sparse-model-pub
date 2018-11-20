#################################################
#  Evaluate a whole dataset, combining the p-vals
# from expression and methylation|expression 
# models, and generate prioritized list of genes
#################################################

library(optparse)
# Command Line Arguments
option_list <- list(
  make_option(c("-d","--data"), type='character', help = 'cached data to use', default = "data/RData/data_LUSC_all_lazyReady.R_CACHE"),
  make_option(c("-s","--cases"), type='character', help = 'comma-separated list of names of sample groups present in the RData file to test against', default = 'ANs'),
  make_option(c("-c","--controls"), type='character', help = 'optional comma-separated list of names of the other sample groups present in the RData', default = "Ts"),
  make_option(c("-f","--folds"), type='integer', default = "4", help = 'number of folds to split the dataset'),
  make_option(c("-r","--fraction"), type='integer', default = "NA", help = 'which fraction of the dataset to hold out and learn prior from'),
  make_option(c("-n","--name"), type='character', help = 'name of the analysis performed', default = "LUSC_Ts_vs_ANs_4foldCV"),
  make_option(c("-p","--parallel"), type='integer', help = 'if desired, provide a number of cores to parallelize the execution of the analysis', default = 1)
)
#opt <- parse_args(OptionParser(option_list = option_list))
opt <- list(data = "data/RData/data_PRAD_all_lazyReady.R_CACHE",  folds = 7, fraction = NA, cases = "ANs", controls = "Ts", name = "PRAD_Ts_vs_ANs_7foldCV_uneqalVariances", parallel = 16)

library(parallel)
library(methods)
library(IrksomeModels)
library(dplyr)
library(pROC)
source("R/utilities_new.R")

grouping <- read_grouping(cache = opt$data, group1 = opt$cases, group2 = opt$controls)
samplenames <- unlist(grouping)
counts <- read_element(cache = opt$data, element_name = "counts")
lib_sizes <- read_element(cache = opt$data, element_name = "lib_sizes")
# normalize expression
expr_logUQ <- normalizeExpr(counts, method = "UQ")
# find consistently lowly expressed genes and filter them out
weights <- exprWeight(expr_logUQ, plot = FALSE)
expr_logUQ <- expr_logUQ[weights < 0.5,]

if (!is.na(opt$fraction)) { # if we set aside a fraction of the dataset for validation purposes
  division <- split(1:length(samplenames), 1:opt$fraction)
  samplenames_l <- samplenames[-division[[1]]]
  
  log_pvals_expr <- apply(expr_logUQ[,samplenames_l], 1, function (x) testDE(x, samplenames_l %in% grouping[[1]], k_win = 2.5))
  log_pvals_meth <- vector(mode="numeric", length=nrow(expr_logUQ))
  for (i in 1:nrow(expr_logUQ)) {
    log_pvals_meth[i] <- testDM(meth = read_methylation_matrix_gene(cache = opt$data, genename = rownames(expr_logUQ)[i]), expr = expr_logUQ[i, samplenames_l], grouping = samplenames_l %in% grouping[[1]])
  }
  log_pvals_fishers <- apply(cbind(exp(log_pvals_expr), exp(log_pvals_meth)), 1, fishersMethod, log.p = TRUE)
} else { # if we run analysis on the whole set
  # Use split to generate CV-sets
  set.seed(1) # Use the same sets accross genes
  g1 <- which(samplenames %in% grouping[[1]])
  g2 <- which(samplenames %in% grouping[[2]])
  g1_cv_sets <- g1 %>% sample %>% split(1:opt$folds)
  g2_cv_sets <- g2 %>% sample %>% split(1:opt$folds)
  
  log_pvals_meth <- auc_cv <- matrix(nrow = nrow(expr_logUQ), ncol = opt$folds)
  fm <- fe <- log_pvals_fishers <- log_pvals_expr_scale <- log_pvals_expr_welch <- log_pvals_expr <- pvals_folds <- joint_logliks_test_list <- meth_logliks_test_list <- expr_logliks_test_list <- ranks <- NULL
  
  # calculate methylation prior
  prior_meth <- methPrior(cache = opt$data, K = nrow(expr_logUQ), grouping = rep(c(TRUE, FALSE), times=c(length(g1), length(g2))), exprMat = expr_logUQ, parallel = opt$parallel)
  # calculate expression prior
  prior_expr <- exprPrior(expr_logUQ, plot = FALSE)
  
  # Run CV
  for(i in 1:opt$folds) {
    train_set_g1 <- setdiff(g1, g1_cv_sets[[i]])
    train_set_g2 <- setdiff(g2, g2_cv_sets[[i]])
    test_set_g1 <- g1_cv_sets[[i]]
    test_set_g2 <- g2_cv_sets[[i]]
    train_labels <- c(rep(TRUE, length(train_set_g1)), rep(FALSE, length(train_set_g2)) )
    test_labels <- c(rep(0, length(test_set_g1)), rep(1, length(test_set_g2)) )
    
    # test the diff. expression
	tmp <- apply(expr_logUQ[,c(train_set_g1, train_set_g2)], 1, function (x) testDE(x, train_labels, k_win = 2.5, prior = prior_expr))
	log_pvals_expr[[i]] <- unlist(lapply(1:length(tmp), FUN = function(gene) tmp[[gene]]$p_location))
	log_pvals_expr_welch[[i]] <- unlist(lapply(1:length(tmp), FUN = function(gene) tmp[[gene]]$p_location_welch))
	log_pvals_expr_scale[[i]] <- unlist(lapply(1:length(tmp), FUN = function(gene) tmp[[gene]]$p_scale))
    
    # Iterate over all considered gene names
    auc_fold <- NULL
    log_pvals_meth_fold <- NULL
    fe_fold <- NULL
    fm_fold <- NULL
    expr_logliks_test_fold <- NULL
    meth_logliks_test_fold <- NULL
    joint_logliks_test_fold <- NULL
	out <- mclapply(1:nrow(expr_logUQ), mc.cores = opt$parallel, FUN = function(j) {
			# read respective methylation data
			meth <- read_methylation_matrix_gene(cache = opt$data, genename = rownames(expr_logUQ)[j])
			# winsorize it
			meth <- apply(meth, 2, methylationWinsor)
			
			# test for the differential methylation | expression
			log_pvals_meth_fold[[j]] <- testDM(meth = meth[c(train_set_g1, train_set_g2),], expr = expr_logUQ[j, c(train_set_g1, train_set_g2)], grouping = train_labels, prior = prior_meth)$p_location
			test <- testDM(meth = meth[c(train_set_g1, train_set_g2),], expr = expr_logUQ[j, c(train_set_g1, train_set_g2)], grouping = train_labels, prior = prior_meth)$p_bartlett
			if (!is.na(test) & test <= log(0.05)) same_cov <- FALSE else same_cov <- TRUE
			if (log_pvals_expr_scale[[i]][j] <= log(0.05)) same_sd <- FALSE else same_sd <- TRUE
			
			# Estimate parameters and calculate likelihood
			fe_fold[[j]] <- fitExpr(expr = expr_logUQ[j, c(train_set_g1, train_set_g2)], grouping = train_labels, k_win = 2.5, prior = prior_expr)
			fm_fold[[j]] <- fitMeth(meth = meth[c(train_set_g1, train_set_g2),], expr = expr_logUQ[j, c(train_set_g1, train_set_g2)], grouping = train_labels, prior = prior_meth)
			
			# Calculate expression and methylation likelihoods separately
			expr_logliks_test_fold[[j]] <- expr_logliks_test <- exprLogLik(expr = expr_logUQ[j, c(test_set_g1, test_set_g2)], param = fe_fold[[j]], k_win = 2.5, same_sd = same_sd)
			meth_logliks_test_fold[[j]] <- meth_logliks_test <- methLogLik(meth = meth[c(test_set_g1, test_set_g2),], expr = expr_logUQ[j, c(test_set_g1, test_set_g2)], param = fm_fold[[j]]$param, same_cov = same_cov)
			# Calculate the joint likelihood
			joint_logliks_test_fold[[j]] <- logliks_test <- expr_logliks_test + meth_logliks_test
			
			# Calculate the test AUC
			auc_fold[[j]] <- tryCatch(auc(response= test_labels, predictor = logliks_test)[1], error = function(e) return(NA))
			list(log_pvals_meth_fold[[j]], fe_fold[[j]], fm_fold[[j]], expr_logliks_test_fold[[j]], meth_logliks_test_fold[[j]], joint_logliks_test_fold[[j]], auc_fold[[j]])
		})
	for (k in 1:length(out)) { # transfer the results from the parallel execution
		log_pvals_meth_fold[[k]] <- out[[k]][[1]]
		fe_fold[[k]] <- out[[k]][[2]]
		fm_fold[[k]] <- out[[k]][[3]]
		expr_logliks_test_fold[[k]] <- out[[k]][[4]]
		meth_logliks_test_fold[[k]] <- out[[k]][[5]]
		joint_logliks_test_fold[[k]] <- out[[k]][[6]]
		auc_fold[[k]] <- out[[k]][[7]]
	}
	rm(out)
	
	fe[[i]] <- fe_fold
	fm[[i]] <- fm_fold
	expr_logliks_test_list[[i]] <- expr_logliks_test_fold
	meth_logliks_test_list[[i]] <- meth_logliks_test_fold
	joint_logliks_test_list[[i]] <- joint_logliks_test_fold
	auc_cv[,i] <- unlist(auc_fold)
	log_pvals_meth[,i] <- unlist(log_pvals_meth_fold)
	# Combine evidence using Fisher's method
	log_pvals_fishers[[i]] <- apply(cbind(log_pvals_expr_scale[[i]], log_pvals_expr_welch[[i]], log_pvals_expr[[i]], log_pvals_meth[,i]), 1, FUN = function(x) if (x[1] <= log(0.05)) fishersMethod(x[c(2,4)], log.p = TRUE, logs = TRUE) else fishersMethod(x[c(3,4)], log.p = TRUE, logs = TRUE)) # variance logic
	pvals_folds[[i]] <- data.frame(genename = rownames(expr_logUQ), log_pvals_expr = log_pvals_expr[[i]], log_pvals_expr_welch = log_pvals_expr_welch[[i]], log_pvals_expr_scale = log_pvals_expr_scale[[i]], log_pvals_meth=log_pvals_meth[,i], log_pvals_fishers = log_pvals_fishers[[i]])
  }
}
save.image(paste("data/latestResults/", opt$name, ".RData", sep=""))

log_pvals_fishers <- matrix(nrow=nrow(expr_logUQ), ncol=opt$folds, data=0)
for (i in 1:opt$folds) log_pvals_fishers[,i] <- pvals_folds[[i]]$log_pvals_fishers
ranks <- apply(log_pvals_fishers, 2, rank)
ranks <- do.call(cbind, lapply(1:opt$folds, FUN=function(x) rank(pvals_folds[[x]]$log_pvals_fishers)))
top_genes <- ranks[(order(apply(ranks, 1, mean))),]
top_genes <- data.frame(geneID = rownames(expr_logUQ)[order(apply(ranks, 1, mean))], top_genes)
# combine evidence
auc_top <- NULL
auc_top_combined <- NULL
predictor_top_combined <- rep(0, length(g1)+length(g2))
for (i in 1:100) {
	min_ranks <- apply(ranks, 2, min)
	if (sum(min_ranks) > opt$folds) {
		cur_rank <- NULL
		for (j in 1:opt$folds) {
			if (i < min_ranks[j]*2) which <- which(ranks[,j]==min_ranks[j]) else which <- which(ranks[,j]==i)
			if (length(which) > 1) which <- which[which(rank(pvals_folds[[j]][which,"log_pvals_expr"]) == i)]
			cur_rank[[j]] <- which
		}
		predictor <- unlist(lapply(1:opt$folds, FUN=function(x) joint_logliks_test_list[[x]][[cur_rank[x]]]))
	} else predictor <- unlist(lapply(1:opt$folds, FUN=function(x) joint_logliks_test_list[[x]][[which(ranks[,x]==i)]]))
  
  predictor_top_combined <- predictor_top_combined + predictor
  response <- unlist(lapply(1:opt$folds, FUN=function(x) c(rep(1, length(g1_cv_sets[[x]])), rep(0, length(g2_cv_sets[[x]])))))
  auc_top[[i]] <- tryCatch(auc(response = response, predictor = predictor)[1], error = function(e) return(NA))
  auc_top_combined[[i]] <- tryCatch(auc(response = response, predictor = predictor_top_combined)[1], error = function(e) return(NA))
}

out_table <- data.frame(top= 1:length(auc_top), auc_top, auc_top_combined, matrix(nrow = length(auc_top), ncol = opt$folds, data = unlist(lapply(1:length(auc_top), FUN=function(x) unlist(lapply(1:opt$folds, FUN=function(y) rownames(expr_logUQ)[which(ranks[,y]==x)])))), byrow=TRUE))

# Write outputs
save.image(paste("data/latestResults/", opt$name, ".RData", sep=""))
write.table(out_table, file=paste("data/latestResults/", opt$name, ".tab", sep=""), quote=FALSE, sep = "\t", row.names = FALSE)
write.table(top_genes, file=paste("data/latestResults/", opt$name, ".tab", sep=""), quote=FALSE, sep = "\t", row.names = FALSE, append=TRUE)