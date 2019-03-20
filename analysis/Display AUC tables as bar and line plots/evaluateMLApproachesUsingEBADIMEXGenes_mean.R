#!/usr/bin/env Rscript
library(optparse)
# Command Line Arguments
option_list <- list(
  make_option(c("-d","--data"), type='character', help = 'cached data to use', default = "data/RData/data_LUSC_all_lazyReady.R_CACHE"),
  make_option(c("-r","--results"), type='character', help = 'model results to use', default = "data/latestResults/top_100_each_seed.RData"),
  make_option(c("-c","--cv_runs"), type='integer', help = 'number of cross-validations runs to perform', default = 50),
  make_option(c("-n","--name"), type='character', help = 'name of the analysis performed', default = "MLs_EBADIMEX_50randomSamples_BRCA_progressing"),
  make_option("--n_top", type='integer', help = 'number of top ranks/genes to train for and combine', default = 100),
  make_option(c("-p","--parallel"), type='integer', help = 'if desired, provide a number of cores to parallelize the execution of the analysis', default = 1)
)

#opt <- parse_args(OptionParser(option_list = option_list))
opt <- list(data = "data/RData/data_BRCA_progressing_lazyReady.R_CACHE", results="data/latestResults/top_100_each_seed.RData",  name = "MLs_EBADIMEX_50randomSamples_BRCA_progressing_mean_UQ", cv_runs = 50, n_top = 100, parallel=16)

library(caret)
library(pROC)
library(dplyr)
library(parallel)
source("R/utilities_new.R")

load(opt$results)
load("UQs_BRCA_progressing.RData")
samplenames <- unlist(read_grouping(cache_name = opt$data, group1="progressed", group2="nonProgressed"))
indiv_predictions_EBADIMEX_cv <- matrix(ncol=length(samplenames), nrow=opt$n_top)
colnames(indiv_predictions_EBADIMEX_cv) <- samplenames
indiv_predictions_EBADIMEX_cv_SVMPol <- indiv_predictions_EBADIMEX_cv_SVMLin <- indiv_predictions_EBADIMEX_cv_LR <- indiv_predictions_EBADIMEX_cv_BLR <- indiv_predictions_EBADIMEX_cv_RF <- indiv_predictions_EBADIMEX_cv

AUCs <- NULL
for (cv_run in 1:50) {
	
	set.seed(cv_run)
	train_set_g1 <- sample(1:38, 10)
	train_set_g2 <- 38+sample(1:43, 15)
	test_set_g1 <- (1:38)[-train_set_g1]
	test_set_g2 <- 38+(1:43)[-(train_set_g2-38)]
	
	train_labels <- as.factor(c(rep("g1", length(train_set_g1)), rep("g2", length(train_set_g2))))
	
	current <- top_100_each_seed[[cv_run]] #  list of top 100 genes for the current c-v run
	if (opt$parallel > 1) {
		out <- mclapply(1:length(current), mc.cores = opt$parallel, FUN = function(gene) {
			genedata <- as.data.frame(read_genedata(cache = opt$data, genename = current[gene]))
			genedata$expr_UQ <- genedata$read_count / UQs
			genedata_training <- genedata %>% 
				dplyr::select(starts_with("pr_"), starts_with("gb_"), matches("expr_UQ")) %>%
				slice(c(train_set_g1, train_set_g2))
			genedata_testing <- genedata %>% 
				dplyr::select(starts_with("pr_"), starts_with("gb_"), matches("expr_UQ")) %>%
				slice(c(test_set_g1, test_set_g2))
			
			preProcValues <- preProcess(genedata_training, method = c("center", "scale"))
			genedata_training_CS <- predict(preProcValues, genedata_training)
			genedata_testing_CS <- predict(preProcValues, genedata_testing)
			set.seed(1)
			# Logistic Regression
			model <- train(x=genedata_training, y=train_labels, method="glm")
			indiv_predictions_EBADIMEX_cv_LR[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing, type="prob")[,1]
			# Bayesian Logistic Regression
			model <- train(x=genedata_training, y=train_labels, method="bayesglm")
			indiv_predictions_EBADIMEX_cv_BLR[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing, type="prob")[,1]
			# Random Forest
			model <- train(x=genedata_training, y=train_labels, method="rf")
			indiv_predictions_EBADIMEX_cv_RF[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing, type="prob")[,1]
			# SVM with linear kernel
			model <- train(x=genedata_training_CS, y=train_labels, method="svmLinear2", trControl=trainControl(classProbs = TRUE))
			indiv_predictions_EBADIMEX_cv_SVMLin[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing_CS, type="prob")[,1]
			# SVM with polynomial kernel
			model <- train(x=genedata_training_CS, y=train_labels, method="svmPoly", trControl=trainControl(classProbs = TRUE))
			indiv_predictions_EBADIMEX_cv_SVMPol[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing_CS, type="prob")[,1]
			list(indiv_predictions_EBADIMEX_cv_LR[gene, c(test_set_g1, test_set_g2)], indiv_predictions_EBADIMEX_cv_RF[gene, c(test_set_g1, test_set_g2)], indiv_predictions_EBADIMEX_cv_BLR[gene, c(test_set_g1, test_set_g2)], indiv_predictions_EBADIMEX_cv_SVMLin[gene, c(test_set_g1, test_set_g2)], indiv_predictions_EBADIMEX_cv_SVMPol[gene, c(test_set_g1, test_set_g2)])
		})
		
		for (gene in 1:length(out)) { # transfer the results from the parallel execution
			indiv_predictions_EBADIMEX_cv_LR[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[1]]
			indiv_predictions_EBADIMEX_cv_RF[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[2]]
			indiv_predictions_EBADIMEX_cv_BLR[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[3]]
			indiv_predictions_EBADIMEX_cv_SVMLin[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[4]]
			indiv_predictions_EBADIMEX_cv_SVMPol[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[5]]
		}
		rm(out)
	} else stop("Serial execution currently unavailable, use --parallel >1")
	
	
	auc_top_SVMLin <- auc_top_SVMPol <- auc_top_RF <- auc_top_LR <- auc_top_BLR <- NULL
	auc_top_combined_SVMLin <- auc_top_combined_SVMPol <- auc_top_combined_RF <- auc_top_combined_LR <- auc_top_combined_BLR <- NULL

	for (i in 1:opt$n_top) {
		response <- as.factor(c(rep(1,length(test_set_g1)), rep(0, length(test_set_g2))))
		auc_top_LR[[i]] <- auc(response = response, predictor = indiv_predictions_EBADIMEX_cv_LR[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_LR[[i]] <- auc(predictor=apply(indiv_predictions_EBADIMEX_cv_LR[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_LR[[i]] <- auc_top_LR[[i]]
		
		auc_top_BLR[[i]] <- auc(response = response, predictor = indiv_predictions_EBADIMEX_cv_BLR[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_BLR[[i]] <- auc(predictor=apply(indiv_predictions_EBADIMEX_cv_BLR[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_BLR[[i]] <- auc_top_BLR[[i]]
		
		auc_top_RF[[i]] <- auc(response = response, predictor = indiv_predictions_EBADIMEX_cv_RF[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_RF[[i]] <- auc(predictor=apply(indiv_predictions_EBADIMEX_cv_RF[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_RF[[i]] <- auc_top_RF[[i]]
		
		auc_top_SVMLin[[i]] <- auc(response = response, predictor = indiv_predictions_EBADIMEX_cv_SVMLin[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_SVMLin[[i]] <- auc(predictor=apply(indiv_predictions_EBADIMEX_cv_SVMLin[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_SVMLin[[i]] <- auc_top_SVMLin[[i]]
		
		auc_top_SVMPol[[i]] <- auc(response = response, predictor = indiv_predictions_EBADIMEX_cv_SVMPol[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_SVMPol[[i]] <- auc(predictor=apply(indiv_predictions_EBADIMEX_cv_SVMPol[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_SVMPol[[i]] <- auc_top_SVMPol[[i]]
	}
	AUCs[[cv_run]] <- list(auc_top_combined_SVMLin, auc_top_combined_SVMPol, auc_top_combined_RF, auc_top_combined_LR, auc_top_combined_BLR)
	names(AUCs[[cv_run]]) <- c("auc_top_combined_SVMLin", "auc_top_combined_SVMPol", "auc_top_combined_RF", "auc_top_combined_LR", "auc_top_combined_BLR")
}
save(AUCs, file=paste(opt$name, ".RData", sep=""))
