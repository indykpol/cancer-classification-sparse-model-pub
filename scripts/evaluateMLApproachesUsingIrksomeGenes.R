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

load(opt$results)
opt <- opts
rownames(ranks) <- rownames(expr_logUQ)
indiv_predictions_Irksome_cv <- matrix(ncol=length(samplenames), nrow=opt$n_top)
colnames(indiv_predictions_Irksome_cv) <- samplenames
indiv_predictions_Irksome_cv_SVMPol <- indiv_predictions_Irksome_cv_SVMLin <- indiv_predictions_Irksome_cv_LR <- indiv_predictions_Irksome_cv_BLR <- indiv_predictions_Irksome_cv_RF <- indiv_predictions_Irksome_cv

for (fold in 1:opt$folds) {
	train_set_g1 <- setdiff(g1, g1_cv_sets[[fold]])
	train_set_g2 <- setdiff(g2, g2_cv_sets[[fold]])
	test_set_g1 <- g1_cv_sets[[fold]]
	test_set_g2 <- g2_cv_sets[[fold]]
	train_labels <- as.factor(c(rep("g1", length(train_set_g1)), rep("g2", length(train_set_g2))))
	
	current <- names(ranks[(which(ranks[,fold] %in% 1:opt$n_top)), fold])
	out <- mclapply(1:length(current), mc.cores = opt$parallel, FUN = function(gene) {
		genedata <- as.data.frame(read_genedata(cache = opt$data, genename = current[gene]))
		genedata_training <- genedata %>% 
			dplyr::select(starts_with("pr_"), starts_with("gb_"), matches("expr_plusOne")) %>%
			slice(c(train_set_g1, train_set_g2))
		genedata_testing <- genedata %>% 
			dplyr::select(starts_with("pr_"), starts_with("gb_"), matches("expr_plusOne")) %>%
			slice(c(test_set_g1, test_set_g2))
		
		preProcValues <- preProcess(genedata_training, method = c("center", "scale"))
		genedata_training_CS <- predict(preProcValues, genedata_training)
		genedata_testing_CS <- predict(preProcValues, genedata_testing)
		set.seed(1)
		# Logistic Regression
		model <- train(x=genedata_training, y=train_labels, method="glm")
		indiv_predictions_Irksome_cv_LR[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing, type="prob")[,1]
		# Bayesian Logistic Regression
		model <- train(x=genedata_training, y=train_labels, method="bayesglm")
		indiv_predictions_Irksome_cv_BLR[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing, type="prob")[,1]
		# Random Forest
		model <- train(x=genedata_training, y=train_labels, method="rf")
		indiv_predictions_Irksome_cv_RF[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing, type="prob")[,1]
		# SVM with linear kernel
		model <- train(x=genedata_training_CS, y=train_labels, method="svmLinear2", trControl=trainControl(classProbs = TRUE))
		indiv_predictions_Irksome_cv_SVMLin[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing_CS, type="prob")[,1]
		# SVM with polynomial kernel
		model <- train(x=genedata_training_CS, y=train_labels, method="svmPoly", trControl=trainControl(classProbs = TRUE))
		indiv_predictions_Irksome_cv_SVMPol[gene, c(test_set_g1, test_set_g2)] <- predict(model, genedata_testing_CS, type="prob")[,1]
		list(indiv_predictions_Irksome_cv_LR[gene, c(test_set_g1, test_set_g2)], indiv_predictions_Irksome_cv_RF[gene, c(test_set_g1, test_set_g2)], indiv_predictions_Irksome_cv_BLR[gene, c(test_set_g1, test_set_g2)], indiv_predictions_Irksome_cv_SVMLin[gene, c(test_set_g1, test_set_g2)], indiv_predictions_Irksome_cv_SVMPol[gene, c(test_set_g1, test_set_g2)])
	})
	
	for (gene in 1:length(out)) { # transfer the results from the parallel execution
		indiv_predictions_Irksome_cv_LR[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[1]]
		indiv_predictions_Irksome_cv_RF[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[2]]
		indiv_predictions_Irksome_cv_BLR[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[3]]
		indiv_predictions_Irksome_cv_SVMLin[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[4]]
		indiv_predictions_Irksome_cv_SVMPol[gene, c(test_set_g1, test_set_g2)] <- out[[gene]][[5]]
	}
	rm(out)
}

auc_top_SVMLin <- auc_top_SVMPol <- auc_top_RF <- auc_top_LR <- auc_top_BLR <- NULL
auc_top_combined_SVMLin <- auc_top_combined_SVMPol <- auc_top_combined_RF <- auc_top_combined_LR <- auc_top_combined_BLR <- NULL

for (i in 1:opt$n_top) {
	response <- as.factor(c(rep(1,length(g1)), rep(0, length(g2))))
	auc_top_LR[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_LR[i,])
	if (i > 1) auc_top_combined_LR[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_LR[1:i,],2,prod), response=response)
	
	auc_top_BLR[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_BLR[i,])
	if (i > 1) auc_top_combined_BLR[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_BLR[1:i,],2,prod), response=response)
	
	auc_top_RF[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_RF[i,])
	if (i > 1) auc_top_combined_RF[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_RF[1:i,],2,prod), response=response)
	
	auc_top_SVMLin[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_SVMLin[i,])
	if (i > 1) auc_top_combined_SVMLin[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_SVMLin[1:i,],2,prod), response=response)
	
	auc_top_SVMPol[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_SVMPol[i,])
	if (i > 1) auc_top_combined_SVMPol[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_SVMPol[1:i,],2,prod), response=response)
}

write.table(x=cbind(top=1:opt$n_top, auc_top_LR, auc_top_combined_LR, auc_top_BLR, auc_top_combined_BLR, auc_top_RF, auc_top_combined_RF, auc_top_SVMLin, auc_top_combined_SVMLin, auc_top_SVMPol, auc_top_combined_SVMPol), row.names=FALSE, sep="\t", quote=FALSE, file=paste(opt$name, ".tab", sep=""))