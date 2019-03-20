#!/usr/bin/env Rscript
library(optparse)
# Command Line Arguments
option_list <- list(
  make_option(c("-d","--data"), type='character', help = 'cached data to use', default = "data/RData/data_LUSC_all_lazyReady.R_CACHE"),
  make_option(c("-r","--results"), type='character', help = 'model results to use', default = "data/latestResults/top_100_each_seed.RData"),
  make_option(c("-f","--folds"), type='integer', help = 'number of folds to split the dataset', default = 4),
  make_option(c("-n","--name"), type='character', help = 'name of the analysis performed', default = "MLs_EBADIMEX_50randomSamples_BRCA_progressing"),
  make_option("--n_top", type='integer', help = 'number of top ranks/genes to train for and combine', default = 100),
  make_option(c("-p","--parallel"), type='integer', help = 'if desired, provide a number of cores to parallelize the execution of the analysis', default = 1)
)
#opt <- parse_args(OptionParser(option_list = option_list))
opt <- list(data = "data/RData/data_BRCA_progressing_lazyReady.R_CACHE", results="data/latestResults/top_100_each_seed.RData",  name = "MLs_EBADIMEX_50randomSamples_BRCA_progressing", folds = 50, n_top = 100, parallel=16)

library(caret)
library(pROC)
library(dplyr)
library(parallel)
source("R/utilities_new.R")

load(opt$results)
samplenames <- unlist(read_grouping(cache_name = opt$data, group1="progressed", group2="nonProgressed"))
indiv_predictions_Irksome_cv <- matrix(ncol=length(samplenames), nrow=opt$n_top)
colnames(indiv_predictions_Irksome_cv) <- samplenames
indiv_predictions_Irksome_cv_SVMPol <- indiv_predictions_Irksome_cv_SVMLin <- indiv_predictions_Irksome_cv_LR <- indiv_predictions_Irksome_cv_BLR <- indiv_predictions_Irksome_cv_RF <- indiv_predictions_Irksome_cv

AUCs <- NULL
for (fold in 1:50) {
	
	set.seed(fold)
	train_set_g1 <- sample(1:38, 10)
	train_set_g2 <- 38+sample(1:43, 15)
	test_set_g1 <- (1:38)[-train_set_g1]
	test_set_g2 <- 38+(1:43)[-(train_set_g2-38)]
	
	train_labels <- as.factor(c(rep("g1", length(train_set_g1)), rep("g2", length(train_set_g2))))
	
	current <- top_100_each_seed[[fold]]
	out <- mclapply(1:length(current), mc.cores = opt$parallel, FUN = function(gene) {
		genedata <- as.data.frame(read_genedata(cache = opt$data, genename = current[gene]))
		genedata$promoter_region <-  genedata %>% 
			select(starts_with("pr_")) %>% 
			rowMeans()
		genedata$geneBody_region <-  genedata %>% 
			select(starts_with("gb_")) %>% 
			rowMeans()
		
		genedata_training <- genedata %>% 
			dplyr::select(matches("expr_plusOne"), matches("promoter_region"), matches("geneBody_region")) %>%
			slice(c(train_set_g1, train_set_g2))
		genedata_testing <- genedata %>% 
			dplyr::select(matches("expr_plusOne"), matches("promoter_region"), matches("geneBody_region")) %>%
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


	auc_top_SVMLin <- auc_top_SVMPol <- auc_top_RF <- auc_top_LR <- auc_top_BLR <- NULL
	auc_top_combined_SVMLin <- auc_top_combined_SVMPol <- auc_top_combined_RF <- auc_top_combined_LR <- auc_top_combined_BLR <- NULL

	for (i in 1:opt$n_top) {
		response <- as.factor(c(rep(1,length(test_set_g1)), rep(0, length(test_set_g2))))
		auc_top_LR[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_LR[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_LR[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_LR[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_LR[[i]] <- auc_top_LR[[i]]
		
		auc_top_BLR[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_BLR[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_BLR[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_BLR[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_BLR[[i]] <- auc_top_BLR[[i]]
		
		auc_top_RF[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_RF[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_RF[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_RF[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_RF[[i]] <- auc_top_RF[[i]]
		
		auc_top_SVMLin[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_SVMLin[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_SVMLin[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_SVMLin[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_SVMLin[[i]] <- auc_top_SVMLin[[i]]
		
		auc_top_SVMPol[[i]] <- auc(response = response, predictor = indiv_predictions_Irksome_cv_SVMPol[i, c(test_set_g1, test_set_g2)])
		if (i > 1) auc_top_combined_SVMPol[[i]] <- auc(predictor=apply(indiv_predictions_Irksome_cv_SVMPol[1:i, c(test_set_g1, test_set_g2)],2,prod), response=response) else auc_top_combined_SVMPol[[i]] <- auc_top_SVMPol[[i]]
	}
	AUCs[[fold]] <- list(auc_top_combined_SVMLin, auc_top_combined_SVMPol, auc_top_combined_RF, auc_top_combined_LR, auc_top_combined_BLR)
	names(AUCs[[fold]]) <- c("auc_top_combined_SVMLin", "auc_top_combined_SVMPol", "auc_top_combined_RF", "auc_top_combined_LR", "auc_top_combined_BLR")
}
save(AUCs, file=paste(opt$name, ".RData", sep=""))

# Analyze AUC tables
auc_top_combined_LR <- data.frame(model="auc_top_combined_LR", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_LR"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_LR"]]), 1, sd))
auc_top_combined_BLR <- data.frame(model="auc_top_combined_BLR", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_BLR"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_BLR"]]), 1, sd))
auc_top_combined_RF <- data.frame(model="auc_top_combined_RF", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_RF"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_RF"]]), 1, sd))
auc_top_combined_SVMLin <- data.frame(model="auc_top_combined_SVMLin", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMLin"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMLin"]]), 1, sd))
auc_top_combined_SVMPol <- data.frame(model="auc_top_combined_SVMPol", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMPol"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMPol"]]), 1, sd))
auc_top_combined_EBADIMEX_expr <- data.frame(model="auc_top_combined_EBADIMEX_expr", rank=1:100, AUC_mean = unlist(lapply(1:100, FUN=function(x) mean(subset(df, n==x)[,1]))), AUC_SD = unlist(lapply(1:50, FUN=function(x) sd(subset(df, n==x)[,1]))))
auc_top_combined_EBADIMEX_meth <- data.frame(model="auc_top_combined_EBADIMEX_meth", rank=1:100, AUC_mean = unlist(lapply(1:100, FUN=function(x) mean(subset(df, n==x)[,2]))), AUC_SD = unlist(lapply(1:50, FUN=function(x) sd(subset(df, n==x)[,2]))))
auc_top_combined_EBADIMEX_joint <- data.frame(model="auc_top_combined_EBADIMEX_joint", rank=1:100, AUC_mean = unlist(lapply(1:100, FUN=function(x) mean(subset(df, n==x)[,3]))), AUC_SD = unlist(lapply(1:50, FUN=function(x) sd(subset(df, n==x)[,3]))))

require(ggplot2)
library(reshape2)
library(dplyr)

plotter <- rbind(auc_top_combined_LR, auc_top_combined_BLR, auc_top_combined_RF, auc_top_combined_SVMLin, auc_top_combined_SVMPol, auc_top_combined_EBADIMEX_expr, auc_top_combined_EBADIMEX_meth, auc_top_combined_EBADIMEX_joint)
plotter$ymin <- plotter$AUC_mean - plotter$AUC_SD
plotter$ymax <- plotter$AUC_mean + plotter$AUC_SD
plotter <- mutate(plotter, section=ifelse(rank>21, "Top 21-100 running ranks combined", "Top 1-20 running ranks combined"))

plot <- ggplot(plotter, aes(x=rank, y=AUC_mean, colour=model)) + geom_line(lwd=1.5) + geom_ribbon(aes(ymin = ymin, ymax = ymax, colour = model), alpha=0.2) + facet_grid(. ~section, scales="free_x") + theme_bw() + ylab("AUC") + xlab("Top ranks combined") + ggtitle("tmp") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"))
plot <- ggplot(plotter, aes(x=rank, y=AUC_mean, colour=model, lwd=model, alpha=model)) + geom_line() + facet_grid(. ~section, scales="free_x") + theme_bw() + ylab("AUC") + xlab("Top ranks combined") + ggtitle("BRCA progression analysis using top 100 EBADIMEX genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white")) + scale_size_manual(values=c(1,1,1,1,1,1.75,1.75,1.75)) + scale_alpha_manual(values=c(0.75,0.75,0.75,0.75,0.75,0.85,0.85, 0.85))
# plot <- ggplot(plotter, aes(x=rank, y=AUC_mean, colour=model, lwd=model, alpha=model)) + geom_line() + facet_grid(. ~section, scales="free_x") + theme_bw() + ylab("AUC") + xlab("Top ranks combined") + ggtitle("tmp") + scale_y_continuous(limits = c(0.45, 0.65)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white")) + scale_size_manual(values=c(2,2,1,1,1,1,1)) + scale_alpha_manual(values=c(0.75,0.75,0.5,0.5,0.5,0.5,0.5))
plot

# write.table(x=cbind(top=1:opt$n_top, auc_top_LR, auc_top_combined_LR, auc_top_BLR, auc_top_combined_BLR, auc_top_RF, auc_top_combined_RF, auc_top_SVMLin, auc_top_combined_SVMLin, auc_top_SVMPol, auc_top_combined_SVMPol), row.names=FALSE, sep="\t", quote=FALSE, file=paste(opt$name, ".tab", sep=""))