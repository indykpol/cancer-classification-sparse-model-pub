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
opt <- list(data = "data/RData/data_BRCA_progressing_lazyReady.R_CACHE", results="data/latestResults/top_100_each_seed.RData",  name = "MLs_EBADIMEX_50randomSamples_BRCA_progressing_mean_UQ_brob", cv_runs = 50, n_top = 100, parallel=16)

library(caret)
library(pROC)
library(dplyr)
library(Brobdingnag)
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
			genedata$mean_pr <- apply(select(genedata, starts_with("pr_")), 1, mean)
			genedata$mean_gb <- apply(select(genedata, starts_with("gb_")), 1, mean)
			genedata_training <- genedata %>% 
				dplyr::select(matches("mean_pr"), matches("mean_gb"), matches("expr_UQ")) %>%
				slice(c(train_set_g1, train_set_g2))
			genedata_testing <- genedata %>% 
				dplyr::select(matches("mean_pr"), matches("mean_gb"), matches("expr_UQ")) %>%
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
	
	# As there was a problem multiplying posterior probs for LR models, I implemented a combination of posterior probs using Brobdingnag package. For consistency and to be on the safe side, I use this implementation throughout all methods now
	replaceInfs <- function(x) { # in case the product or a prob. is actual 0, the exponent becomes -Inf, which is not compatible with AUC calculation; replace it with the lowest number
		x[is.infinite(x)] <- min(x[!is.infinite(x)]) / 2
		x
	}
	for (i in 1:opt$n_top) {
		response <- as.factor(c(rep(1,length(test_set_g1)), rep(0, length(test_set_g2))))
		auc_top_LR[[i]] <- auc(response = response, predictor = replaceInfs(as.brob(indiv_predictions_EBADIMEX_cv_LR[i, c(test_set_g1, test_set_g2)])@x))
		if (i > 1) auc_top_combined_LR[[i]] <- auc(predictor=replaceInfs(apply(apply(indiv_predictions_EBADIMEX_cv_LR[1:i, c(test_set_g1, test_set_g2)], 1, FUN=function(x) as.brob(x)@x), 1, sum)), response=response) else auc_top_combined_LR[[i]] <- auc_top_LR[[i]]
		
		auc_top_BLR[[i]] <- auc(response = response, predictor = replaceInfs(as.brob(indiv_predictions_EBADIMEX_cv_BLR[i, c(test_set_g1, test_set_g2)])@x))
		if (i > 1) auc_top_combined_BLR[[i]] <- auc(predictor=replaceInfs(apply(apply(indiv_predictions_EBADIMEX_cv_BLR[1:i, c(test_set_g1, test_set_g2)], 1, FUN=function(x) as.brob(x)@x), 1, sum)), response=response) else auc_top_combined_BLR[[i]] <- auc_top_BLR[[i]]
		
		auc_top_RF[[i]] <- auc(response = response, predictor = replaceInfs(as.brob(indiv_predictions_EBADIMEX_cv_RF[i, c(test_set_g1, test_set_g2)])@x))
		if (i > 1) auc_top_combined_RF[[i]] <- auc(predictor=replaceInfs(apply(apply(indiv_predictions_EBADIMEX_cv_RF[1:i, c(test_set_g1, test_set_g2)], 1, FUN=function(x) as.brob(x)@x), 1, sum)), response=response) else auc_top_combined_RF[[i]] <- auc_top_RF[[i]]
		
		auc_top_SVMLin[[i]] <- auc(response = response, predictor = replaceInfs(as.brob(indiv_predictions_EBADIMEX_cv_SVMLin[i, c(test_set_g1, test_set_g2)])@x))
		if (i > 1) auc_top_combined_SVMLin[[i]] <- auc(predictor=replaceInfs(apply(apply(indiv_predictions_EBADIMEX_cv_SVMLin[1:i, c(test_set_g1, test_set_g2)], 1, FUN=function(x) as.brob(x)@x), 1, sum)), response=response) else auc_top_combined_SVMLin[[i]] <- auc_top_SVMLin[[i]]
		
		auc_top_SVMPol[[i]] <- auc(response = response, predictor = replaceInfs(as.brob(indiv_predictions_EBADIMEX_cv_SVMPol[i, c(test_set_g1, test_set_g2)])@x))
		if (i > 1) auc_top_combined_SVMPol[[i]] <- auc(predictor=replaceInfs(apply(apply(indiv_predictions_EBADIMEX_cv_SVMPol[1:i, c(test_set_g1, test_set_g2)], 1, FUN=function(x) as.brob(x)@x), 1, sum)), response=response) else auc_top_combined_SVMPol[[i]] <- auc_top_SVMPol[[i]]
	}
	AUCs[[cv_run]] <- list(auc_top_combined_SVMLin, auc_top_combined_SVMPol, auc_top_combined_RF, auc_top_combined_LR, auc_top_combined_BLR)
	names(AUCs[[cv_run]]) <- c("auc_top_combined_SVMLin", "auc_top_combined_SVMPol", "auc_top_combined_RF", "auc_top_combined_LR", "auc_top_combined_BLR")
}
save(AUCs, file=paste(opt$name, ".RData", sep=""))

# Analyze AUC tables
auc_top_combined_LR <- data.frame(model="auc_top_combined_LR", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_LR"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_LR"]]), 1, sd))
auc_top_combined_BLR <- data.frame(model="auc_top_combined_BLR", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_BLR"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_BLR"]]), 1, sd))
auc_top_combined_RF <- data.frame(model="auc_top_combined_RF", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_RF"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_RF"]]), 1, sd))
auc_top_combined_SVMLin <- data.frame(model="auc_top_combined_SVMLin", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMLin"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMLin"]]), 1, sd))
auc_top_combined_SVMPol <- data.frame(model="auc_top_combined_SVMPol", rank=1:100, AUC_mean =apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMPol"]]), 1, mean) , AUC_SD = apply(sapply(1:50, FUN = function(x) AUCs[[x]][["auc_top_combined_SVMPol"]]), 1, sd))
# write.table(x=cbind(top=1:opt$n_top, auc_top_LR, auc_top_combined_LR, auc_top_BLR, auc_top_combined_BLR, auc_top_RF, auc_top_combined_RF, auc_top_SVMLin, auc_top_combined_SVMLin, auc_top_SVMPol, auc_top_combined_SVMPol), row.names=FALSE, sep="\t", quote=FALSE, file=paste(opt$name, ".tab", sep=""))

require(ggplot2)
library(reshape2)
library(dplyr)

plotter <- rbind(auc_top_combined_LR, auc_top_combined_BLR, auc_top_combined_RF, auc_top_combined_SVMLin, auc_top_combined_SVMPol)
plotter$ymin <- plotter$AUC_mean - plotter$AUC_SD
plotter$ymax <- plotter$AUC_mean + plotter$AUC_SD
plotter <- mutate(plotter, section=ifelse(rank>21, "Top 21-100 running ranks combined", "Top 1-20 running ranks combined"))

plot <- ggplot(plotter, aes(x=rank, y=AUC_mean, colour=model)) + geom_line(lwd=1.5) + geom_ribbon(aes(ymin = ymin, ymax = ymax, colour = model), alpha=0.2) + facet_grid(. ~section, scales="free_x") + theme_bw() + ylab("AUC") + xlab("Top ranks combined") + ggtitle("tmp") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white"))

# with EBADIMEX AUCs
plotter <- rbind(auc_top_combined_LR, auc_top_combined_BLR, auc_top_combined_RF, auc_top_combined_SVMLin, auc_top_combined_SVMPol, auc_top_combined_EBADIMEX_expr, auc_top_combined_EBADIMEX_meth, auc_top_combined_EBADIMEX_joint, auc_top_combined_PINCAGE, auc_top_combined_PINCAGE_LS)
plotter$ymin <- plotter$AUC_mean - plotter$AUC_SD
plotter$ymax <- plotter$AUC_mean + plotter$AUC_SD
plotter <- mutate(plotter, section=ifelse(rank>21, "Top 21-100 running ranks combined", "Top 1-20 running ranks combined"))


plot <- ggplot(plotter, aes(x=rank, y=AUC_mean, colour=model, lwd=model, alpha=model)) + geom_line() + facet_grid(. ~section, scales="free_x") + theme_bw() + ylab("AUC") + xlab("Top ranks combined") + ggtitle("BRCA progression analysis using top 100 EBADIMEX genes") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background =element_rect(fill="white")) + scale_size_manual(values=c(rep(0.75,5),rep(1.75,5))) + scale_alpha_manual(values=c(0.75,0.75,0.75,0.75,0.75,0.85,0.85, 0.85, 0.85, 0.85))
plot

