#!/usr/bin/env Rscript
library(optparse)
option_list <- list(
  make_option(c("-d","--data"), type='character', help = 'cached data to work on', default = ' data/data_BRCA_progressing_lazyReady.R_CACHE'),
  make_option(c("-g","--gene_symbols"), type='character', help = 'a comma-separated list of genes to plot in given order'),
  make_option(c("-n","--name"), type='character', default = NULL, help = 'desired name of the output'),
  make_option(c("-s","--cases"), type='character', help = 'comma-separated list of names of sample groups present in the RData file to test against', default = 'progressed'),
  make_option(c("-c","--controls"), type='character', help = 'optional comma-separated list of names of the other sample groups present in the RData', default = "nonProgressed"),
  make_option("--alpha_cases", type='numeric', help = 'alpha value for cases', default = 1),
  make_option("--alpha_controls", type='numeric', help = 'alpha value for cases', default = 1),
  make_option("--expression_winsor", type='logical', help = 'should expression measurements be winsorized?', default = TRUE)
)
opt <- parse_args(OptionParser(option_list = option_list))
#opt <- list(data = "data/data_BRCA_progressing_lazyReady.R_CACHE", name="PINCAGE_BRCA_progressing", cases = "progressed", controls = "nonProgressed", gene_symbols="HEMK1,CCDC72")
library(ggplot2)
library(dplyr)
library(reshape2)
library(gridExtra)
source("/project/iCancerGenomics/faststorage/R/utilities_PINCAGE.R")
source("/project/iCancerGenomics/faststorage/R/utilities_new.R")

ptm <- proc.time()[3]

genes <- unlist(strsplit(opt$gene_symbols, ","))
grouping <- read_grouping(cache = opt$data, group1 = opt$cases, group2 = opt$controls)
G1 <- grouping[[1]]
G2 <- grouping[[2]]
colour_palette <- c("#56B4E9","#D55E00")
pdf(file=paste(opt$name, "_selectedGenes.pdf", sep=""),width=11.7,height=8.27)
for (i in 1:length(genes)) {
	cat(paste("Plotting", genes[i], "\n"))
	genedata <- read_genedata(cache_name = opt$data, genename = genes[i]) %>% as.data.frame() %>% cbind(cases=c(G1,G2) %in% G1)
	genedata$promoter_region <-  genedata %>% 
		select(starts_with("pr_")) %>% 
		rowMeans()
	genedata$geneBody_region <-  genedata %>% 
		select(starts_with("gb_")) %>% 
		rowMeans()
	if (opt$expression_winsor) genedata$expr <- Winsorize(genedata$expr)
	
	plotter_promoter_CpGs <- genedata %>% select(starts_with("pr_")) %>% cbind(cases = (c(G1,G2) %in% G1)) %>% cbind(region="Promoter region") %>% melt(id.vars=c("cases", "region"))
	plotter_geneBody_CpGs <- genedata %>% select(starts_with("gb_")) %>% cbind(cases = (c(G1,G2) %in% G1)) %>% cbind(region="Gene body region") %>% melt(id.vars=c("cases", "region"))
	n_CpG_pr <- genedata %>% 
		select(starts_with("pr_")) %>% 
		ncol()
	n_CpG_gb <- genedata %>% 
		select(starts_with("gb_")) %>% 
		ncol()
	plotter_region_CpGs <- rbind(plotter_promoter_CpGs, plotter_geneBody_CpGs)
	max_val <- max(plotter_region_CpGs$value)+0.25
	plot_CpGs <- ggplot(plotter_region_CpGs, aes(x=variable, y=value, colour=cases)) + geom_boxplot() + geom_segment(x=0 ,y=max_val,xend=n_CpG_pr+0.5,yend=max_val, colour="darkgreen") + geom_segment(x=n_CpG_pr+0.5, y=max_val, xend=n_CpG_pr+n_CpG_gb+1, yend=max_val, colour="darkred") + theme_bw()+ scale_colour_manual(values=colour_palette)+ theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8), legend.position="bottom") + xlab("CpG sites") + ylab("Methylation level (M-value)") + ggtitle(paste(genes[i],": regional methylation landscape"))
	plot_expression <- ggplot(genedata, aes(x=cases, y=expr*10^6, colour=cases)) + geom_boxplot(outlier.alpha=0) + theme_bw() + scale_colour_manual(values=colour_palette) + scale_alpha_discrete(range=c(opt$alpha_cases, opt$alpha_controls)) + geom_jitter(data=genedata, width=0.10, aes(alpha=cases,x=cases, y=expr*10^6)) + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8), legend.position="bottom") + ylab("Expression (reads per million)") + ggtitle(paste(genes[i],": gene expression"))
	
	correlation_g1 <- signif(cor(genedata[G1,"promoter_region"], genedata[G1,"geneBody_region"], method="spearman"),digits=2)
	correlation_g2 <- signif(cor(genedata[G2,"promoter_region"], genedata[G2,"geneBody_region"], method="spearman"),digits=2)
	plot_pr_gb <-  ggplot(genedata, aes(x=promoter_region,y=geneBody_region,colour=cases)) + stat_density2d(alpha=0.5)  + geom_point(alpha=0.25) + scale_colour_manual(values=colour_palette) + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8), legend.position="bottom") + xlab("Pr. meth. (M-value)") + ylab("GB. meth. (M-value)") + ggtitle(paste(genes[i],"; cor(cases)=", correlation_g1,", cor(controls)=", correlation_g2,", Spearman's rho",sep=""))
	
	correlation_g1 <- signif(cor(genedata[G1,"promoter_region"], genedata[G1,"expr"], method="spearman"),digits=2)
	correlation_g2 <- signif(cor(genedata[G2,"promoter_region"], genedata[G2,"expr"], method="spearman"),digits=2)
	plot_pr_expr <- ggplot(genedata, aes(x=promoter_region,y=expr*10^6,colour=cases)) + stat_density2d(alpha=0.5)  + geom_point(alpha=0.25) + scale_colour_manual(values=colour_palette) + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8), legend.position="bottom") + xlab("Pr. meth. (M-value)") + ylab("Expression (reads per million)") + ggtitle(paste(genes[i],"; cor(cases)=", correlation_g1,", cor(controls)=", correlation_g2,", Spearman's rho",sep=""))
	
	correlation_g1 <- signif(cor(genedata[G1,"geneBody_region"], genedata[G1,"expr"], method="spearman"),digits=2)
	correlation_g2 <- signif(cor(genedata[G2,"geneBody_region"], genedata[G2,"expr"], method="spearman"),digits=2)
	plot_gb_expr <- ggplot(genedata, aes(x=geneBody_region,y=expr*10^6,colour=cases)) + stat_density2d(alpha=0.5)  + geom_point(alpha=0.25) + scale_colour_manual(values=colour_palette) + theme_bw() + theme(panel.grid=element_line(linetype = 0), plot.title=element_text(size=8), legend.position="bottom") + xlab("GB. meth. (M-value)") + ylab("Expression (reads per million)") + ggtitle(paste(genes[i],"; cor(cases)=", correlation_g1,", cor(controls)=", correlation_g2,", Spearman's rho",sep=""))
	grid.arrange(arrangeGrob(plot_CpGs, plot_expression, ncol=2, widths=c(3,1)), arrangeGrob(plot_pr_expr, plot_gb_expr, plot_pr_gb, ncol=3), nrow=2, heights=c(3,2))
}
dev.off()
cat(paste("Done plotting in", round(proc.time()[3]-ptm), "seconds\n"))