expr_readCounts <- expr * genedata$lib_size
UQs <- apply(expr_readCounts, 2, FUN= function(x) quantile(x=x[x!=0], probs=0.75))
