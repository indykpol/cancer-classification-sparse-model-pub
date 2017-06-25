##################################################
# list of utility functions
##################################################

read_element <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", element_name = "counts") {
	require(SOAR)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- get(element_name)
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

read_methylation_matrix_gene <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", genename = "A1BG") {
	require(SOAR)
	require(dplyr)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- as.data.frame(get(genename)) %>% select(starts_with("pr_"), starts_with("gb_")) %>% as.matrix()
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

read_grouping <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", group1 = "ANs", group2 = "Ts") {
	require(SOAR)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- list(get(group1), get(group2))
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

read_genedata <- function(cache_name = "data/RData/data_BRCA_all_lazyReady.R_CACHE", genename = "A1BG") {
	require(SOAR)
	oldLC <- Sys.getenv("R_LOCAL_CACHE", unset = ".R_Cache")
	Sys.setenv(R_LOCAL_CACHE=cache_name)
	Objects()
	tmp <- get(genename)
	Sys.setenv(R_LOCAL_CACHE=oldLC)
	return(tmp)
}

fishersMethod <- function(x, df=2*length(x[!is.na(x)]), log.p = FALSE, logs = FALSE) {
	x <- x[!is.na(x)]
	if (!logs) pchisq(-2*sum(log(x)), df, lower.tail = FALSE, log.p = log.p) else pchisq(-2*sum(x), df, lower.tail = FALSE, log.p = log.p)
}

# From DescTools package
Winsorize <- function(x, minval = NULL, maxval = NULL,
                      probs=c(0.05, 0.95), na.rm = TRUE) {

  # following an idea from Gabor Grothendieck
  # http://r.789695.n4.nabble.com/how-to-winsorize-data-td930227.html

  # don't eliminate NAs in x, moreover leave them untouched,
  # just calc quantile without them...

  if(is.null(minval) || is.null(maxval)){
    xq <- quantile(x=x, probs=probs, na.rm = na.rm)
    if(is.null(minval)) minval <- xq[1]
    if(is.null(maxval)) maxval <- xq[2]
  }

  x[x<minval] <- minval
  x[x>maxval] <- maxval

  return(x)
}

methylationWinsor <- function(x) {return(Winsorize(x, probs = c(0.01,0.99, na.rm=TRUE)))}