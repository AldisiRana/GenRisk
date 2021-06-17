#!/usr/bin/Rscript
# #if (!require("pacman")) install.packages("pacman", repos = "https://CRAN.R-project.org/")
#pacman::p_load(data.table, 'betareg', scales, parallel, optparse, tidyverse)

if (!require("data.table")) install.packages("data.table", repos = "https://CRAN.R-project.org/")
library(data.table)
if (!require(scales)) install.packages("scales", repos = "https://CRAN.R-project.org/")
library(scales)
if (!require(optparse)) install.packages("optparse", repos = "https://CRAN.R-project.org/")
library(optparse)
if (!require("betareg")) install.packages("betareg", repos = "https://CRAN.R-project.org/")
library('betareg')
if (!require(purrr)) install.packages("purrr", repos = "https://CRAN.R-project.org/")
library(purrr)
if (!require(utils)) install.packages("utils", repos = "https://CRAN.R-project.org/")
library(utils)
if (!require(parallel)) install.packages("parallel", repos = "https://CRAN.R-project.org/")
library(parallel)
if (!require(lmtest)) install.packages("lmtest", repos = "https://CRAN.R-project.org/")


option_list = list(
  make_option(c("-s", "--scoresfile"), type="character",
              help="scores file path", metavar="character", default = NA),
  make_option(c("-o", "--outputfile"), type="character", default="output.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--phenofile"), type="character",
              help="phenotypes file path", default = NULL),
  make_option(c("--samplescol"), type='character', help="name of samples column", default="IID"),
  make_option(c("--casescol"), type='character', help="name of cases column", default="cases"),
  make_option(c("--processes"), type='integer', help="number of processes", default=1),
  make_option(c("--covariates"), type='character', help="all covariates for calculation, seperated by comma", default="PC1,PC2,age")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (any(is.na(opt))) stop("Please make sure to include all required args.")

message("Reading files...")

mydata=read.table(opt$scoresfile, header=TRUE)

mydata[is.na(mydata)] = 0

pheno=fread(opt$phenofile)

covariates = c(opt$casescol, strsplit(opt$covariates, ",")[[1]])

epsilon=0.001

normalize <- function(gene)
{
    gene = tryCatch(as.numeric(gene), error=function(err) gene, warning=function(w) gene)
    rescaled=tryCatch(rescale(gene,c(0+epsilon,1-epsilon)), error=function(err) gene)
    return(rescaled)
}


message("rescaling scores ...")
output=apply(mydata,2,normalize)
output=as.data.table(output)
output[[opt$samplescol]] = mydata[[opt$samplescol]]
rm(mydata)

message("merging dataframes")
completed=merge(output,pheno,by=opt$samplescol)
rm(pheno)


genes_list <- names(completed)[2:ncol(output)]
cols <- c(opt$samplescol, genes_list, covariates)
completed <- completed[, ..cols]
completed<-completed[complete.cases(completed),]
rm(output)

output <- opt$outputfile
write("genes\tcoeff\tp_value\tstderr", file=output)


apply_betareg <- function(x){
  cols = c(x, covariates)
  data=completed[, ..cols]
  form <- paste(x, paste(" ."), sep = " ~")
  betaMod <- betareg(form, data=data)
  coefficient=tryCatch(betaMod$coefficients$mean[2], error=function(err) NA)
  pval=tryCatch(coeftest(betaMod)[2,4], error=function(err) NA)
  stderr=tryCatch(coeftest(betaMod)[2,2], error=function(err) NA)
  results = c(x,coefficient,pval,stderr)
  write(paste(results, collapse = "\t"), file=output, append=TRUE)
  return(results)
}

message("Calculating pvalues ...")
cl <- makeCluster(opt$processes)

clusterEvalQ(cl, {
    library(data.table)
    library('betareg')
    library(purrr)
    library(utils)
    library(lmtest)
})

clusterExport(cl, c("covariates", "completed", "output"))

rm(completed)

models <- parLapply(cl, genes_list, possibly(apply_betareg,NA_real_))
#models = mclapply(genes_list, possibly(apply_betareg,NA_real_), mc.cores=opt$processes)

stopCluster(cl)

message("Done")
