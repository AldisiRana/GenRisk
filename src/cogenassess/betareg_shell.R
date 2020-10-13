
if (!require("pacman")) install.packages("pacman", repos = "https://CRAN.R-project.org/")
pacman::p_load(data.table, betareg, scales, parallel, optparse, tidyverse)

library(data.table)
library(betareg)
library(scales)
library(parallel)
library(optparse)
library(tidyverse)

option_list = list(
  make_option(c("-s", "--scoresfile"), type="character",
              help="scores file path", metavar="character", default = NA),
  make_option(c("-o", "--outputfile"), type="character", default="output.txt", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("--phenofile"), type="character",
              help="phenotypes file path", default = NULL),
  make_option(c("--pcfile"), type='character', help="covariates file path", default = NULL),
  make_option(c("--samplescol"), type='character', help="name of samples column", default="IID"),
  make_option(c("--casescol"), type='character', help="name of cases column", default="cases"),
  make_option(c("--covariates"), type='character', help="all covariates for calculation, seperated by comma", default="PC1,PC2,age.y"),
  make_option(c("--nprocesses"), type='integer', default=5)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (any(is.na(opt))) stop("Please make sure to include all required args.")

message("Reading files...")

mydata=fread(opt$scoresfile)
mydata = mydata[, .SD, .SDcols = unique(names(mydata))]
pheno=fread(opt$phenofile)
pc=fread(opt$pcfile)

covariates = c(opt$casescol, as.list(strsplit(opt$covariates, ",")[[1]]))

epsilon=0.001

mydata[is.na(mydata)] = 0

normalize <- function(gene)
{
  rescaled=rescale(gene,c(0+epsilon,1-epsilon))
  return(rescaled)
}



get_beta_pvals <- function(x) {

  form <- as.formula(paste(x, paste(covariates, collapse = "+"), sep = "~"))
  betaMod <- betareg(form, data = completed)
  coefficient=betaMod$coefficients$mean[2]
  pval=coef(summary(betaMod))$mean[2,4]
  stderr=coef(summary(betaMod))$mean[2,2]
  return(c(x,coefficient,pval,stderr))
}

message("rescaling scores ...")
output=apply(mydata,2,normalize)
output=as.data.frame(output)
output[, 1] = mydata[,1]

message("merging dataframes")
merge=merge(output,pheno,by=opt$samplescol)
completed=merge(merge,pc,by=opt$samplescol)
completed<-completed[complete.cases(completed),]

message("Calculating pvalues ...")
varlist <- names(completed)[2:ncol(output)]
cl <- makeCluster(opt$nprocesses)


clusterEvalQ(cl, {
  library(scales)
  library(data.table)
  library(betareg)
})


clusterExport(cl, c("completed", "covariates"))

models = parLapply(cl,varlist,possibly(get_beta_pvals,NA_real_))

message("saving to output file ...")
pvals_df = data.frame(Reduce(rbind, models))
colnames(pvals_df) <- c("gene", "coeff","pval","stderr")
write.table(pvals_df, file=opt$outputfile, quote=FALSE, sep='\t', row.names = FALSE)
message("Done")
