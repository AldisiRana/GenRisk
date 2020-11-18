#! /usr/bin/Rscript
# #if (!require("pacman")) install.packages("pacman", repos = "https://CRAN.R-project.org/")
#pacman::p_load(data.table, 'betareg', scales, parallel, optparse, tidyverse)

if (!require("data.table")) install.packages("data.table", repos = "https://CRAN.R-project.org/")
library(data.table)
if (!require(scales)) install.packages("scales", repos = "https://CRAN.R-project.org/")
library(scales)
if (!require(parallel)) install.packages("parallel", repos = "https://CRAN.R-project.org/")
library(parallel)
if (!require(optparse)) install.packages("optparse", repos = "https://CRAN.R-project.org/")
library(optparse)
if (!require("betareg")) install.packages("betareg", repos = "https://CRAN.R-project.org/")
library('betareg')
if (!require(purrr)) install.packages("purrr", repos = "https://CRAN.R-project.org/")
library(purrr)
if (!require(utils)) install.packages("utils", repos = "https://CRAN.R-project.org/")
library(utils)


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
  make_option(c("--covariates"), type='character', help="all covariates for calculation, seperated by comma", default="PC1,PC2,age"),
  make_option(c("--nprocesses"), type='integer', default=3)
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (any(is.na(opt))) stop("Please make sure to include all required args.")

message("Reading files...")

mydata=read.table(opt$scoresfile, header=TRUE)
#mydata = mydata[, .SD, .SDcols = unique(names(mydata))]
#cols <- colnames(mydata)[grep('^[0-9]', colnames(mydata))]
#new_cols <- paste0("g", cols)
#colnames(mydata)[cols] <- new_cols

mydata[is.na(mydata)] = 0
mydata = Filter(var, mydata)
mydata = mydata[colMeans(mydata == 0) <= 0.9]

pheno=fread(opt$phenofile)
pc=fread(opt$pcfile)

covariates = c(opt$casescol, strsplit(opt$covariates, ",")[[1]])

epsilon=0.001



normalize <- function(gene)
{
  rescaled=rescale(gene,c(0+epsilon,1-epsilon))
  return(rescaled)
}



get_beta_pvals <- function(x, data) {
  
  #cols = c(x, covariates)
  #form <- as.formula(paste(x, paste(covariates, collapse = "+"), sep = "~"))
  #x = gsub(" ", "", x, fixed = TRUE)
  form <- paste(x, paste(" ."), sep = " ~")
  betaMod <- betareg(form, data=data)
  coefficient=betaMod$coefficients$mean[2]
  pval=coef(summary(betaMod))$mean[2,4]
  stderr=coef(summary(betaMod))$mean[2,2]
  return(c(x,coefficient,pval,stderr))
}

message("rescaling scores ...")
output=apply(mydata,2,normalize)
output=as.data.table(output)
output[, 1] = mydata[,1]
rm(mydata)

message("merging dataframes")
merged=merge(output,pheno,by=opt$samplescol)
rm(pheno)
completed=merge(merged,pc,by=opt$samplescol)
rm(merged)
rm(pc)
completed<-completed[complete.cases(completed),]
message("Calculating pvalues ...")
varlist <- names(completed)[2:ncol(output)]
rm(output)
#cl <- makeCluster(opt$nprocesses)

write("gene\tcoeff\tpval\stderr", file=opt$outputfile)

apply_betareg <- function(x){
  cols = c(x, covariates)
  data=completed[, ..cols]
  model=get_beta_pvals(x, data)
  write(paste(model, collapse = "\t"), file=opt$outputfile, append=TRUE)
  return(model)
  
}

#clusterEvalQ(cl, {
#  library(scales)
#  library(data.table)
#  library(betareg)
#})
#i <- 1
#pb = txtProgressBar(min = 0, max = length(varlist), initial = 0) 
#models <- c()
#for (x in varlist){
#  message(setTxtProgressBar(pb,i))
#  cols = c(x, covariates)
#  data=completed[, ..cols]
#  model = possibly(get_beta_pvals(x, data), otherwise = NA, quiet = TRUE)
#  models <- c(models,model)
#  i <- i+1
#}

#clusterExport(cl, c("completed", "covariates"))
#rm(completed)
models = mclapply(varlist, possibly(apply_betareg,NA_real_), mc.cores = opt$nprocesses, mc.preschedule = FALSE)
#models = parLapply(cl,varlist,possibly(get_beta_pvals,NA_real_))
#models = lapply(varlist,possibly(get_beta_pvals,NA_real_))

#message("saving to output file ...")
#pvals_df = data.frame(Reduce(rbind, models))
#head(pvals_df)
#colnames(pvals_df) <- c("gene", "coeff","pval","stderr")
#write.table(pvals_df, file=opt$outputfile, quote=FALSE, sep='\t', row.names = FALSE)
message("Done")
