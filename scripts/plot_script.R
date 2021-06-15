#!/usr/bin/Rscript
# if (!require("pacman")) install.packages("pacman", repos = "https://CRAN.R-project.org/")
#pacman::p_load(data.table, qqman, ggplot2, grid, optparse, gridGraphics)

if (!require("qqman")) install.packages("qqman", repos = "https://CRAN.R-project.org/")
library(qqman)
if (!require("data.table")) install.packages("data.table", repos = "https://CRAN.R-project.org/")
library(data.table)
if (!require("ggplot2")) install.packages("ggplot2", repos = "https://CRAN.R-project.org/")
library(ggplot2)
if (!require("grid")) install.packages("grid", repos = "https://CRAN.R-project.org/")
library(grid)
if (!require("gridGraphics")) install.packages("gridGraphics", repos = "https://CRAN.R-project.org/")
library(gridGraphics)
if (!require("optparse")) install.packages("optparse", repos = "https://CRAN.R-project.org/")
library(optparse)

option_list = list(
  make_option(c("-p", "--pvals_file"), type="character",
              help="pvalues file path", metavar="character"),
  make_option(c("-i", "--info_file"), type="character",
              help="info file path", metavar="character"),
  make_option(c("--pvalcol"), default='p_value', type="character", help='name of p_values column'),
  make_option(c("--genescol_1"), type='character', help="name of genes column in pvals file", default="gene"),
  make_option(c("--genescol_2"), type='character', help="name of genes column in info file", default="Gene.refGene"),
  make_option(c("--manhattan_output"), type='character', help="path to output manhattan plot image, should end with .jpeg"),
  make_option(c("--qq_output"), type='character', help="path to output manhattan plot image, should end with .jpeg"),
  make_option(c("--chr_col"), type='character', default='Chr', help='name of column with Chr number'),
  make_option(c("--pos_col"), type='character', default='Start', help='name of column with start position')
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (any(is.na(opt))) stop("Please make sure to include all required args.")

message("Reading files...")

pvals_df=fread(opt$pvals_file)
info_df=fread(opt$info_file)
complete_df=merge(pvals_df,info_df,by.x=opt$genescol_1,by.y=opt$genescol_2)
complete_df = unique(complete_df, by=opt$genescol_1)
complete_df[complete_df == 'X'] = 23
complete_df[complete_df == 'Y'] = 24
complete_df$Chr <- as.integer(complete_df$Chr)

plot_max = -log10(1e-09)
p = opt$pvalcol
if (-log10(min(complete_df[, ..p], na.rm=T)) > -log10(1e-09)){
  plot_max = -log10(min(complete_df[, ..p], na.rm=T))
}
  
jpeg(opt$manhattan_output, res=300, width = 12, height = 6, units = 'in')
manhattan(complete_df,chr=opt$chr_col, bp=opt$pos_col, snp=opt$genescol_1, p=opt$pvalcol,
          ylim = c(0, plot_max), chrlabs = NULL,
          suggestiveline = -log10(1.00e-05), genomewideline = -log10(5.00e-08), logp = TRUE, main=opt$pvals_file, highlight=T, annotatePval=1, annotateTop=T)
dev.off()
lambda=median(qchisq(na.omit(complete_df[[opt$pvalcol]]), df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
jpeg(opt$qq_output, res=300, width = 6, height = 6, units = 'in')
qq(na.omit(complete_df[[opt$pvalcol]]),main=toString(lambda))
dev.off()


