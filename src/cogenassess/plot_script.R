
if (!require("pacman")) install.packages("pacman", repos = "https://CRAN.R-project.org/")
pacman::p_load(data.table, qqman, ggplot2, grid, optparse, gridGraphics)

library(qqman)
library(data.table)
library(ggplot2)
library(grid)
library(gridGraphics)
library(optparse)

option_list = list(
  make_option(c("-p", "--pvals_file"), type="character",
              help="pvalues file path", metavar="character"),
  make_option(c("-i", "--info_file"), type="character",
              help="info file path", metavar="character"),
  make_option(c("--genescol_1"), type='character', help="name of genes column in pvals file", default="gene"),
  make_option(c("--genescol_2"), type='character', help="name of genes column in info file", default="Gene.refGene"),
  make_option(c("-o", "--output_file"), type='character')
); 


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (any(is.na(opt))) stop("Please make sure to include all required args.")

pdf(opt$output_file)

message("Reading files...")

pvals_df=fread(opt$pvals_file)
info_df=fread(opt$info_file)
complete_df=merge(pvals_df,info_df,by.x=opt$genescol_1,by.y=opt$genescol_2)
complete_df = unique(complete_df, by='gene')
complete_df[complete_df == 'X'] = 23
complete_df[complete_df == 'Y'] = 24
complete_df$Chr <- as.integer(complete_df$Chr)

manhattan(complete_df, chr="Chr", bp="Start", snp="gene", p="pval", 
          ylim = c(0, -log10(1e-06)), chrlabs = NULL,
          suggestiveline = -log10(1e-03), genomewideline = -log10(1e-05),
          highlight = NULL, logp = TRUE, annotatePval = 0.005,
          annotateTop = TRUE, main=opt$pvals_file)

qq(complete_df$pval, main = opt$pvals_file)

