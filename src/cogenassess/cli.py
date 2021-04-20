# -*- coding: utf-8 -*-
import os
import re

import click
import shutil

import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from .utils import get_gene_info, plink_process, combine_scores, create_model
from .pipeline import normalize_gene_len, find_pvalue, betareg_pvalues, r_visualize


@click.group()
def main():
    """Handle cogenassess functions."""


@main.command()
@click.option('-a', '--annotated-file', required=True, help='the annotated file')
@click.option('--bed', required=True, help="text file for genotype")
@click.option('--bim', required=True, help="file with variant information")
@click.option('--fam', required=True, help="text file for pedigree information")
@click.option('--plink', default='plink', help="the directory of plink, if not set in environment")
@click.option('-t', '--temp-dir', required=True, help="a temporary directory to save temporary files before merging.")
@click.option('-o', '--output-file', required=True, help="the final output scores matrix.")
@click.option('--beta-param', default=(1.0, 25.0), nargs=2, type=float,
              help="the parameters from beta weight function.")
@click.option('--weight-func', default='beta', type=click.Choice(['beta', 'log10']),
              help="the weighting function used in score calculation.")
@click.option('--variant-col', default='SNP', help="the column containing the variant IDs.")
@click.option('--gene-col', default='Gene.refGene', help="the column containing gene names.")
@click.option('--af-col', default='MAF', help="the column containing allele frequency.")
@click.option('--del-col', default='CADD_raw', help="the column containing the deleteriousness score.")
@click.option('--alt-col', default='Alt', help="the column containing the alternate base.")
@click.option('--maf-threshold', default=0.01, help="the threshold for minor allele frequency.")
@click.option('--remove-temp', is_flag=True,
              help="if flagged the temporary directory will be deleted after process completion.")
def score_genes(
    *,
    annotated_file,
    bed,
    bim,
    fam,
    plink,
    beta_param,
    temp_dir,
    output_file,
    weight_func,
    variant_col,
    gene_col,
    af_col,
    del_col,
    alt_col,
    maf_threshold,
    remove_temp,
):
    """

    :param annotated_file: a file containing
    :param bed: text file for genotype.
    :param bim: file with variant information
    :param fam: text file for pedigree information
    :param plink: the directory of plink, if not set in environment
    :param beta_param: the parameters from beta weight function.
    :param temp_dir: a temporary directory to save temporary files before merging.
    :param output_file: the final output scores matrix.
    :param weight_func: the weighting function used in score calculation.
    :param variant_col: the column containing the variant IDs.
    :param gene_col: the column containing gene names.
    :param af_col: the column containing allele frequency.
    :param del_col: the column containing deleteriousness score.
    :param alt_col: the column containing alternate base.
    :param maf_threshold: the threshold for minor allele frequency.
    :param remove_temp: if True temporary directory will be deleted after process completion.
    :return: the final dataframe information.
    """
    click.echo('getting information from vcf files')
    genes_folder = get_gene_info(
        annotated_file=annotated_file,
        output_dir=temp_dir,
        beta_param=beta_param,
        weight_func=weight_func,
        del_col=del_col,
        maf_threshold=maf_threshold,
        genes_col=gene_col,
        variant_col=variant_col,
        af_col=af_col,
        alt_col=alt_col,
    )
    click.echo('calculating gene scores ...')
    plink_process(genes_folder=genes_folder, plink=plink, bed=bed, bim=bim, fam=fam)
    click.echo('combining score files ...')
    df = combine_scores(input_path=temp_dir, output_path=output_file)
    if remove_temp:
        shutil.rmtree(temp_dir)
    click.echo('process is complete.')
    return df.info()


@main.command()
@click.option('--bed', required=True, help="text file for genotype")
@click.option('--bim', required=True, help="file with variant information")
@click.option('--fam', required=True, help="text file for pedigree information")
@click.option('--plink', default='plink')
@click.option('--genes-folder', required=True, help="a folder that contains two files for each gene, w and v files.")
def run_plink(*, genes_folder, plink, bed, bim, fam):
    """
    Get the genes' scores from a folder of genes info.
    :param genes_folder: a folder that contains two files for each gene,
    one containing gene and ID (.v) and the other containing the rest of the information (.w)
    :param plink: the directory of plink, if not set in environment
    :param bed: test file for genotype.
    :param bim: file with variant information
    :param fam: text file for pedigree information
    :return:
    """
    click.echo('staring plink processing ...')
    plink_process(genes_folder=genes_folder, plink=plink, bed=bed, bim=bim, fam=fam)
    click.echo('plink processing is complete.')


@main.command()
@click.option('-s', '--scores-file', required=True, help="The scoring file of genes across a population.")
@click.option('-i', '--genotype-file', required=True, help="File containing information about the cohort.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('-g', '--genes',
              help="a list containing the genes to calculate. if not provided all genes will be used.")
@click.option('-t', '--test', required=True,
              type=click.Choice(['ttest_ind', 'mannwhitneyu', 'logit', 'glm', 'betareg']),
              help='statistical test for calculating P value.')
@click.option('-c', '--cases-column', required=True, help="the name of the column that contains the case/control type.")
@click.option('-m', '--samples-column', required=True, help="the name of the column that contains the samples.")
@click.option('-p', '--pc-file', default=None,
              help="Principle components values for logistic regression. if not in genotype file")
@click.option('--adj-pval', type=click.Choice(
    ['bonferroni', 'sidak', 'holm-sidak', 'holm',
     'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']))
@click.option('--covariates', default='PC1,PC2', help="the covariates used for calculation")
def calculate_pval(
    *,
    scores_file,
    genotype_file,
    output_path,
    genes,
    cases_column,
    samples_column,
    test,
    pc_file,
    adj_pval,
    covariates,
):
    """
    Calculate the P-value between two given groups.
    :param scores_file: the file containing gene scores.
    :param genotype_file: file containing the phenotype.
    :param output_path: the path for final output.
    :param genes: a list of genes to calculate. if not, all genes in scoring file will be used.
    :param cases_column: the name of the column with phenotypes.
    :param samples_column: the name of the column with sample IDs. All files need to have the same format.
    :param test: the test used to calculate pvalue.
    :param pc_file: the file with PC (alternatively the file with covariates to use in test).
    :param adj_pval: the adjustment method used (if any).
    :param covariates: the column names of covariates to use, with comma in between. (e.g: PC1,PC2,age)
    :return:
    """
    if test == 'betareg':
        betareg_pvalues(
            scores_file=scores_file,
            pheno_file=genotype_file,
            cases_col=cases_column,
            samples_col=samples_column,
            output_path=output_path,
            covariates=covariates
        )
    else:
        scores_df = pd.read_csv(scores_file, sep=r'\s+')

        click.echo("The process for calculating the p_values will start now.")
        df = find_pvalue(
            scores_df=scores_df,
            output_file=output_path,
            genotype_file=genotype_file,
            genes=genes,
            cases_column=cases_column,
            samples_column=samples_column,
            test=test,
            pc_file=pc_file,
            adj_pval=adj_pval,
        )
        click.echo('Process is complete.')
        click.echo(df.info())


@main.command()
@click.option('-i', '--input-path', required=True, help="The directory that contains the matrices to merge.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('-r', '--remove-input', is_flag=True, help='if flagged will remove input folder')
def merge(
    *,
    output_path,
    input_path,
    remove_input,
):
    """
    This command merges all matrices in a directory into one big matrix.
    :param output_path: the path of final merged matrix
    :param input_path: the directory containing the matrices to merge.
    :param remove_input: if True, the input directory will be removed after merge.
    :return:
    """
    click.echo("Starting the merging process")
    df = combine_scores(input_path=input_path, output_path=output_path)
    click.echo(df.info())
    if remove_input:
        shutil.rmtree(input_path)
    click.echo("Merging is done.")


@main.command()
@click.option('-m', '--matrix-file', required=True, help="The scoring matrix to normalize.")
@click.option('-g', '--genes-lengths-file',
              help="The file containing the lengths of genes. If not provided it will be produced.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('-s', '--samples-col', default='IID', help='the name of the samples column')
def normalize(
    *,
    matrix_file,
    genes_lengths_file=None,
    output_path=None,
    samples_col
):
    """This command normalizes the scoring matrix by gene length."""
    click.echo("Normalization in process.")
    normalize_gene_len(
        matrix_file=matrix_file,
        genes_lengths_file=genes_lengths_file,
        output_path=output_path,
        samples_col=samples_col
    )


@main.command()
@click.option('--first-file', required=True, help="the path to the first scores file.")
@click.option('--second-file', required=True, help="the path to the second scores file.")
@click.option('--samples-col', default='IID', help="the column containing the samples IDs.")
@click.option('--output-file', required=True, help="the path to the output file with correlation values.")
def calc_corr(
    *,
    first_file,
    second_file,
    samples_col,
    output_file,
):
    """
    Calculate the pearson's correlation between same genes in two scoring matices.
    :param first_file: the path to the first scores file.
    :param second_file: the path to the second scores file.
    :param samples_col: the column containing the samples IDs.
    :param output_file: the path to the output file with correlation values.
    :return:
    """
    with open(first_file) as f:
        genes_01 = re.split('\s+', f.readline().strip('\n'))
        genes_01.remove(samples_col)
    with open(second_file) as f:
        genes_02 = re.split('\s+', f.readline().strip('\n'))
        genes_02.remove(samples_col)
    as_set = set(genes_01)
    common_genes = as_set.intersection(genes_02)
    genes = list(common_genes)
    corr_info = []
    first_df = pd.read_csv(first_file, sep=r'\s+', index_col=False)
    second_df = pd.read_csv(second_file, sep=r'\s+', index_col=False)
    for gene in tqdm(genes, desc='calculating correlation'):
        gene_df = pd.merge(first_df[[samples_col, gene]], second_df[[samples_col, gene]], on=samples_col)
        gene_df.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
        corr, pval = pearsonr(gene_df[gene + '_x'], gene_df[gene + '_y'])
        corr_info.append([gene, corr, pval])
    corr_df = pd.DataFrame(corr_info, columns=['genes', 'corr', 'p_value']).sort_values(by=['p_value'])
    corr_df.to_csv(output_file, sep='\t', index=False)
    click.echo('Process is complete.')
    click.echo(corr_df.info())


@main.command()
@click.option('--pvals-file', required=True, help="the file containing p-values.")
@click.option('--info-file', required=True, help="file containing variant/gene info.")
@click.option('--genescol-1', default='gene', help="the name of the genes column in pvals file.")
@click.option('--genescol-2', default='Gene.refGene', help="the name of the genes column in info file.")
@click.option('--qq-output', required=True, help="the name of the qq plot file.")
@click.option('--manhattan-output', required=True, help="the name of the manhatten plot file.")
@click.option('--pvalcol', default='p_value', help="the name of the pvalues column.")
def visualize(
    pvals_file,
    info_file,
    genescol_1,
    genescol_2,
    qq_output,
    manhattan_output,
    pvalcol,
):
    """
    Visualize manhatten plot and qqplot for the data.
    :param pvals_file: the file containing p-values.
    :param info_file: file containing variant/gene info.
    :param genescol_1: the name of the genes column in pvals file.
    :param genescol_2: the name of the genes column in info file.
    :param qq_output: the name of the qq plot file.
    :param manhattan_output: the name of the manhatten plot file.
    :param pvalcol: the name of the pvalues column.
    :return:
    """
    r_visualize(
        pvals_file=pvals_file,
        info_file=info_file,
        genescol_1=genescol_1,
        genescol_2=genescol_2,
        qq_output=qq_output,
        manhattan_output=manhattan_output,
        pvalcol=pvalcol,
    )


@main.command()
@click.option('--data-file', required=True, help='file with all features and target for training model.')
@click.option('--output-folder', required=True, help='path of folder that will contain all outputs.')
@click.option('--test-size', default=0.25, help='test size for cross validation and evaluation.')
@click.option('--test', is_flag=True,
              help='if flagged, a test set will be created for evaluating the final model.')
@click.option('--model-name', required=True, help='name of model file.')
@click.option('--model-type', required=True, type=click.Choice(['reg', 'classifier']), help='type of prediction model.')
@click.option('--target-col', required=True, help='name of target column in data_file.')
@click.option('--imbalanced', is_flag=True, help='if flagged methods will be used to account for the imbalance.')
@click.option('--normalize', is_flag=True, help='if flagged the data will be normalized before training.')
@click.option('--folds', default=10, type=int, help='number of cross-validation folds in training.')
@click.option('--metric', help='the metric used to choose best model after training.')
def prediction_model(
    *,
    data_file,
    output_folder,
    test_size=None,
    test,
    model_name,
    model_type,
    target_col,
    imbalanced,
    normalize,
    folds,
    metric,
):
    """
    Create a predicition model with given dataset.
    :param data_file: file containing features and target.
    :param output_folder: a folder path to save all outputs.
    :param test_size: the size of testing set.
    :param test: if True the dataset will be split into training and testing for extra evaluation after finalization.
    :param model_name: the name of the model to be saved.
    :param model_type: the type of model (reg or classifier).
    :param target_col: the column of the target.
    :param imbalanced: if true methods will be used to account for the imbalance.
    :param normalize:  if true the data will be normalized before training
    :param folds: the number of folds used for cross validation
    :param metric: the metric used to choose best model after training.
    :return: the final model
    """
    training_set = pd.read_csv(data_file, sep='\s+', index=False)
    if test:
        training_set, testing_set = train_test_split(training_set, test_size=test_size)
    else:
        testing_set = None
    os.mkdir(output_folder)
    os.chdir(output_folder)
    model = create_model(
        model_name=model_name,
        model_type=model_type,
        imbalanced=imbalanced,
        normalize=normalize,
        folds=folds,
        metric=metric,
        y_col=target_col,
        training_set=training_set,
        testing_set=testing_set,
        test_size=test_size,
    )
    return model


if __name__ == '__main__':
    main()
