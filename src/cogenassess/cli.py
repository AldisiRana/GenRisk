# -*- coding: utf-8 -*-
import os

import click
import shutil

import pandas as pd
from sklearn.model_selection import train_test_split

from .utils import get_gene_info, plink_process, combine_scores
from .pipeline import find_pvalue, betareg_pvalues, r_visualize, create_prediction_model


@click.group()
def main():
    """Handle cogenassess functions."""


@main.command()
@click.option('-a', '--annotated-vcf', required=True, help='the annotated vcf')
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
    annotated_vcf,
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
    Calculate the gene-based scores for a given dataset.
    :param annotated_vcf: an annotated containing variant IDs, alt, info and samples genotypes.
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
        annotated_vcf=annotated_vcf,
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
    plink_process(genes_folder=genes_folder, plink=plink, annotated_vcf=annotated_vcf)
    click.echo('combining score files ...')
    df = combine_scores(input_path=temp_dir, output_path=output_file)
    if remove_temp:
        shutil.rmtree(temp_dir)
    click.echo('process is complete.')
    return df.info()


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
@click.option('--model-type', required=True, type=click.Choice(['regressor', 'classifier']),
              help='type of prediction model.')
@click.option('--target-col', required=True, help='name of target column in data_file.')
@click.option('--imbalanced', is_flag=True, help='if flagged methods will be used to account for the imbalance.')
@click.option('--normalize', is_flag=True, help='if flagged the data will be normalized before training.')
@click.option('--folds', default=10, type=int, help='number of cross-validation folds in training.')
@click.option('--metric', help='the metric used to choose best model after training.')
def create_model(
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
    Create a machine learning model with given dataset.
    :param data_file: file containing features and target.
    :param output_folder: a folder path to save all outputs.
    :param test_size: the size of testing set.
    :param test: if True the dataset will be split into training and testing for extra evaluation after finalization.
    :param model_name: the name of the model to be saved.
    :param model_type: the type of model (regressor or classifier).
    :param target_col: the name of the target column in data file.
    :param imbalanced: if true methods will be used to account for the imbalance.
    :param normalize:  if true the data will be normalized before training
    :param folds: the number of folds used for cross validation
    :param metric: the metric used to choose best model after training.
    :return: the final model
    """
    training_set = pd.read_csv(data_file, sep='\s+', index_col=False)
    testing_set = pd.DataFrame()
    if test:
        training_set, testing_set = train_test_split(training_set, test_size=test_size)
    os.mkdir(output_folder)
    os.chdir(output_folder)
    model = create_prediction_model(
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
