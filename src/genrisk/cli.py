# -*- coding: utf-8 -*-
import os
import random
import shutil
import time

import click
import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split

from .gene_scoring import calculate_gbrs, pathway_scoring
from .helpers import create_logger
from .pipeline import find_pvalue, betareg_pvalues, create_prediction_model, model_testing, scoring_process, \
    normalize_data
from .prs_scoring import prs_prompt
from .utils import draw_qqplot, draw_manhattan, merge_files

SAMPLES_COL = click.option('-m', '--samples-col',  default='IID', show_default=True,
                           help="the name of the column that contains the samples.")
OUTPUT_FILE = click.option('-o', '--output-file', required=True, help="the final output path")
logger = create_logger()


@click.group()
def main():
    """Handle genrisk functions."""


@main.command()
@click.option('-a', '--annotated-vcf', required=True, type=click.Path(exists=True), help='an annotated containing variant IDs, alt, info and samples genotypes.')
@click.option('-b', '--bfiles', default=None,
              help='provide binary files if annotated vcf does not contain the samples info')
@click.option('--plink', default='plink', help="the directory of plink, if not set in environment")
@click.option('-t', '--temp-dir', required=True, help="a temporary directory to save temporary files before merging.")
@OUTPUT_FILE
@click.option('-p', '--beta-param', show_default=True, default=(1.0, 25.0), nargs=2, type=float,
              help="the parameters from beta weight function.")
@click.option('-w', '--weight-func', show_default=True, default='beta', type=click.Choice(['beta', 'log10']),
              help="the weighting function used in score calculation.")
@click.option('-v', '--variant-col', show_default=True, default='SNP', help="the column containing the variant IDs.")
@click.option('-g', '--gene-col', show_default=True, default='Gene.refGene', help="the column containing gene names.")
@click.option('-f', '--af-col', show_default=True, default='MAF', help="the column containing allele frequency.")
@click.option('-d', '--del-col', show_default=True, default='CADD_raw', help="the column containing the deleteriousness score.")
@click.option('-l', '--alt-col', show_default=True, default='Alt', help="the column containing the alternate base.")
@click.option('-m', '--maf-threshold', show_default=True, default=0.01, help="the threshold for minor allele frequency.")
def score_genes(
        *,
        annotated_vcf,
        bfiles,
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
):
    """
    Calculate the gene-based scores for a given dataset.

    Example
    ---------
    ::

        $ genrisk score-genes --annotated-vcf annotated_vcf_toy.vcf --temp-dir test/
        --output-file test.tsv --weight-func beta --maf-threshold 0.01 --alt-col ALT
        --variant-col ID --af-col AF --del-col CADD --gene-col Gene

    \f

    Parameters
    ----------
    annotated_vcf : str
        an annotated containing variant IDs, alt, info and samples genotypes.
    bfiles : str
        the binary files for plink process.
        this arg is not needed if the annotated vcf contains all information.
    plink : str
        the location of plink, if not set in environment
    beta_param : tuple
        the parameters from beta weight function.
    temp_dir : str
        a temporary directory to save temporary files before merging.
    output_file : str
        the location and name of the final output scores matrix.
    weight_func : str
        the weighting function used on allele frequency in score calculation. [beta| log10]
    variant_col : str
        the column containing the variant IDs.
    gene_col : str
        the column containing gene names.
        If the genes are in the INFO column, use the identifier of the value (i.e gene=IF, identifier is 'gene')
    af_col : str
        the column containing allele frequency. If in INFO, follow previous example
    del_col : str
        the column containing deleteriousness score (functional annotation). If in INFO, follow previous example
    alt_col : str
        the column containing alternate base.
    maf_threshold : float
        the threshold for minor allele frequency.

    Returns
    -------
    DataFrame information
        the final scores dataframe information
        the DataFrame is saved into the output path indicated in the arguments

    """
    confirm = click.confirm('Would you like us to delete the temporary files when process is done?')
    logger.info('Gene-based scoring')
    logger.info('Score genes process is starting now...')
    logger.info(locals())
    logger.info('getting information from vcf files')
    start_time = time.time()
    df = scoring_process(
        logger=logger,
        annotated_vcf=annotated_vcf,
        bfiles=bfiles,
        plink=plink,
        beta_param=beta_param,
        temp_dir=temp_dir,
        output_file=output_file,
        weight_func=weight_func,
        variant_col=variant_col,
        gene_col=gene_col,
        af_col=af_col,
        del_col=del_col,
        alt_col=alt_col,
        maf_threshold=maf_threshold,
    )
    if confirm:
        logger.info('The temporary files will be removed now.')
        shutil.rmtree(temp_dir)
    time.sleep(1)
    end_time = time.time()
    logger.info(f"Runtime of the program is {end_time - start_time} Sec")
    return df.info()


@main.command()
@click.option('-s', '--scores-file', required=True, type=click.Path(exists=True),
              help="The scoring file of genes across a population.")
@click.option('-i', '--info-file', required=True, type=click.Path(exists=True),
              help="File containing information about the cohort.")
@OUTPUT_FILE
@click.option('-g', '--genes',
              help="a file containing the genes to calculate. if not provided all genes will be used.")
@click.option('-t', '--test', required=True,
              type=click.Choice(['ttest_ind', 'mannwhitneyu', 'logit', 'betareg', 'linear']),
              help='statistical test for calculating P value.')
@click.option('-c', '--cases-col', required=True,
              help="the name of the column that contains the case/control or quantitative vals.")
@SAMPLES_COL
@click.option('-a', '--adj-pval', type=click.Choice(
    ['bonferroni', 'sidak', 'holm-sidak', 'holm',
     'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']))
@click.option('-v', '--covariates', default='', help="the covariates used for calculation")
@click.option('-p', '--processes', show_default=True, type=int, default=1, help='number of processes for parallelization')
def find_association(
        *,
        scores_file,
        info_file,
        output_file,
        genes,
        cases_col,
        samples_col,
        test,
        adj_pval,
        covariates,
        processes,
):
    """
    Calculate the P-value between two given groups.

    Example
    ---------
    ::

        $ genrisk find-association --scores-file toy_example/toy_dataset_scores --info-file
        toy_example/toy.pheno --cases-column trait1 --samples-column IID --test betareg --output-file
        toy_dataset_betareg.tsv --covariates age,sex --adj-pval bonferroni
    \f

    Parameters
    ----------
    scores_file : str
        the file containing gene-based scores.
    info_file : str
        file containing the phenotype.
    output_file : str
        path to the final output.
    genes : str
        a file that contains a list of genes to calculate p-values. if not, all genes in scoring file will be used.
    cases_col : str
        the name of the column with phenotypes. Phenotypes can be either binary or quantitative.
    samples_col : str
         the name of the column with sample IDs. All files need to have the same format.
    test : str
        the statistical test used for calculating p-values.
    adj_pval : str, optional
        the method used to adjust the p-values.
    covariates : str, optional
        the covariates used for calculation. Not all tests are able to include covariates.
        (e.g. Mann Whinteny U doesn't allow for covariates)
    processes : int, optional
        if more than 1 processer is selected, the function will be parallelized.

    Returns
    -------
    DataFrame information
        the final dataframe information
        the DataFrame is saved into the output path indicated in the arguments

    """
    logger.info('Finding associations')
    logger.info(locals())
    logger.info("The process for calculating the p_values will start now.")
    start_time = time.time()
    if test == 'betareg':
        betareg_pvalues(
            scores_file=scores_file,
            pheno_file=info_file,
            cases_col=cases_col,
            samples_col=samples_col,
            output_path=output_file,
            covariates=covariates,
            processes=processes,
            genes=genes,
            logger=logger,
        )
        end_time = time.time()
        logger.info(f"Runtime of the program is {end_time - start_time}")
    else:
        if genes:
            with open(genes) as f:
                content = f.readlines()
            genes = [x.strip() for x in content]
        df = find_pvalue(
            scores_file=scores_file,
            output_file=output_file,
            info_file=info_file,
            genes=genes,
            cases_column=cases_col,
            samples_column=samples_col,
            test=test,
            adj_pval=adj_pval,
            covariates=covariates,
            processes=processes,
            logger=logger
        )
        end_time = time.time()
        logger.info(f"Runtime of the program is {end_time - start_time}")
        return df.info()


@main.command()
@click.option('-p', '--pvals-file', required=True, type=click.Path(exists=True), help="the file containing p-values.")
@click.option('-i', '--info-file', type=click.Path(exists=True), help="file containing variant/gene info.")
@click.option('--genescol-1', show_default=True, default='gene', help="the name of the genes column in pvals file.")
@click.option('--genescol-2', show_default=True, default='Gene.refGene', help="the name of the genes column in info file.")
@click.option('-q', '--qq-output', default=None, help="the name of the qq plot file.")
@click.option('-m', '--manhattan-output', default=None, help="the name of the manhatten plot file.")
@click.option('-v', '--pval-col', show_default=True, default='p_value', help="the name of the pvalues column.")
@click.option('-c', '--chr-col', show_default=True, default='Chr', help='the name of the chromosomes column')
@click.option('-s', '--pos-col', show_default=True, default='Start', help='the name of the position/start of the gene column')
def visualize(
        *,
        pvals_file,
        info_file,
        genescol_1,
        genescol_2,
        qq_output,
        manhattan_output,
        pval_col,
        chr_col,
        pos_col,
):
    """
    Visualize manhatten plot and qqplot for the data.

    Example
    --------
    ::

        $ genrisk visualize --pvals-file toy_example/toy_dataset_scores
        --info-file annotated_toy_dataset.vcf --qq-output toy_example/toy_dataset_qqplot.jpg
        --manhattan-output toy_example/toy_dataset_manhattanplot.jpg
    \f

    Parameters
    ----------
    pvals_file : str
        the file containing the calculated p-values.
    info_file : str
        file containing variant/gene info.
    genescol_1 : str
        the name of the genes column in pvals file.
    genescol_2 : str
        the name of the genes column in info file.
    qq_output : str
        the name of the qq plot file. If left empty no file will be produced.
    manhattan_output : str
        the name of the manhatten plot file. If left empty no file will be produced
    pval_col : str
        the name of the pvalues column.
    chr_col : str
        the name of chromosomes column.
    pos_col : str
        the name of the position/start column.

    Returns
    -------

    """
    logger.info('Creating plots for data')
    logger.info(locals())
    logger.info('Reading p_values file...')
    pvals_df = pd.read_csv(pvals_file, sep='\t', index_col=False)
    if qq_output:
        logger.info('Creating QQ-plot...')
        try:
            draw_qqplot(pvals=pvals_df[pval_col], qq_output=qq_output)
        except Exception as arg:
            logger.exception(arg)
            raise
    if manhattan_output:
        logger.info('Creating Manhattan plot...')
        if not info_file:
            logger.exception('Please provide a file with gene information to generate manhattan plot.')
        info_df = pd.read_csv(info_file, sep="\t", index_col=False)
        merged = pd.merge(pvals_df, info_df, left_on=genescol_1, right_on=genescol_2, how='left')
        try:
            draw_manhattan(
                data=merged,
                chr_col=chr_col,
                pos_col=pos_col,
                pvals_col=pval_col,
                genes_col=genescol_1,
                manhattan_output=manhattan_output
            )
        except Exception as arg:
            logger.exception(arg)
            raise


@main.command()
@click.option('-d', '--data-file', type=click.Path(exists=True), required=True,
              help='file with all features and target for training model.')
@click.option('-o', '--output-folder', required=True, help='path of folder that will contain all outputs.')
@click.option('-i', '--test-size', show_default=True, default=0.25, help='test size for cross validation and evaluation.')
@click.option('-t', '--test', is_flag=True,
              help='if flagged, a test set will be created for evaluating the final model.')
@click.option('-n', '--model-name', required=True, help='name of model file.')
@click.option('--model-type', required=True, type=click.Choice(['regressor', 'classifier']),
              help='type of prediction model.')
@click.option('-l', '--target-col', required=True, help='name of target column in data_file.')
@click.option('-b', '--imbalanced', is_flag=True, help='if flagged methods will be used to account for the imbalance.')
@click.option('--normalize', is_flag=True, help='if flagged the data will be normalized before training.')
@click.option('--normalize-method', show_default=True, default='zscore', type=click.Choice(['zscore', 'minmax', 'maxabs', 'robust']),
              help='features normalization method.')
@click.option('-f', '--folds', show_default=True, default=10, type=int, help='number of cross-validation folds in training.')
@click.option('--metric', help='the metric used to choose best model after training.')
@SAMPLES_COL
@click.option('--seed', default=random.randint(1, 2147483647),
              help='add number to create reproduciple train_test splitting.')
@click.option('--include-models', default=None,
              help='choose specific models to compare with comma in between. e.g lr,gbr,dt')
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
        normalize_method,
        folds,
        metric,
        samples_col,
        seed,
        include_models,
):
    """
    Create a prediction model with given dataset.

    Example
    --------
     ::

        $ genrisk create-model --data-file toy_example_regressor_features.tsv --model-type regressor
        --output-folder toy_regressor  --test-size 0.25 --test --model-name toy_regressor
        --target-col trait1 --imbalanced --normalize

    Notes
    -------
    The types of models available for training can be found :ref:`model_types`
    \f

    Parameters
    ----------
    data_file : str
        file containing features and target.
    output_folder : str
        a folder path to save all outputs.
    test_size : float
        the size of testing set.
    test : bool
        if True the dataset will be split into training and testing for extra evaluation after finalization.
    model_name : str
        the name of the model to be saved.
    model_type : str
        the type of model [regressor| classifier].
    target_col : str
        the name of the target column in data file.
    imbalanced : bool
        if true methods will be used to account for the imbalance.
    normalize : bool
        if true the data will be normalized before training
    normalize_method : str
        method used to normalize data. [zscore| minmax| maxab| robust]
    folds : int
        the number of folds used for cross validation
    metric : str
        the metric used to choose best model after training.
    samples_col : str
        the name of the column with samples IDs.
    seed : int
        random seed number to run the machine learning models.
    include_models : str
        list of specific models to compare. more information in the documentations

    Returns
    -------
    Final prediction model

    """
    logger.info('Create prediction model')
    logger.info(locals())
    logger.info('Reading and preparing data ...')
    start_time = time.time()
    training_set = pd.read_csv(data_file, sep='\t', index_col=samples_col)
    training_set.dropna(subset=[target_col], inplace=True)
    training_set.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
    testing_set = pd.DataFrame()
    if test:
        training_set, testing_set = train_test_split(training_set, test_size=test_size, random_state=int(seed))
    os.mkdir(output_folder)
    os.chdir(output_folder)
    logger.info('Model generation starting ...')
    if include_models:
        include_models = include_models.split(',')
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
        seed=int(seed),
        include_models=include_models,
        normalize_method=normalize_method,
    )
    logger.info('Model is generated.')
    end_time = time.time()
    logger.info(f"Runtime of the program is {end_time - start_time}")
    return model


@main.command()
@click.option('-t', '--model-type', required=True, type=click.Choice(['regressor', 'classifier']),
              help='type of prediction model.')
@click.option('-i', '--input-file', required=True, type=click.Path(exists=True), help='testing dataset')
@click.option('-l', '--label-col', required=True, help='the target/phenotype/label column')
@click.option('-m', '--model-path', required=True, type=click.Path(exists=True), help='path to the trained model.')
@click.option('-s', '--samples-col', show_default=True, default='IID', help='the samples column.')
@OUTPUT_FILE
def test_model(
        *,
        model_path,
        input_file,
        model_type,
        label_col,
        samples_col,
        output_file,
):
    """
    Evaluate a prediction model with a given dataset.

    Example
    --------
     ::

        $ genrisk test-model --model-path regressor_model.pkl --input-file testing_dataset.tsv
        --model-type regressor --labels-col target --samples-col IID
    \f

    Parameters
    ----------
    model_path : str
        the path to the ML model.
    input_file : str
        the testing (independent) dataset.
    model_type : str
        the type of model [classifier|regressor].
    label_col : str
        the labels/target column.
    samples_col : str
        the sample ids column.
    output_file : str
        the path to the dataframe with the prediction results.

    Returns
    -------
    DataFrame
        dataframe with the prediction results.

    """
    logger.info('Testing prediction model')
    logger.info(locals())
    testing_df = model_testing(
        model_path=model_path,
        input_file=input_file,
        label_col=label_col,
        samples_col=samples_col,
        model_type=model_type,
    )
    logger.info('saving test predictions')
    testing_df.to_csv(output_file, sep='\t')
    return testing_df


@main.command()
@click.option('-p', '--plink', default='plink')
def get_prs(
        *,
        plink,
):
    """
    Calculate PRS. This command is interactive.
    This command gets a pgs file (provided by the user or downloaded) then calculates the PRS for dataset.

    Example
    -------
    This function is performed using commandline interface::

        $ genrisk get-prs
    \f

    Parameters
    ----------
    plink : str
        provide plink path if not default in environment.

    Returns
    -------

    """

    download = click.confirm('Do you want to download PGS file?')
    return prs_prompt(plink=plink, download=download)


@main.command()
@click.option('-f', '--files', required=True,
              help='input all files to merge with a comma in between. E.g: file1,file2,file3')
@click.option('-s', '--sep', show_default=True, default='\t', help='the column seperator in files.')
@click.option('-b', '--by', show_default=True, default='IID', help='the common column between all files to merge.')
@click.option('-c', '--cols', default=None,
              help='if desired, a list of columns can be chosen to save final file, e.g: col1,col2,col5')
@OUTPUT_FILE
def merge(
        *,
        files,
        sep,
        by,
        cols,
        output_file
):
    """
    Merge all files given into one dataframe. \f

    Parameters
    ----------
    files : str
         all files to merge with a comma in between. E.g: file1,file2,file3
    sep : str
        the column seperator in files.
    by : str
        the common column between all files to merge.
    cols : str
        a list of columns can be chosen to save final file, e.g: col1,col2,col5. if empty all columns will be merged
    output_file : str
        the name and path to save final dataframe

    Returns
    -------
    DataFrame
        final merged dataframe.

    """
    logger.info('GenRisk - merging dataframes')
    logger.info(locals())
    files_lst = files.split(',')
    df = merge_files(
        files_lst=files_lst,
        sep=sep,
        by=by,
        cols=cols
    )
    logger.info('saving merged dataframe.')
    df.to_csv(output_file, sep=sep)
    return df


@main.command()
@click.option('-s', '--scores-file', required=True, type=click.Path(exists=True), help='the gene-based scores file.')
@click.option('-w', '--weights-file', default=None, type=click.Path(exists=True),
              help='a file containg the weights for each gene.')
@click.option('-p', '--pheno-file', default=None, type=click.Path(exists=True),
              help='if no weights are given, pheno file is used to calculate them.')
@click.option('-c', '--pheno-col', default=None, help='the column containing the phenotype.')
@click.option('-v', '--covariates', show_default=True, default='sex,age,bmi,PC1,PC2,PC3,PC4',
              help='the covariates to use in the linear model.')
@click.option('-g', '--genes-col', show_default=True, default='genes', help='the column containing gene names.')
@click.option('-e', '--weights-col', show_default=True, default='zscore', help='the name and path to output results.')
@OUTPUT_FILE
@click.option('--split-size', show_default=True, default=0.25, help='the size ratio to split dataset for weight calculation.')
@click.option('--method', type=click.Choice(['sum', 'pathways']),
              help='method for presenting the risk scores.')
@SAMPLES_COL
@click.option('--pathway-file', required=True, help='.gmt file containing the pathway and its genes.')
def get_gbrs(
        *,
        method,
        scores_file,
        weights_file,
        pheno_file,
        pheno_col,
        covariates,
        genes_col,
        samples_col,
        weights_col,
        output_file,
        split_size,
        pathway_file
):
    """
    Calculate gene-based risk score for individuals. 
    \f

    Parameters
    ----------
    method : str
        for presenting the risk score, either sum of all genes or pathways.
    scores_file : str
        the gene-based scores file
    weights_file : str
        a file containing the weight for each gene.
    pheno_file : str
        if no weights are given, pheno file is used to calculate them
    pheno_col : str
        the column containing the phenotype.
    covariates : str
        the covariates to use in the linear model.
    genes_col : str
        column containing the gene names.
    samples_col : str
        column containing the samples ids.
    weights_col : str
        column containing the weights for each gene (in the weights file)
    output_file : str
        the name and path to output results.
    split_size : str
        the proportion of the dataset to use for calculating the weights. should be between 0.0 and 1.0.
    pathway_file : str
        if pathways method is chosen, a file with the pathways is needed

    Returns
    -------
    DataFrame
        df contianing the final risk scores.

    """
    logger.info('GenRisk - Calculating gene-based risk scores')
    logger.info(locals())
    logger.info("Reading files...")
    scores_df = pd.read_csv(scores_file, sep=r'\s+', index_col=False)
    if weights_file:
        weights_df = pd.read_csv(weights_file, sep='\t')
    else:
        logger.info("No weights available, weights will be caluclated now.")
        logger.info("Creating temporary files...")
        scores_df, scores_temp = train_test_split(scores_df, test_size=split_size, random_state=0)
        scores_temp.to_csv('scores_temp.tsv', sep='\t', index=False)
        logger.info("The process for calculating the p_values will start now.")
        weights_df = find_pvalue(
            scores_file='scores_temp.tsv',
            info_file=pheno_file,
            output_file='association_' + output_file,
            cases_column=pheno_col,
            samples_column=samples_col,
            test='linear',
            covariates=covariates,
            logger=logger,
        )
        logger.info("Exluding samples used in weights calculation.")
        scores_df.reset_index(drop=True, inplace=True)
        logger.info("Remove temporary files.")
        del scores_temp
        os.remove('scores_temp.tsv')
        weights_df['zscore'] = weights_df['beta_coef'] / weights_df['std_err']
    logger.info("Calculating GBRS now ...")
    weights_df[weights_col] = weights_df[weights_df[weights_col].apply(lambda x: x.isnumeric())]
    weights_df[weights_col] = pd.to_numeric(weights_df[weights_col])
    df = calculate_gbrs(
        scores_df=scores_df,
        weights_df=weights_df,
        weights_col=weights_col,
        genes_col=genes_col,
        method=method,
        pathway_file=pathway_file,
        samples_col=samples_col,
        logger=logger,
    )
    df[samples_col] = scores_df[samples_col]
    logger.info("GBRS dataframe is being saved ...")
    df.to_csv(output_file, sep='\t', index=False)
    logger.info("Process is complete.")
    return df


@main.command()
@OUTPUT_FILE
@SAMPLES_COL
@click.option('-p', '--pathway-file', required=True, type=click.Path(exists=True),
              help='.gmt file containing the pathway and its genes.')
@click.option('-s', '--scores-file', required=True, type=click.Path(exists=True),
              help='genes scores file to calculate pathway scores.')
def calculate_pathways(
        *,
        output_file,
        pathway_file,
        scores_file,
        samples_col
):
    """
    Calculate pathway scores using gene-based scores and gmt pathway file. \f

    Parameters
    ----------
    output_file : str
        the path to final output.
    pathway_file : str
        .gmt file containing the pathway and its genes.
    scores_file : str
        genes scores file to calculate pathway scores.
    samples_col : str
        column containing samples ids.

    Returns
    -------
    DataFrame
        the pathway scores df

    """
    logger.info('GenRisk - calculating pathway scores')
    logger.info(locals())
    pathways = {line.strip().split('\t')[0]: line.strip().split('\t')[2:] for line in open(pathway_file, 'r')}
    all_genes = [item for sublist in list(pathways.values()) for item in sublist]
    scored_genes = open(scores_file).readline().rstrip().split()
    combined_genes = list(set(all_genes) & set(scored_genes))
    df = pathway_scoring(pathways=pathways, genes=combined_genes, scores_file=scores_file, samples_col=samples_col,
                         logger=logger)
    df.to_csv(output_file, sep='\t', index=False)
    logger.info('Process is done.')
    return df


@main.command()
@click.option('--method', required=True, type=click.Choice(['gene_length', 'zscore', 'minmax','maxabs', 'robust']))
@click.option('--data-file', required=True, type=click.Path(exists=True))
@click.option('--genes-info', default=None)
@SAMPLES_COL
@click.option('--genes-col', show_default=True, default='HGNC symbol', type=str)
@click.option('--lengths-col', show_default=True, default='gene_length', type=str)
@OUTPUT_FILE
def normalize(
        *,
        method,
        data_file,
        genes_info,
        samples_col,
        genes_col,
        lengths_col,
        output_file,

):
    """
    Normalize/standarize data.

    Example
    -------
    ::

        $ genrisk normalize --data-file toy_example/toy_dataset_scores --method gene_length --samples-col IID
        --output-file toy_dataset_scores_normalized.tsv

    \f

    Parameters
    ----------
    genes_info : str
        the file containing genes names and length. if not provided ensembl database is used to retrieve data.
    method : str
        the method of normalizing data. [gene_length|zscore|minmax|maxabs|robust]
    data_file : str
        the file containg data to be normalized.
    samples_col : str
        the column containing sample ids.
    genes_col : str
        the column containing gene names. ignore if genes_info file is not provided.
    lengths_col : str
        the column containing gene lengths. ignore if genes_info file is not provided.
    output_file : str
        the name of the file for final output

    Returns
    -------
    DataFrame with normalized data.

    """
    logger.info('GenRisk - normalizing data using '+ method)
    logger.info(locals())
    df = normalize_data(
        method=method, data_file=data_file, samples_col=samples_col, genes_col=genes_col, length_col=lengths_col,
        genes_info=genes_info,
    )
    df.to_csv(output_file, sep='\t', index=False)
    logger.info('Process is done.')
    return df

if __name__ == '__main__':
    main()
