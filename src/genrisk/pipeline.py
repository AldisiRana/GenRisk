# -*- coding: utf-8 -*-
import os
import subprocess

import joblib
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .association_analysis import run_mannwhitneyu, run_ttest, get_pvals_logit, get_pvals_linear
from .gene_scoring import get_gene_info, plink_process, combine_scores
from .prediction_models import regression_model, classification_model, test_classifier, test_regressor
from .utils import gene_length_normalize

PATH = os.path.abspath(os.path.join((__file__), os.pardir, os.pardir, os.pardir))
BETAREG_SHELL = os.path.join(PATH, 'scripts', 'betareg_shell.R')
association_functions = {
    'mannwhitneyu': run_mannwhitneyu,
    'ttest_ind': run_ttest,
    'logit': get_pvals_logit,
    'linear': get_pvals_linear,
}


def scoring_process(
    *,
    logger,
    annotated_vcf,
    temp_dir,
    beta_param,
    weight_func,
    del_col,
    maf_threshold,
    gene_col,
    variant_col,
    af_col,
    alt_col,
    bfiles,
    plink,
    output_file
):
    """
    Calculate gene-based scores.
    This is calculated by a weighted sum of the variants in a gene.

    Parameters
    ----------
    logger
        an object that logs function outputs.
    annotated_vcf : str
        an annotated containing variant IDs, alt, info and samples genotypes.
    temp_dir : str
        a temporary directory to save temporary files before merging.
    beta_param : tuple
        the parameters from beta weight function. ignore if log10 function is chosen.
    weight_func : str
        the weighting function used in score calculation.
    del_col : str
        the column containing deleteriousness score or functional annotation.
    maf_threshold : float
        the threshold for minor allele frequency. between [0.0-1.0]
    gene_col : str
        the column containing gene names.
    variant_col : str
        the column containing variant IDs.
    af_col : str
        the column containing allele frequency.
    alt_col : str
        the column containing alternate allele.
    bfiles : str
        the binary files for plink process. if vcf contains all information bfiles are not needed.
    plink : str
        the directory of plink, if not set in environment
    output_file : str
        the path to save the final output scores matrix.

    Returns
    -------
    DataFrame
        final scores matrix

    """
    try:
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
    except Exception as arg:
        logger.exception(arg)
        raise
    logger.info('calculating gene scores ...')
    try:
        plink_process(genes_folder=genes_folder, plink=plink, annotated_vcf=annotated_vcf, bfiles=bfiles)
    except Exception as arg:
        logger.exception(arg)
        raise
    logger.info('combining score files ...')
    try:
        df = combine_scores(input_path=temp_dir, output_path=output_file)
    except Exception as arg:
        logger.exception(arg)
        raise
    logger.info('process is complete.')
    return df


def find_pvalue(
    *,
    scores_file,
    info_file,
    output_file,
    genes=None,
    cases_column,
    samples_column,
    test='mannwhitneyu',
    adj_pval=None,
    covariates=None,
    cases=None,
    controls=None,
    processes=1,
    logger
):
    """
    Calculate the significance of a gene in a population using different statistical analyses [mannwhitneyu, logit, linear, ttest_ind]

    Parameters
    ----------
    scores_file : str
        dataframe containing the scores of genes across samples.
    info_file : str
        a file containing the information of the sample. this includes target phenotype and covariates.
    output_file : str
         a path to save the summary statistics.
    genes : list
        a list of the genes to calculate the significance. if None will calculate for all genes.
    cases_column : str
        the name of the column containing phenotype information.
    samples_column : str
        the name of the column contining samples IDs.
    test : str
        the type of statistical test to use, choices are: ttest_ind, mannwhitenyu, linear, logit.
    adj_pval : str
        the method for pvalue adjustment. Choices: ['bonferroni', 'sidak', 'holm-sidak', 'holm', 'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']
    covariates : str
        the list of covariates used in the calculation.
    cases : str
        the cases category. if binary phenotype.
    controls : str
        the controls category. if binary phenotype.
    processes : int
        number of processes used, for parallel computing.
    logger
        an object that logs function outputs.

    Returns
    -------
    DataFrame
        dataframe with genes and their p_values

    """
    logger.info("Reading scores file...")
    if genes:
        scores_df = pd.read_csv(scores_file, sep=r'\s+', index_col=samples_column,
                                usecols=lambda x: x in genes+[samples_column])
    else:
        scores_df = pd.read_csv(scores_file, sep=r'\s+', index_col=samples_column)
    scores_df.replace([np.inf, -np.inf], 0, inplace=True)
    scores_df.fillna(0, inplace=True)
    scores_df = scores_df.loc[:, scores_df.var() != 0.0].reset_index()
    logger.info("Reading info file...")
    if covariates:
        covariates = covariates.split(',')
        genotype_df = pd.read_csv(info_file, sep='\t', usecols=covariates+[cases_column, samples_column])
        genotype_df.fillna(genotype_df.mean(), inplace=True)
        genotype_df.dropna(inplace=True)
    else:
        genotype_df = pd.read_csv(info_file, sep='\t')
        genotype_df.dropna(subset=[cases_column], inplace=True)
    logger.info("Processing files...")
    merged_df = pd.merge(scores_df, genotype_df, how='inner', on=samples_column)
    merged_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    merged_df.apply(lambda x: x.fillna(x.mean()), axis=1)
    genes = scores_df.columns.tolist()[1:]
    del scores_df
    args = {
        'processes': processes, 'cases': cases, 'controls': controls, 'covariates': covariates, 'logger': logger,
    }
    logger.info("Calculating p_values using the following test: " + test)
    try:
        p_values_df = association_functions.get(test)(df=merged_df, genes=genes, cases_column=cases_column, **args)
        p_values_df.dropna(subset=['p_value'], inplace=True)
    except Exception as arg:
        logger.exception(arg)
        raise
    if adj_pval:
        logger.info("Calculating the adjusted p_values...")
        adjusted = multipletests(list(p_values_df['p_value']), method=adj_pval)
        p_values_df[adj_pval + '_adj_pval'] = list(adjusted)[1]
    p_values_df.to_csv(output_file, sep='\t', index=False)
    logger.info("Process is complete. The association analysis results have been saved.")
    return p_values_df


def betareg_pvalues(
    *,
    scores_file,
    pheno_file,
    samples_col,
    cases_col,
    output_path,
    covariates,
    processes,
    genes,
    logger,
):
    """
    Calculate association significance using betareg. This function runs in Rscript.


    Parameters
    ----------
    scores_file : str
        the path to the scores file.
    pheno_file : str
        the path to the phenotypes and covariates file.
    samples_col : str
        column containing samples ids.
    cases_col : str
        column containing the phenotype.
    output_path : str
        a path to save the summary statistics.
    covariates : str
        the list of covariates used in the calculation.
    processes : int
        number of processes used, for parallel computing.
    genes : str
        a list of the genes to calculate the significance. if None will calculate for all genes.
    logger
        an object that logs function outputs.

    Returns
    -------

    """
    logger.info("The p_values will be calculated using beta regression.")
    try:
        p = subprocess.call(
            ["Rscript", BETAREG_SHELL,
             "-s", scores_file,
             "--phenofile", pheno_file,
             "--samplescol", samples_col,
             "--casescol", cases_col,
             "-o", output_path,
             "--covariates", covariates,
             "--processes", str(processes),
             "--genes", genes]
        )
    except Exception as arg:
        logger.exception(arg)
        raise
    logger.info("Process is complete. The association analysis results have been saved.")


def create_prediction_model(
    *,
    model_name='final_model',
    model_type='regressor',
    y_col,
    imbalanced=True,
    normalize=True,
    folds=10,
    training_set,
    testing_set=pd.DataFrame(),
    test_size=0.25,
    metric=None,
    seed,
    include_models,
    normalize_method,
):
    """
    Create a prediction model (classifier or regressor) using the provided dataset.

    Parameters
    ----------
    model_name : str
        the name of the prediction model.
    model_type : str
        type of model [regressor|classifier]
    y_col : str
        the column containing the target (qualitative or quantitative).
    imbalanced : bool
        True means data is imbalanced.
    normalize : bool
        True if data needs normalization.
    folds : int
        how many folds for cross-validation.
    training_set : pd.DataFrame
        the training set for the model.
    testing_set : pd.DataFrame
        if exists an extra evaluation step will be done using the testing set.
    test_size : float
         the size to split the training set for cross-validation.
    metric : str
        the metric to evaluate the best model.
    seed : int
        the intilization state random number
    include_models : list
        a list of models that the user wants to test. if None all models will be used.
    normalize_method : str
        the method to normalize the data. Choices [zscore, minmax, maxabs, robust]

    Returns
    -------
    Final model

    """
    try:
        model_func = {'classifier': classification_model, 'regressor': regression_model}
        final_model = model_func.get(model_type)(
            y_col=y_col,
            training_set=training_set,
            normalize=normalize,
            test_size=test_size,
            folds=folds,
            metric=metric,
            model_name=model_name,
            testing_set=testing_set,
            imbalanced=imbalanced,
            seed=seed,
            include_models=include_models,
            normalize_method=normalize_method
        )
    except Exception:
        raise Exception('Model requested is not available. Please choose regressor or classifier.')
    return final_model


def model_testing(
    *,
    model_path,
    input_file,
    samples_col,
    label_col,
    model_type
):
    """
    Load a prediction model and use it to predict label values in an independent dataset.

    Parameters
    ----------
    model_path : str
        the path to saved model.
    input_file : str
        the file with the test dataset.
    samples_col : str
        the name of samples column.
    label_col : str
        the name of the target column
    model_type :
        regressor or classifier

    Returns
    -------
    DataFrame
        the results of predictions.

    """
    model = joblib.load(model_path)
    testing_df = pd.read_csv(input_file, sep='\t', index_col=samples_col)
    testing_df = testing_df.dropna(subset=[label_col])
    testing_df.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
    x_set = testing_df.drop(columns=label_col)
    model_func = {'classifier': test_classifier, 'regressor': test_regressor}
    unseen_predictions = model_func.get(model_type)(
        y_col=testing_df[label_col],
        output=input_file.split('.')[0],
        model=model,
        x_set=x_set
    )
    prediction_df = unseen_predictions[['Label']]
    return prediction_df


def normalize_data(
    *,
    method='gene_length',
    genes_info=None,
    genes_col='HGNC symbol',
    length_col='gene_length',
    data_file,
    samples_col,
):
    """
    Normalize dataset using gene_length, minmax, maxabs, zscore or robust

    Parameters
    ----------
    method : str
        the normalization method. [zscore, gene_length, minmax, maxabs, robust]
    genes_info : str
        file containing the genes and their lengths. if gene_length method chosen with no file, info will be retrieved from ensembl database.
    genes_col : str
        the column containing genes (if genes_info file is provided)
    length_col : str
        the columns containing genes length (if genes_info file is provided)
    data_file : str
        file containing dataset to be normalized.
    samples_col : str
        the column containing samples ids.

    Returns
    -------
    DataFrame
        a df with the normalized dataset.

    """
    scores_df = pd.read_csv(data_file, sep=r'\s+')
    if method == 'gene_length':
        scores_df = gene_length_normalize(
            genes_info=genes_info, genes_col=genes_col, length_col=length_col, scores_df=scores_df,
            samples_col=samples_col
        )
    elif method == 'maxabs':
        for col in scores_df.columns:
            scores_df[col] = scores_df[col]/scores_df[col].abs().max()
    elif method == 'minmax':
        for col in scores_df.columns:
            scores_df[col] = (scores_df[col] - scores_df[col].min()) / (scores_df[col].max() - scores_df[col].min())
    elif method == 'zscore ':
        scores_df.std(ddof=0)
        for col in scores_df.columns:
            scores_df[col] = (scores_df[col] - scores_df[col].mean()) / scores_df[col].std()
    elif method == 'robust':
        for col in scores_df.columns:
            scores_df[col] = (scores_df[col] - scores_df[col].median()) / (scores_df[col].quantile(0.75) - scores_df[col].quantile(0.25))
    else:
        raise Exception(
            'This function does not support the normalization method you selected. Methods: [zscore, gene_length, minmax, maxabs, robust]'
        )
    return scores_df