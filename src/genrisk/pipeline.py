# -*- coding: utf-8 -*-
import os
import subprocess

import joblib
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from .prediction_models import regression_model, classification_model, test_classifier, test_regressor
from .association_analysis import run_mannwhitneyu, run_ttest, get_pvals_logit, get_pvals_linear
from .gene_scoring import get_gene_info, plink_process, combine_scores

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
    adj_pval,
    covariates=None,
    cases=None,
    controls=None,
    processes=1,
    logger
):
    """
    Calculate the significance of a gene in a population using Mann-Whitney-U test.
    \f

    :param test: the type of statistical test to use, choices are: t-test, mannwhitenyu, GLM, logit.
    :param scores_file: dataframe containing the scores of genes across samples.
    :param info_file: a file containing the information of the sample.
    :param output_file: a path to save the output file.
    :param genes: a list of the genes to calculate the significance. if None will calculate for all genes.
    :param cases_column: the name of the column containing cases and controls information.
    :param samples_column: the name of the column contining samples IDs.
    :param adj_pval: the method for pvalue adjustment.
    :param controls: the name of the controls category.
    :param cases:  the name of the cases category.
    :param covariates: the covariates of the phenotype.
    :param processes: number of processes working in parallel.

    :return: dataframe with genes and their p_values
    """
    logger.info("Reading scores file...")
    scores_df = pd.read_csv(scores_file, sep=r'\s+', index_col=samples_column)
    scores_df.replace([np.inf, -np.inf], 0, inplace=True)
    scores_df.fillna(0, inplace=True)
    scores_df = scores_df.loc[:, scores_df.var() != 0.0].reset_index()
    logger.info("Reading info file...")
    genotype_df = pd.read_csv(info_file, sep='\t')
    genotype_df.dropna(subset=[cases_column], inplace=True)
    logger.info("Processing files...")
    merged_df = genotype_df.merge(scores_df, how='inner', on=samples_column)
    merged_df.replace([np.inf, -np.inf], 0, inplace=True)
    merged_df.fillna(0, inplace=True)
    if genes is None:
        genes = scores_df.columns.tolist()[1:]
    del scores_df
    if covariates:
        covariates = covariates.split(',')
    args = {
        'processes': processes, 'cases': cases, 'controls': controls, 'covariates': covariates,
    }
    logger.info("Calculating p_values using the following test: " + test)
    try:
        p_values_df = association_functions.get(test)(df=merged_df, genes=genes, cases_column=cases_column, **args)
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
    logger,
):
    """
    Calculate association significance between two groups using betareg.
    \f

    :param scores_file: the path to the scores file.
    :param pheno_file:  the path to the phenotypes and covariates file.
    :param samples_col: the name of the column containing the samples IDs.
    :param cases_col: the name of the column containing the case/controls.
    :param output_path: the path to the output file.
    :param covariates: the covariates used in calculations, written with no space and comma in between (e.g PC1,PC2)
    :param processes: number of processes to parallelize.

    :return:
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
             "--processes", str(processes)]
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
):
    """
    Create a prediction model (classifier or regressor) using the provided dataset.
    \f

    :param model_name: the name of the prediction model.
    :param model_type: type of model (reg or classifier).
    :param y_col: the column containing the target (qualitative or quantitative).
    :param imbalanced: True means data is imbalanced.
    :param normalize: True if data needs normalization.
    :param folds: how many folds for cross-validation.
    :param training_set: the training set for the model.
    :param testing_set: if exists an extra evaluation step will be done using the testing set.
    :param test_size: the size to split the training/testing set.
    :param metric: the metric to evaluate the best model.
    :return: the metrics.
    """
    if model_type == 'regressor':
        if not metric:
            metric = 'RMSE'
        final_model = regression_model(
            y_col=y_col,
            training_set=training_set,
            normalize=normalize,
            test_size=test_size,
            folds=folds,
            metric=metric,
            model_name=model_name,
            testing_set=testing_set
        )
    elif model_type == 'classifier':
        if not metric:
            metric = 'AUC'
        final_model = classification_model(
            y_col=y_col,
            training_set=training_set,
            normalize=normalize,
            test_size=test_size,
            folds=folds,
            metric=metric,
            model_name=model_name,
            testing_set=testing_set,
            imbalanced=imbalanced
        )
    else:
        return Exception('Model requested is not available. Please choose regressor or classifier.')
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
    Load a prediction model and use it to predict label values in a dataset.
    \f

    :param model_path: the path to saved model.
    :param input_file: the file with features and label.
    :param samples_col: the name of samples column.
    :param label_col: the name of the target column.
    :param model_type: regressor or classifier

    :return: a dataframe with prediction values.
    """
    model = joblib.load(model_path)
    testing_df = pd.read_csv(input_file, sep='\t', index_col=samples_col)
    x_set = testing_df.drop(columns=label_col)
    if model_type == 'classifier':
        unseen_predictions = test_classifier(
            y_col=testing_df[label_col],
            output=input_file.split('.')[0],
            model=model,
            x_set=x_set,
        )
    else:
        unseen_predictions = test_regressor(
            y_col=testing_df[label_col],
            output=input_file.split('.')[0],
            model=model,
            x_set=x_set,
        )
    prediction_df = unseen_predictions[['Label']]
    return prediction_df
