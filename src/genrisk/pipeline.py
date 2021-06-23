# -*- coding: utf-8 -*-
import os
import subprocess
import random
import multiprocessing

import joblib
import matplotlib.pyplot as plt
from functools import partial
import numpy as np

import pandas as pd
import pycaret.classification as cl
import pycaret.regression as pyreg
import scipy.stats as stats
import statsmodels.api as sm
from pycaret.utils import check_metric
import sklearn.metrics as metrics
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

PATH = os.path.abspath(os.path.join((__file__), os.pardir, os.pardir, os.pardir))
BETAREG_SHELL = os.path.join(PATH, 'scripts', 'betareg_shell.R')


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
):
    """
    Calculate the significance of a gene in a population using Mann-Whitney-U test.

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

    :return: dataframe with genes and their p_values
    """
    scores_df = pd.read_csv(scores_file, sep='\t')
    scores_df.replace([np.inf, -np.inf], 0, inplace=True)
    scores_df.fillna(0, inplace=True)
    scores_df = scores_df.loc[:, scores_df.var() != 0.0]
    genotype_df = pd.read_csv(info_file, sep='\t')
    genotype_df.dropna(subset=[cases_column], inplace=True)
    merged_df = genotype_df.merge(scores_df, how='inner', on=samples_column)
    merged_df.replace([np.inf, -np.inf], 0, inplace=True)
    merged_df.fillna(0, inplace=True)
    if genes is None:
        genes = scores_df.columns.tolist()[1:]
    del scores_df
    df_by_cases = merged_df.groupby(cases_column)
    if covariates:
        covariates = covariates.split(',')
    if cases and controls:
        cc = [cases, controls]
    else:
        cc = list(df_by_cases.groups.keys())
    p_values = []
    if test == 'mannwhitneyu':
        if len(cc) > 2:
            Warning('There are more than two categories here. We will only consider the first two categories.')
        for gene in tqdm(genes, desc='Calculating p_values for genes'):
            case_0 = df_by_cases.get_group(cc[0])[gene].tolist()
            case_1 = df_by_cases.get_group(cc[1])[gene].tolist()
            try:
                u_statistic, p_val = stats.mannwhitneyu(case_0, case_1, alternative='greater')
            except:
                continue
            p_values.append([gene, u_statistic, p_val])
        p_values_df = pd.DataFrame(p_values, columns=['genes', 'statistic', 'p_value']).sort_values(by=['p_value'])
    elif test == 'ttest_ind':
        if len(cc) > 2:
            Warning('There are more than two categories here. We will only consider the first two categories.')
        for gene in tqdm(genes, desc='Calculating p_values for genes'):
            case_0 = df_by_cases.get_group(cc[0])[gene].tolist()
            case_1 = df_by_cases.get_group(cc[1])[gene].tolist()
            try:
                statistic, p_val = stats.ttest_ind(case_0, case_1)
            except:
                continue
            p_values.append([gene, statistic, p_val])
        p_values_df = pd.DataFrame(p_values, columns=['genes', 'statistic', 'p_value']).sort_values(by=['p_value'])
    elif test == 'logit':
        Y = merged_df[[cases_column]]
        X = merged_df[covariates]
        genes_df = merged_df[genes]
        del merged_df
        pool = multiprocessing.Pool(processes=processes)
        partial_func = partial(run_logit, X=X, Y=Y)
        p_values = list(pool.imap(partial_func, genes_df.iteritems()))
        cols = ['genes', 'const_pval'] + covariates + ['p_value' + 'std_err']
        p_values_df = pd.DataFrame(p_values, columns=cols).sort_values(by=['p_value'])
    elif test == 'linear':
        Y = merged_df[[cases_column]]
        X = merged_df[covariates]
        genes_df = merged_df[genes]
        del merged_df
        pool = multiprocessing.Pool(processes=processes)
        partial_func = partial(run_linear, X=X, Y=Y)
        p_values = list(pool.imap(partial_func, genes_df.iteritems()))
        cols = ['genes', 'const_pval'] + covariates + ['p_value', 'beta_coef', 'std_err']
        p_values_df = pd.DataFrame(p_values, columns=cols).sort_values(by=['p_value'])
    elif test == 'glm':
        Y = merged_df[[cases_column]]
        X = merged_df[covariates]
        genes_df = merged_df[genes]
        del merged_df
        pool = multiprocessing.Pool(processes=processes)
        partial_func = partial(run_glm, X=X, Y=Y)
        p_values = list(pool.imap(partial_func, genes_df.iteritems()))
        cols = ['genes', 'const_pval'] + covariates + ['p_value', 'beta_coef', 'std_err']
        p_values_df = pd.DataFrame(p_values, columns=cols).sort_values(by=['p_value'])
    else:
        raise Exception("The test you selected is not available.")
    if adj_pval:
        adjusted = multipletests(list(p_values_df['p_value']), method=adj_pval)
        p_values_df[adj_pval + '_adj_pval'] = list(adjusted)[1]
    p_values_df.to_csv(output_file, sep='\t', index=False)
    return p_values_df


def run_linear(gene_col, X, Y):
    """
    Helper function to run linear regression association.

    :param gene_col: a tuple from df.iteritems()
    :param X: the covariates.
    :param Y: the target.

    :return: a list with gene name, pvalues, coefs and std err.
    """
    X[gene_col[0]] = gene_col[1]
    X = sm.add_constant(X)
    linear_model = sm.OLS(Y, X)
    result = linear_model.fit()
    pval = list(result.pvalues)
    beta_coef = list(result.params)[-1]
    std_err = result.bse[-1]
    return [gene_col[0]] + pval + [beta_coef, std_err]


def run_glm(gene_col, X, Y):
    """
    Helper function to run Generalize linear model regression association.

    :param gene_col: a tuple from df.iteritems()
    :param X: the covariates.
    :param Y: the target.

    :return: a list with gene name, pvalues, coefs and std err.
    """
    X[gene_col.name] = gene_col
    X = sm.add_constant(X)
    glm_model = sm.GLM(Y, X)
    result = glm_model.fit()
    pval = list(result.pvalues)
    beta_coef = list(result.params)[-1]
    std_err = result.bse[-1]
    return [gene_col.name] + pval + [beta_coef, std_err]


def run_logit(gene_col, X, Y):
    """
    Helper function to run logistic regression association.

    :param gene_col: a tuple from df.iteritems()
    :param X: the covariates.
    :param Y: the target.

    :return: a list with gene name, pvalues, and std err.
    """
    X[gene_col.name] = gene_col
    X = sm.add_constant(X)
    logit_model = sm.Logit(Y, X)
    result = logit_model.fit()
    pval = list(result.pvalues)
    std_err = result.bse[-1]
    return [gene_col.name] + pval + [std_err]


def betareg_pvalues(
    *,
    scores_file,
    pheno_file,
    samples_col,
    cases_col,
    output_path,
    covariates,
    processes
):
    """
    Calculate association significance between two groups using betareg.

    :param scores_file: the path to the scores file.
    :param pheno_file:  the path to the phenotypes and covariates file.
    :param samples_col: the name of the column containing the samples IDs.
    :param cases_col: the name of the column containing the case/controls.
    :param output_path: the path to the output file.
    :param covariates: the covariates used in calculations, written with no space and comma in between (e.g PC1,PC2)
    :param processes: number of processes to parallelize.

    :return:
    """
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
    metrics = None
    if model_type == 'regressor':
        if not metric:
            metric = 'RMSE'
        setup = pyreg.setup(target=y_col, data=training_set, normalize=normalize, train_size=1 - test_size, fold=folds,
                            silent=True, session_id=random.randint(1, 2147483647))
        best_model = pyreg.compare_models(sort=metric)
        pyreg.pull().to_csv(model_name + '_compare_models.tsv', sep='\t', index=False)
        reg_model = pyreg.create_model(best_model)
        reg_tuned_model = pyreg.tune_model(reg_model, optimize=metric)
        pyreg.pull().to_csv(model_name + '_tuned_model.tsv', sep='\t', index=False)
        final_model = pyreg.finalize_model(reg_tuned_model)
        pyreg.plot_model(final_model, save=True)
        pyreg.plot_model(final_model, plot='feature', save=True)
        pyreg.plot_model(final_model, plot='error', save=True)
        if len(testing_set.index) != 0:
            unseen_predictions = pyreg.predict_model(final_model, data=testing_set)
            r2 = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'R2')
            rmse = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'RMSE')
            metrics = ['R2: ' + str(r2), 'RMSE: ' + str(rmse)]
            unseen_predictions.to_csv(model_name + '_external_testing_results.tsv', sep='\t', index=False)
        pyreg.save_model(final_model, model_name)
    elif model_type == 'classifier':
        if not metric:
            metric = 'AUC'
        setup = cl.setup(target=y_col, fix_imbalance=imbalanced, data=training_set, train_size=1 - test_size,
                         silent=True, fold=folds, session_id=random.randint(1, 2147483647))
        best_model = cl.compare_models(sort=metric)
        cl_model = cl.create_model(best_model)
        cl_tuned_model = cl.tune_model(cl_model, optimize=metric)
        final_model = cl.finalize_model(cl_tuned_model)
        cl.plot_model(final_model, plot='pr', save=True)
        cl.plot_model(final_model, plot='confusion_matrix', save=True)
        cl.plot_model(final_model, plot='feature', save=True)
        if len(testing_set.index) != 0:
            unseen_predictions = cl.predict_model(final_model, data=testing_set)
            auc = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'AUC')
            accuracy = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'Accuracy')
            metrics = ['AUC: ' + str(auc), 'Accuracy: ' + str(accuracy)]
            unseen_predictions.to_csv(model_name + '_external_testing_results.tsv', sep='\t', index=False)
        cl.save_model(final_model, model_name)
        cl.pull().to_csv(model_name + '_evaluation.tsv', sep='\t', index=False)
        setup_list = list(setup)
        setup_list[:] = [x for x in setup_list if type(x) == list]
        setup[3][0][1].to_csv(model_name + '_setup.tsv', sep='\t')
        setup[2][1].to_csv(model_name + '_compare_models.tsv', sep='\t')
        setup[0][-1].to_csv(model_name + '_tuned_model.tsv', sep='\t')
    else:
        return Exception('Model requested is not available. Please choose regressor or classifier.')
    return metrics, final_model


def model_testing(
    *,
    model_path,
    input_file,
    samples_col,
    label_col,
    model_type
):
    model = joblib.load(model_path)
    testing_df = pd.read_csv(input_file, sep='\t', index_col=samples_col)
    x_set = testing_df.drop(columns=label_col)
    if model_type == 'classifier':
        unseen_predictions = cl.predict_model(model, data=testing_df)
        report = metrics.classification_report(unseen_predictions[label_col], unseen_predictions.Label)
        acc = metrics.accuracy_score(unseen_predictions[label_col], unseen_predictions.Label)
        auc = metrics.auc(unseen_predictions[label_col], unseen_predictions.Label)
        confusion = metrics.plot_confusion_matrix(x_set, unseen_predictions[label_col])
        confusion.ax_.set_title('Classifier confusion matrix')
        plt.show()
        plt.savefig(input_file.split('.')[0] + '_classifier_confusion_matrix.png')
        textfile = open(input_file.split('.')[0] + "_report.txt", "w")
        textfile.write('Testing model report: \n')
        textfile.write(report + '\n')
        textfile.write('AUC = ' + str(auc) + '\n')
        textfile.write('Accuracy = ' + str(acc) + '\n')
        textfile.close()
    else:
        unseen_predictions = pyreg.predict_model(model, data=testing_df)
        r2 = metrics.r2_score(unseen_predictions[label_col], unseen_predictions.Label)
        rmse = metrics.mean_squared_error(unseen_predictions[label_col], unseen_predictions.Label, squared=False)
        plt.scatter(unseen_predictions.Label, unseen_predictions[label_col], alpha=0.5)
        m, b = np.polyfit(unseen_predictions.Label, unseen_predictions[label_col], 1)
        plt.plot(unseen_predictions.Label, m * unseen_predictions.Label + b, 'r')
        plt.title('Actual vs predicted scatterplot')
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.savefig(input_file.split('.')[0] + '_regressor_scatterplot.png')
        textfile = open(input_file.split('.')[0] + "_report.txt", "w")
        textfile.write('Testing model report: \n')
        textfile.write('R^2 = ' + str(r2) + '\n')
        textfile.write('RMSE = ' + str(rmse) + '\n')
        textfile.close()
    return unseen_predictions

