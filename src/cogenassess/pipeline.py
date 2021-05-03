# -*- coding: utf-8 -*-
import os
import subprocess

import pandas as pd
import pycaret.classification as cl
import pycaret.regression as pyreg
from pycaret.utils import check_metric
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

PATH = os.path.abspath(os.path.join((__file__), os.pardir, os.pardir, os.pardir))
BETAREG_SHELL = os.path.join(PATH, 'scripts', 'betareg_shell.R')
PLOT_SHELL = os.path.join(PATH, 'scripts', 'plot_script.R')


def find_pvalue(
    *,
    scores_df,
    genotype_file,
    output_file,
    genes=None,
    cases_column,
    samples_column,
    pc_file=None,
    test='mannwhitneyu',
    adj_pval,
):
    """
    Calculate the significance of a gene in a population using Mann-Whitney-U test.
    :param pc_file: if there is an extra file with PC values.
    :param test: the type of statistical test to use, choices are: t-test, mannwhitenyu, GLM, logit.
    :param scores_df: dataframe containing the scores of genes across samples.
    :param genotype_file: a file containing the information of the sample.
    :param output_file: a path to save the output file.
    :param genes: a list of the genes to calculate the significance. if None will calculate for all genes.
    :param cases_column: the name of the column containing cases and controls information.
    :param samples_column: the name of the column contining samples IDs.
    :param adj_pval: the method for pvalue adjustment.
    :return: dataframe with genes and their p_values
    """
    genotype_df = pd.read_csv(genotype_file, sep=r'\s+', usecols=[samples_column, cases_column])
    merged_df = pd.merge(genotype_df, scores_df, on=samples_column)
    df_by_cases = merged_df.groupby(cases_column)
    cases = list(df_by_cases.groups.keys())
    p_values = []
    if genes is None:
        genes = scores_df.columns.tolist()[1:]
    if test == 'mannwhitneyu':
        for gene in tqdm(genes, desc='Calculating p_values for genes'):
            case_0 = df_by_cases.get_group(cases[0])[gene].tolist()
            case_1 = df_by_cases.get_group(cases[1])[gene].tolist()
            try:
                u_statistic, p_val = stats.mannwhitneyu(case_0, case_1, alternative='greater')
            except:
                continue
            p_values.append([gene, u_statistic, p_val])
        p_values_df = pd.DataFrame(p_values, columns=['genes', 'statistic', 'p_value']).sort_values(by=['p_value'])
    elif test == 'ttest_ind':
        for gene in tqdm(genes, desc='Calculating p_values for genes'):
            case_0 = df_by_cases.get_group(cases[0])[gene].tolist()
            case_1 = df_by_cases.get_group(cases[1])[gene].tolist()
            try:
                statistic, p_val = stats.ttest_ind(case_0, case_1)
            except:
                continue
            p_values.append([gene, statistic, p_val])
        p_values_df = pd.DataFrame(p_values, columns=['genes', 'statistic', 'p_value']).sort_values(by=['p_value'])
    elif test == 'logit':
        if not pc_file:
            raise Exception("Need principle components file.")
        pc_df = pd.read_csv(pc_file, sep=r'\s+', index_col=False)
        merged_df = pd.merge(merged_df, pc_df, on=samples_column)
        for gene in tqdm(genes, desc='Calculating p_values for genes'):
            X = merged_df[[gene, 'PC1', 'PC2', 'PC3']]
            X = sm.add_constant(X)
            Y = merged_df[[cases_column]]
            try:
                logit_model = sm.Logit(Y, X)
                result = logit_model.fit()
            except:
                continue
            pval = list(result.pvalues)
            # add std err
            p_values.append([gene] + pval)
        p_values_df = pd.DataFrame(
            p_values, columns=['genes', 'const_pval', 'p_value', 'PC1_pval', 'PC2_pvcal', 'PC3_pval']
        ).sort_values(by=['p_value'])
    elif test == 'glm':
        if not pc_file:
            raise Exception("Need principle components file.")
        pc_df = pd.read_csv(pc_file, sep=r'\s+', index_col=False)
        merged_df = pd.merge(merged_df, pc_df, on=samples_column)
        for gene in tqdm(genes, desc='Calculating p_values for genes'):
            X = merged_df[[gene, 'PC1', 'PC2', 'PC3']]
            X = sm.add_constant(X)
            Y = merged_df[[cases_column]]
            try:
                glm_model = sm.GLM(Y, X)
                result = glm_model.fit()
            except:
                continue
            pval = list(result.pvalues)
            beta_coef = list(result.params)[1]
            p_values.append([gene] + pval + [beta_coef])
        p_values_df = pd.DataFrame(
            p_values, columns=['genes', 'const_pval', 'p_value', 'PC1_pval', 'PC2_pvcal', 'PC3_pval', 'beta_coef']
        ).sort_values(by=['p_value'])
    else:
        raise Exception("The test you selected is not valid.")
    if adj_pval:
        adjusted = multipletests(list(p_values_df['p_value']), method=adj_pval)
        p_values_df[adj_pval + '_adj_pval'] = list(adjusted)[1]
    p_values_df.to_csv(output_file, sep='\t', index=False)
    return p_values_df


def betareg_pvalues(
    *,
    scores_file,
    pheno_file,
    samples_col,
    cases_col,
    nprocesses=3,
    output_path,
    covariates
):
    """
    Calculate association significance between two groups using betareg.
    :param scores_file: the path to the scores file.
    :param pheno_file:  the path to the phenotypes and covariates file.
    :param samples_col: the name of the column containing the samples IDs.
    :param cases_col: the name of the column containing the case/controls.
    :param nprocesses: the number of processes used in parallel.
    :param output_path: the path to the output file.
    :param covariates: the covariates used in calculations, written with no space and comma in between (e.g PC1,PC2)
    :return:
    """
    p = subprocess.call(
        ["Rscript", BETAREG_SHELL,
         "-s", scores_file,
         "--phenofile", pheno_file,
         "--samplescol", samples_col,
         "--casescol", cases_col,
         "--nprocesses", str(nprocesses),
         "-o", output_path,
         "--covariates", covariates], shell=True
    )


def r_visualize(
    *,
    genescol_1,
    genescol_2,
    info_file,
    pvals_file,
    qq_output,
    manhattan_output,
    pvalcol
):
    """
    Visualize the results of association test. Manhattan plot and QQ-plot.
    :param genescol_1: the name of the column containing the genes in the pvals_file.
    :param genescol_2: the name of the column containing the genes in the info_file.
    :param info_file: the file containing gene information.
    :param pvals_file: file containing the pvalues.
    :param qq_output: the name for the qq plot file.
    :param manhattan_output: the name for manhattan plot file.
    :param pvalcol: the column containing the pvalues.
    :return:
    """
    subprocess.run(
        ["Rscript", PLOT_SHELL,
         "-p", pvals_file,
         "-i", info_file,
         "--qq_output", qq_output,
         "--manhattan_output", manhattan_output,
         "--genescol_1", genescol_1,
         "--genescol_2", genescol_2,
         "--pvalcol", pvalcol]
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
    testing_set=None,
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
        reg = pyreg.setup(target=y_col, data=training_set, normalize=normalize, train_size=1-test_size, fold=folds,
                          silent=True)
        best_model = pyreg.compare_models(sort=metric)
        reg_model = pyreg.create_model(best_model)
        reg_tuned_model = pyreg.tune_model(reg_model, optimize=metric)
        pyreg.plot_model(reg_tuned_model, save=True)
        pyreg.plot_model(reg_tuned_model, plot='feature', save=True)
        pyreg.plot_model(reg_tuned_model, plot='error', save=True)
        final_model = pyreg.finalize_model(reg_tuned_model)
        if len(testing_set.index) != 0:
            unseen_predictions = pyreg.predict_model(final_model, data=testing_set)
            r2 = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'R2')
            rmse = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'RMSE')
            metrics = ['R2: ' + str(r2), 'RMSE: ' + str(rmse)]
        pyreg.save_model(final_model, model_name)
        with open('model.log', 'w') as f:
            f.writelines("%s\n" % output for output in reg)
            f.writelines("%s\n" % output for output in metrics)
    elif model_type == 'classifier':
        if not metric:
            metric = 'AUC'
        classifier = cl.setup(target=y_col, fix_imbalance=imbalanced, data=training_set, train_size=1-test_size,
                              silent=True, fold=folds)
        best_model = cl.compare_models(sort=metric)
        cl_model = cl.create_model(best_model)
        cl_tuned_model = pyreg.tune_model(cl_model, optimize=metric)
        cl.plot_model(cl_tuned_model, plot='pr', save=True)
        cl.plot_model(cl_tuned_model, plot='confusion_matrix', save=True)
        cl.plot_model(cl_tuned_model, plot='feature', save=True)
        final_model = cl.finalize_model(cl_tuned_model)
        if len(testing_set.index) != 0:
            unseen_predictions = cl.predict_model(final_model, data=testing_set)
            auc = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'AUC')
            accuracy = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'Accuracy')
            metrics = ['AUC: ' + str(auc), 'Accuracy: ' + str(accuracy)]
        cl.save_model(final_model, model_name)
        with open('model.log', 'w') as f:
            f.writelines("%s\n" % output for output in classifier)
            f.writelines("%s\n" % output for output in metrics)
    else:
        return Exception('Model requested is not available. Please choose regressor or classifier.')
    return final_model
