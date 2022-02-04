import multiprocessing
import re
from functools import partial

import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm
from scipy.stats import pearsonr
from tqdm import tqdm


def calc_corr(
    *,
    first_file,
    second_file,
    samples_col,
    output_file,
):
    """
    Calculate the pearson's correlation between same genes in two scoring matices.

    Parameters
    ----------
    first_file : str
        the path to the first scores file.
    second_file : str
        the path to the second scores file.
    samples_col : str
        the column containing the samples IDs.
    output_file : str
        the path to the output file with correlation values.

    Returns
    -------
    pd.DataFrame
        the correlation dataframe.

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
    first_df = pd.read_csv(first_file, sep='\t', index_col=False)
    second_df = pd.read_csv(second_file, sep='\t', index_col=False)
    for gene in tqdm(genes, desc='calculating correlation'):
        gene_df = pd.merge(first_df[[samples_col, gene]], second_df[[samples_col, gene]], on=samples_col)
        gene_df.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
        corr, pval = pearsonr(gene_df[gene + '_x'], gene_df[gene + '_y'])
        corr_info.append([gene, corr, pval])
    corr_df = pd.DataFrame(corr_info, columns=['target', 'corr', 'p_value']).sort_values(by=['p_value'])
    corr_df.to_csv(output_file, sep='\t', index=False)
    return corr_df


def run_linear(gene_col, x_set, y_set):
    """
    Helper function to run linear regression association.

    Parameters
    ----------
    gene_col : tuple
        a tuple from df.iteritems()
    x_set : pd.DataFrame
        a dataframe with the covariates.
    y_set : pd.Series
        the target.

    Returns
    -------
    List
         a list with gene name, pvalues, coefs and std err.

    """
    x_set[gene_col[0]] = gene_col[1]
    x_set = sm.add_constant(x_set)
    linear_model = sm.OLS(y_set, x_set)
    result = linear_model.fit()
    pval = list(result.pvalues)
    beta_coef = list(result.params)[-1]
    std_err = result.bse[-1]
    return [gene_col[0]] + pval + [beta_coef, std_err]


def run_logit(gene_col, x_set, y_set):
    """
    Helper function to run logistic regression association.

    Parameters
    ----------
    gene_col : tuple
        a tuple from df.iteritems()
    x_set : pd.DataFrame
        a dataframe with the covariates.
    y_set : pd.Series
        the target.

    Returns
    -------
    List
         a list with gene name, pvalues, coefs and std err.

    """
    try:
        x_set[gene_col[0]] = gene_col[1]
        x_set = sm.add_constant(x_set)
        logit_model = sm.Logit(y_set, x_set)
        result = logit_model.fit()
        pval = list(result.pvalues)
        std_err = result.bse[-1]
    except:
        return
    return [gene_col[0]] + pval + [std_err]


def run_mannwhitneyu(*, df, genes, cases_column, **kwargs):
    """
    Helper function to run mann whitney u association.

    Parameters
    ----------
    df : df.DataFrame
        the dataframe with all info.
    genes : List
        the list of genes to run the association
    cases_column : str
        the taget column
    kwargs

    Returns
    -------
    pd.DataFrame
        dataframe containing the p_values.

    """
    p_values = []
    df_by_cases = df.groupby(cases_column)
    kwargs['logger'].info(df_by_cases.size())
    if kwargs['cases'] and kwargs['controls']:
        cc = [kwargs['cases'], kwargs['controls']]
    else:
        cc = list(df_by_cases.groups.keys())
    if len(cc) > 2:
        Warning('There are more than two categories here. We will only consider the first two categories.')
    for gene in tqdm(genes, desc='Calculating p_values for genes'):
        case_0 = df_by_cases.get_group(cc[0])[gene].tolist()
        case_1 = df_by_cases.get_group(cc[1])[gene].tolist()
        try:
            u_statistic, p_val = stats.mannwhitneyu(case_0, case_1)
        except:
            continue
        p_values.append([gene, u_statistic, p_val])
    cols = ['genes', 'statistic', 'p_value']
    p_values_df = pd.DataFrame(p_values, columns=cols).sort_values(by=['p_value'])
    return p_values_df


def run_ttest(*, df, genes, cases_column, **kwargs):
    p_values = []
    df_by_cases = df.groupby(cases_column)
    kwargs['logger'].info(df_by_cases.size())
    if kwargs['cases'] and kwargs['controls']:
        cc = [kwargs['cases'], kwargs['controls']]
    else:
        cc = list(df_by_cases.groups.keys())
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
    cols = ['genes', 'statistic', 'p_value']
    p_values_df = pd.DataFrame(p_values, columns=cols).sort_values(by=['p_value'])
    return p_values_df


def get_pvals_logit(*, df, genes, cases_column, **kwargs):
    covariates = kwargs['covariates']
    df[cases_column] = np.interp(
        df[cases_column], (df[cases_column].min(), df[cases_column].max()), (0, 1))
    kwargs['logger'].info(df.groupby(cases_column).size())
    y_set = df[[cases_column]]
    x_set = df[covariates]
    genes_df = df[genes]
    pool = multiprocessing.Pool(processes=kwargs['processes'])
    partial_func = partial(run_logit, x_set=x_set, y_set=y_set)
    p_values = list(pool.imap(partial_func, genes_df.iteritems()))
    cols = ['genes', 'const_pval'] + covariates + ['p_value', 'std_err']
    p_values_df = pd.DataFrame(list(filter(None, p_values)), columns=cols).sort_values(by=['p_value'])
    return p_values_df


def get_pvals_linear(*, df, genes, cases_column, **kwargs):
    covariates = kwargs.get('covariates')
    y_set = df[[cases_column]]
    x_set = df[covariates]
    genes_df = df[genes]
    pool = multiprocessing.Pool(processes=kwargs['processes'])
    partial_func = partial(run_linear, x_set=x_set, y_set=y_set)
    p_values = list(pool.imap(partial_func, genes_df.iteritems()))
    cols = ['genes', 'const_pval'] + covariates + ['p_value', 'beta_coef', 'std_err']
    p_values_df = pd.DataFrame(p_values, columns=cols).sort_values(by=['p_value'])
    return p_values_df
