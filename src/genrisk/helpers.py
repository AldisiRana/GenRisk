# -*- coding: utf-8 -*-
import re

import numpy as np
import pandas as pd
import statsmodels.api as sm


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


def uni_profiles(df, f):
    """
    Merge two dataframes.

    :param df: the main dataframe with all the scores.
    :param f: the file containing the scores of one gene.

    :return: the merged dataframe.
    """
    df2 = pd.read_csv(str(f), usecols=['IID', 'SCORESUM'], sep=r'\s+').astype({'SCORESUM': np.float32})
    r = re.compile("([a-zA-Z0-9_.-]*).profile$")
    gene2 = r.findall(str(f))
    df2.rename(columns={'SCORESUM': gene2[0]}, inplace=True)
    df = pd.merge(df, df2, on='IID')
    return df

