# -*- coding: utf-8 -*-
import os
import subprocess

import pandas as pd
from pybiomart import Dataset
import scipy.stats as stats
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

PATH = os.path.abspath(os.path.join((__file__), os.pardir, os.pardir, os.pardir))
BETAREG_SHELL = os.path.join(PATH, 'scripts', 'betareg_shell.R')
PLOT_SHELL = os.path.join(PATH, 'scripts', 'plot_script.R')


def normalize_gene_len(
    *,
    genes_lengths_file=None,
    matrix_file,
    samples_col,
    output_path
):
    """
    Normalize matrix by gene length.

    :param genes_lengths_file: a file containing genes, and their start and end bps.
    :param matrix_file: a tsv file containing a matrix of samples and their scores across genes.
    :param output_path: the path to save the normalized matrix.
    :return: a normalized dataframe.
    """
    if genes_lengths_file:
        genes_df = pd.read_csv(genes_lengths_file, sep=r'\s+')
    else:
        gene_dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
        genes_df = gene_dataset.query(
            attributes=['external_gene_name', 'start_position', 'end_position'],
            only_unique=False,
        )
    genes_lengths = {
        row['Gene name']: round((row['Gene end (bp)'] - row['Gene start (bp)']) / 1000, 3)
        for _, row in genes_df.iterrows()
    }
    scores_df = pd.read_csv(matrix_file, sep=r'\s+')
    unnormalized = []
    for (name, data) in tqdm(scores_df.iteritems(), desc="Normalizing genes scores"):
        if name == samples_col:
            continue
        if name not in genes_lengths.keys():
            unnormalized.append(name)
            continue
        # normalize genes by length
        scores_df[name] = round(scores_df[name] / genes_lengths[name], 5)
    # drop genes with unknown length
    scores_df = scores_df.drop(unnormalized, axis=1)
    if output_path:
        scores_df.to_csv(output_path, sep='\t', index=False)
    return scores_df


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
    :param pc_file:
    :param test:
    :param scores_df: dataframe containing the scores of genes across samples.
    :param genotype_file: a file containing the information of the sample.
    :param output_file: a path to save the output file.
    :param genes: a list of the genes to calculate the significance. if None will calculate for all genes.
    :param cases_column: the name of the column containing cases and controls information.
    :return: dataframe with genes and their p_values
    """
    genotype_df = pd.read_csv(genotype_file, sep=r'\s+', usecols=[samples_column, cases_column])
    merged_df = pd.merge(genotype_df, scores_df, on=samples_column)
    df_by_cases = merged_df.groupby(cases_column)
    cases = list(df_by_cases.groups.keys())
    p_values = []
    if genes is None:
        genes = scores_df.columns.tolist()[1:]
    # improvement: make it so that more than one statistical test can be used
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
    pc_file,
    samples_col,
    cases_col,
    nprocesses=3,
    output_path,
    covariates
):
    subprocess.call(
        ["Rscript", BETAREG_SHELL,
         "-s", scores_file,
         "--phenofile", pheno_file,
         "--pcfile", pc_file,
         "--samplescol", samples_col,
         "--casescol", cases_col,
         "--nprocesses", str(nprocesses),
         "-o", output_path,
         "--covariates", covariates]
    )


def r_visualize(
    *,
    genescol_1,
    genescol_2,
    info_file,
    pvals_file,
    qq_output,
    manhattan_output
):
    subprocess.run(
        ["Rscript", PLOT_SHELL,
         "-p", pvals_file,
         "-i", info_file,
         "--qq_output", qq_output,
         "--manhattan_output", manhattan_output,
         "--genescol-1", genescol_1,
         "--genescol-2", genescol_2]
    )


def merge_matrices(
    *,
    directory,
    output_path,
    samples_col,
    scores_col,
    file_suffix='.tsv'
):
    """
    Merges multiple files in a directory, each file should contain the score of a gene across the samples.

    :param directory: the directory that contains files to merge.
    :param output_path: the path for the merged tsv file.
    :return: a dataframe combining all information from files.
    """
    full_data = pd.DataFrame(data=None, columns=samples_col)
    for filename in tqdm(os.listdir(directory), desc="merging matrices"):
        if not filename.endswith(file_suffix):
            continue
        data = pd.read_csv(os.path.join(directory, filename), sep=r'\s+', usecols=samples_col + [scores_col])
        gene_name = filename.split('.')[0]
        data = data.rename(columns={scores_col: gene_name})
        full_data = pd.merge(data, full_data, on=samples_col, how='left')
    full_data.to_csv(output_path, sep='\t', index=False)
    return full_data

def merge_files_fun(
    *,
    input_dir,
    samples_col,
):
    folder = os.listdir(input_dir)
    df = pd.read_csv(os.path.join(input_dir, folder[0]), sep='\t', index_col=False).sort_values(by=[samples_col])
    for filename in tqdm(folder[1:], desc="merging files"):
        new_df = pd.read_csv(os.path.join(input_dir, filename), sep='\t', index_col=False).sort_values(by=[samples_col])
        new_df.drop(columns=[samples_col])
        df = pd.concat([df, new_df])
    return df
