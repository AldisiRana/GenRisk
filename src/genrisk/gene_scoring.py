# -*- coding: utf-8 -*-
import gzip
import os
import re
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import beta
from tqdm import tqdm

from .helpers import uni_profiles


def get_gene_info(
    *,
    annotated_vcf,
    variant_col,
    af_col,
    alt_col='Alt',
    del_col,
    output_dir,
    genes_col,
    maf_threshold=0.01,
    beta_param,
    weight_func='beta'
):
    """
    Create temporary files with variant information for each gene, plus the weights calculated.

    Parameters
    ----------
    annotated_vcf : str
        a file containing the variant, AF, ALT, Gene, and deleterious score.
    variant_col : str
        the name of the variant column.
    af_col : str
        the name of the Allele Frequency column.
    alt_col : str
        the name of the alternate allele column.
    del_col : str
        the name of functional annotation column.
    output_dir : str
        directory to save in temporary files.
    genes_col : str
        the name of genes column.
    maf_threshold : float
        between [0.0-1.0]. the minor allele frequency threshold, default is 0.01.
    beta_param : tuple
        the parameters of the beta function, if chosen for weighting.
    weight_func : str
        the weighting function, beta or log10.

    Returns
    -------
        output directory with all the temporary files.

    """
    skip = 0
    if annotated_vcf.endswith('.gz'):
        with gzip.open(annotated_vcf, 'r') as fin:
            for line in fin:
                if line.decode('utf-8').startswith('##'):
                    skip += 1
    else:
        with open(annotated_vcf, 'r') as file:
            for line in file:
                if line.startswith('##'):
                    skip += 1
    df = pd.read_csv(annotated_vcf, usecols=[variant_col, alt_col, 'INFO'], skiprows=skip, sep=r'\s+', index_col=False)
    info = df['INFO'].str.split(pat=';', expand=True)
    missing_info = info[info.isnull().any(axis=1)].index
    df.drop(missing_info, inplace=True)
    df.reset_index(drop=True, inplace=True)
    info.drop(missing_info, inplace=True)
    info.reset_index(drop=True, inplace=True)
    for col in info.columns:
        val = info[col][0].split('=')
        if len(val) == 1:
            continue
        info.rename(columns={col: val[0]}, inplace=True)
        info[val[0]] = info[val[0]].str.replace(val[0] + '=', r'')
    df = pd.concat([df, info], axis=1)
    df = df[df[af_col].values.astype(float) < maf_threshold]
    df.replace('.', 0.0, inplace=True)
    if weight_func == 'beta':
        df[weight_func] = beta.pdf(df[af_col].values.astype(float), beta_param[0], beta_param[1])
    elif weight_func == 'log10':
        df[weight_func] = -np.log10(df[af_col].values.astype(float))
        df[weight_func].replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
    df['score'] = df[weight_func].values.astype(float) * df[del_col].values.astype(float)
    genes = list(set(df[genes_col]))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    gene_file = output_dir + '.genes'
    with open(gene_file, 'w') as f:
        f.writelines("%s\n" % gene for gene in genes)
    [df[df[genes_col] == gene][[variant_col, alt_col, 'score', genes_col]].to_csv(os.path.join(output_dir, (
        str(gene) + '.w')), index=False, sep='\t') for gene in tqdm(genes, desc="writing w gene files")]
    [df[df[genes_col] == gene][[variant_col, alt_col]].to_csv(os.path.join(output_dir, (str(gene) + '.v')),
                                                              index=False, sep='\t') for gene in
     tqdm(genes, desc="writing v gene files")]
    return output_dir


def combine_scores(
    *,
    input_path,
    output_path,
):
    """
    Combine the files that contain the scores into one file.

    Parameters
    ----------
    input_path : str
        the directory containing scores files.
    output_path : str
        the name of the output file.

    Returns
    -------
    pd.DataFrame
        dataframe with all the scores.

    """
    all_files = [os.path.join(path, name) for path, subdirs, files in os.walk(input_path) for name in files]
    profile_files = [f for f in all_files if re.match(r'.+profile$', f)]
    df = pd.read_csv(str(profile_files[0]), usecols=['IID', 'SCORESUM'], sep=r'\s+').astype({'SCORESUM': np.float32})
    r = re.compile("([a-zA-Z0-9_.-]*).profile$")
    gene = r.findall(str(profile_files[0]))
    df.rename(columns={'SCORESUM': gene[0]}, inplace=True)
    pf = profile_files
    for i in tqdm(range(1, len(pf) - 1), desc='merging in process'):
        df = uni_profiles(df, pf[i])
    df.to_csv(output_path, sep='\t', index=False)
    return df


def plink_process(*, genes_folder, plink, annotated_vcf, bfiles=None, method=""):
    """
    Use plink to calculate and sum the scores for each gene.

    Parameters
    ----------
    genes_folder : str
        the folder containing the temporary genes files.
    plink : str
        the directory of plink (if not default).
    annotated_vcf : str
        vcf with samples information
    bfiles : str
        the name of the binary files, not needed if annotated file is given.
    method : str
        if empty method is the average. if 'sum', then the variants scores will be summed.

    Returns
    -------

    """
    genes = [line.strip() for line in open(genes_folder + '.genes', 'r')]
    if bfiles:
        input_files = " --bfile " + bfiles
    else:
        input_files = " --vcf " + annotated_vcf
    for gene in tqdm(genes, desc='calculating genes scores'):
        v_file = os.path.join(genes_folder, (gene + '.v'))
        w_file = os.path.join(genes_folder, (gene + '.w'))
        p = subprocess.call(
            plink + input_files + " --double-id" + " --extract " + v_file + " --score " + w_file +
            " 1 2 3 " + method + "--out " + os.path.join(genes_folder, gene), shell=True
        )


def calculate_gbrs(
    *,
    scores_df,
    weights_df,
    weights_col,
    genes_col,
    pathway_file=None,
    samples_col,
    method,
    logger
):
    """
    Calculate a gene-based risk score for each individual in a given dataset.

    Parameters
    ----------
    scores_df : pd.DataFrame
        the matrix with gene-based scores.
    weights_df : pd.DataFrame
        the matrix with weights for each gene.
    weights_col : str
        the name of the column with the weights.
    genes_col : str
        the name of the column with the genes.
    pathway_file : str
        if pathway method chosen, the file contains a list of pathways
    samples_col : str
        the name of the column with samples IDs.
    method : str
        the method used to calulate the gene-based risk scores. [sum | pathways]
    logger
        the logger to log the process

    Returns
    -------
    pd.DataFrame
        a dataframe with the gene-based risk scores.

    """
    genes = set(weights_df[genes_col]).intersection(set(scores_df.columns))
    df = scores_df[genes]
    df = df.reindex(sorted(df.columns), axis=1)
    df *= list(weights_df.sort_values(by=genes_col)[weights_col].values)
    if method == 'sum':
        df['gbrs'] = df.sum(axis=1)
        df = df[['gbrs']]
    elif method == 'pathways':
        pathways = {line.strip().split('\t')[0]: line.strip().split('\t')[2:] for line in open(pathway_file, 'r')}
        pathway_scores = pd.DataFrame(columns=[samples_col] + list(pathways))
        logger.info('calculating pathway scores ...')
        for path, path_genes in tqdm(pathways.items(), desc='calculating pathway scores'):
            selected_genes = list(set(genes) & (set(path_genes)))
            if len(selected_genes) == 0:
                pathway_scores.drop(columns=[path], inplace=True)
            pathway_scores[path] = df[selected_genes].sum(axis=1)
        df = pathway_scores
    return df


def pathway_scoring(
    *,
    pathways,
    genes,
    scores_file,
    samples_col,
    logger
):
    """
    sum the gene-based scores to pathway-based scores.

    Parameters
    ----------
    pathways : dict
        the pathways (keys) and their genes (values)
    genes : str
        the genes to select from the gene-based scores.
    scores_file : str
        the gene-based scores file.
    samples_col : str
        the name the column of the samples IDs.
    logger
        the logger to log the process

    Returns
    -------
    pd.DataFrame
        dataframe containing the pathway-based scores.

    """
    logger.info('reading scores file ...')
    scores_df = pd.read_csv(scores_file, sep=r'\s+', usecols=[samples_col] + genes)
    pathway_scores = pd.DataFrame(columns=[samples_col] + list(pathways))
    pathway_scores[samples_col] = scores_df[samples_col]
    logger.info('calculating pathway scores ...')
    for path, path_genes in tqdm(pathways.items(), desc='calculating pathway scores'):
        selected_genes = list(set(genes) & (set(path_genes)))
        if len(selected_genes) == 0:
            pathway_scores.drop(columns=[path], inplace=True)
        pathway_scores[path] = scores_df[selected_genes].sum(axis=1)
    return pathway_scores
