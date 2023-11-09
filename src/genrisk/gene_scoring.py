# -*- coding: utf-8 -*-
import os
import re
import subprocess
import warnings

import numpy as np
import pandas as pd
from scipy.stats import beta
from tqdm import tqdm

from .helpers import uni_profiles


def get_gene_info(
        *,
        annotation_file,
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
    annotation_file : str
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
    df = pd.read_csv(annotation_file, sep=r'\s+', index_col=False, on_bad_lines='warn')
    df.drop_duplicates(inplace=True)
    df[af_col].replace('.', 3.98e-6, inplace=True)  # 1 allele out 125,748 indiv in gnomADexome (251496 alleles)
    df[af_col].replace('nan', 3.98e-6, inplace=True)  # 1 allele out 125,748 indiv in gnomADexome (251496 alleles)
    df[del_col].replace('.', 0, inplace=True)
    df[del_col].replace('nan', 0, inplace=True)
    df = df[df[af_col].values.astype(float) <= maf_threshold]
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
        plink_extension,
):
    """
    Combine the files that contain the scores into one file.

    Parameters
    ----------
    plink_extension: str
        either sscore or profile depending on the plink version
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
    score_files = [f for f in all_files if re.match(r'.+' + plink_extension + '$', f)]
    try:
        df = pd.read_csv(str(score_files[0]), sep=r'\s+', on_bad_lines='warn').iloc[:, [1, -1]]
    except:
        raise Exception("It seems that there is a problem with your score files, please make sure that plink is "
                        "running correctly.")
    scores_col = df.columns[1]
    df.astype({scores_col: np.float32})
    r = re.compile("([a-zA-Z0-9_.-]*)." + plink_extension + "$")
    gene = r.findall(str(score_files[0]))
    df.rename(columns={scores_col: gene[0]}, inplace=True)
    for i in tqdm(range(1, len(score_files)), desc='merging in process'):
        df = uni_profiles(df, score_files[i], plink_extension=plink_extension)
    df.to_csv(output_path, sep='\t', index=False)
    return df


def plink_process(*, genes_folder, plink, bfiles=None, method="", vcf=None):
    """
    Use plink to calculate and sum the scores for each gene.

    Parameters
    ----------
    genes_folder : str
        the folder containing the temporary genes files.
    plink : str
        the directory of plink (if not default).
    bfiles : str
        the name of the binary files, not needed if a vcf is given.
    vcf : str
        the name of vcf file, not needed if binary files are given.
    method : str
        if empty method is the average. if 'sum', then the variants scores will be summed.

    Returns
    -------

    """
    genes = [line.strip() for line in open(genes_folder + '.genes', 'r')]
    if 'plink2' in plink:
        dosage = ""
    else:
        dosage = "double-dosage "
    if bfiles:
        for gene in tqdm(genes, desc='calculating genes scores'):
            v_file = os.path.join(genes_folder, (gene + '.v'))
            w_file = os.path.join(genes_folder, (gene + '.w'))
            try:
                p = subprocess.run(
                    plink + " --bfile " + bfiles + " --double-id --extract " + v_file + " --score " + w_file +
                    " 1 2 3 " + method + dosage + "--out " + os.path.join(genes_folder, gene), shell=True, check=True
                )
            except subprocess.CalledProcessError as e:
                if e.returncode == 127:
                    raise Exception(
                        "Problem running plink, please make sure that you installed plink and provide the correct path")
                else:
                    with open("genes.error", "a") as f:
                        f.write(gene + "\n")
                    warnings.warn(
                        "The score for %s was not calculated because of some issue in the plink process" % gene)
                    continue
    elif vcf:
        for gene in tqdm(genes, desc='calculating genes scores'):
            v_file = os.path.join(genes_folder, (gene + '.v'))
            w_file = os.path.join(genes_folder, (gene + '.w'))
            try:
                p = subprocess.run(
                    plink + " --vcf " + vcf + " --double-id --extract " + v_file + " --score " + w_file +
                    " 1 2 3 " + method + dosage + "--out " + os.path.join(genes_folder, gene), shell=True, check=True
                )
            except subprocess.CalledProcessError as e:
                if e.returncode == 127:
                    raise Exception(
                        "Problem running plink, please make sure that you installed plink and provide the correct path")
                else:
                    with open("genes.error", "a") as f:
                        f.write(gene + "\n")
                    warnings.warn(
                        "The score for %s was not calculated because of some issue in the plink process" % gene)
                    continue
    else:
        raise Exception("Gene scoring failed, no binary or VCF files are provided!")


def calculate_gbrs(
        *,
        scores_df,
        weights_df,
        weights_col,
        genes_col,
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

    Returns
    -------
    pd.DataFrame
        a dataframe with the gene-based risk scores.

    """
    genes = set(weights_df[genes_col]).intersection(set(scores_df.columns))
    df = scores_df[genes]
    df = df.reindex(sorted(df.columns), axis=1)
    df *= list(weights_df.sort_values(by=genes_col)[weights_col].values)
    df['gbrs'] = df.sum(axis=1)
    df = df[['gbrs']]
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
    genes : List
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
    scores_df = pd.read_csv(scores_file, sep=r'\s+', usecols=[samples_col] + genes, on_bad_lines='warn')
    pathway_scores = pd.DataFrame(columns=[samples_col] + list(pathways))
    pathway_scores[samples_col] = scores_df[samples_col]
    logger.info('calculating pathway scores ...')
    for path, path_genes in tqdm(pathways.items(), desc='calculating pathway scores'):
        selected_genes = list(set(genes) & (set(path_genes)))
        if len(selected_genes) == 0:
            pathway_scores.drop(columns=[path], inplace=True)
        pathway_scores[path] = scores_df[selected_genes].sum(axis=1) / len(selected_genes)
    return pathway_scores
