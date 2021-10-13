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

    :param annotated_vcf: a file containing the variant, AF, ALT, Gene, and deleterious score.
    :param variant_col: the name of the variant column.
    :param af_col: the name of the AF column.
    :param alt_col: the name of the ALT column.
    :param del_col: the name of deleterious score column.
    :param output_dir: directory to save in temporary files.
    :param genes_col: the name of genes column.
    :param maf_threshold: the minor allele frequency threshold, default is 0.01.
    :param beta_param: the parameters of the beta function, if chosen for weighting.
    :param weight_func: the weighting function, beta or log10.

    :return: returns the output directory with all the temporary files.
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

    :param input_path: the directory containing scores files.
    :param output_path: the name of the output file.

    :return: dataframe with all the scores.
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


def plink_process(*, genes_folder, plink, annotated_vcf, bfiles=None):
    """
    Use plink to calculate and sum the scores for each gene.

    :param genes_folder: the folder containing the temporary genes files.
    :param plink: the directory of plink (if not default).
    :param annotated_vcf: vcf with samples information

    :return:
    """
    genes = [line.strip() for line in open(os.path.join(genes_folder, (genes_folder + '.genes')), 'r')]
    if bfiles:
        input_files = " --bfile " + bfiles
    else:
        input_files = " --vcf " + annotated_vcf
    for gene in tqdm(genes, desc='calculating genes scores'):
        v_file = os.path.join(genes_folder, (gene + '.v'))
        w_file = os.path.join(genes_folder, (gene + '.w'))
        p = subprocess.call(
            plink + input_files + " --double-id" + " --extract " + v_file + " --score " + w_file +
            " 1 2 3 sum --out " + os.path.join(genes_folder, gene), shell=True
        )


def calculate_gbrs(
    *,
    scores_df,
    weights_df,
    weights_col,
    genes_col,
    sum=True
):
    """
    Calculate a gene-based risk score for each individual in a given dataset.

    :param scores_df: the matrix with gene-based scores.
    :param weights_df: the matrix with weights for each gene.
    :param weights_col: the name of the column with the weights.
    :param genes_col: the name of the column with the genes.
    :param sum: if True it will sum the gene-based risk scores into one value.

    :return: a dataframe with the gene-based risk scores.
    """
    genes = set(weights_df[genes_col]).intersection(set(scores_df.columns))
    df = scores_df[genes]
    df = df.reindex(sorted(df.columns), axis=1)
    df *= list(weights_df.sort_values(by=genes_col)[weights_col].values)
    if sum:
        df['gbrs'] = df.sum(axis=1)
        df = df[['gbrs']]
    return df


def pathway_scoring(
    *,
    pathway_file,
    output_file,
    scores_file,
):
    pathways = {}
    with open(pathway_file, "r") as file:
        for line in file:
            line = line.strip().split("\t")
            pathways[line[0]] = line[2:]
    all_genes = [item for sublist in list(pathways.values()) for item in sublist]
    scores_df = pd.read_csv(scores_file, sep='\t', usecols=['IID']+all_genes)

