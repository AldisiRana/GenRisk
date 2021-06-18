# -*- coding: utf-8 -*-

import gzip
import os
import re
import subprocess
import urllib.request as urllib

from adjustText import adjust_text
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from qmplot import qqplot
import requests
from pybiomart import Dataset
import seaborn as sns
from scipy.stats import beta, pearsonr
from tqdm import tqdm


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
    with open(os.path.join(output_dir, gene_file), 'w') as f:
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


def download_pgs(
    *,
    prs_id,
):
    """
    Get PGS from pgscatalog

    :param prs_id: the PRS ID in the pgscatalog.
    :return:
    """
    # make sure that the columns are present and matching
    resp = requests.get('https://www.pgscatalog.org/rest/score/%s' % prs_id)
    prs_info = resp.json()
    if resp.status_code != 200 or not prs_info:
        raise Exception('The PRS score might be wrong!')
    url = prs_info['ftp_scoring_file']
    prs_file = prs_id + '.gz'
    urllib.urlretrieve(url, prs_file)
    return prs_file


def calc_corr(
    *,
    first_file,
    second_file,
    samples_col,
    output_file,
):
    """
    Calculate the pearson's correlation between same genes in two scoring matices.

    :param first_file: the path to the first scores file.
    :param second_file: the path to the second scores file.
    :param samples_col: the column containing the samples IDs.
    :param output_file: the path to the output file with correlation values.
    :return:
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
    first_df = pd.read_csv(first_file, sep=r'\s+', index_col=False)
    second_df = pd.read_csv(second_file, sep=r'\s+', index_col=False)
    for gene in tqdm(genes, desc='calculating correlation'):
        gene_df = pd.merge(first_df[[samples_col, gene]], second_df[[samples_col, gene]], on=samples_col)
        gene_df.replace([np.inf, -np.inf, np.nan], 0.0, inplace=True)
        corr, pval = pearsonr(gene_df[gene + '_x'], gene_df[gene + '_y'])
        corr_info.append([gene, corr, pval])
    corr_df = pd.DataFrame(corr_info, columns=['genes', 'corr', 'p_value']).sort_values(by=['p_value'])
    corr_df.to_csv(output_file, sep='\t', index=False)
    return corr_df


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
    :param samples_col: column containing the samples IDs.
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


def draw_manhattan(*, data, chr_col, pos_col, pvals_col, genes_col, manhattan_output):
    data.drop_duplicates(subset=[genes_col], inplace=True)
    data['-logp'] = - np.log10(data[pvals_col])
    data = data.dropna(how="any", axis=0)
    data[chr_col].replace('X', 23, inplace=True)
    data[chr_col] = data[chr_col].astype('int64')
    data = data.sort_values([chr_col, pos_col])
    data.reset_index(inplace=True, drop=True)
    data['i'] = data.index
    # Generate Manhattan plot: (#optional tweaks for relplot: linewidth=0, s=9)
    sns.set_style("white")
    custom_palette = sns.color_palette(
        ["#000000", "#808080", "#000000", "#808080", "#000000", "#808080", "#000000", "#808080", "#000000", "#808080",
         "#000000", "#808080", "#000000", "#808080", "#000000", "#808080", "#000000", "#808080", "#000000", "#808080",
         "#000000", "#808080", "#000000"])
    plot = sns.relplot(data=data, x='i', y='-logp', aspect=3.7,
                       hue=chr_col, palette=custom_palette, kind='scatter', legend=None)
    plot.fig.suptitle(manhattan_output.split('.')[0])
    plot.ax.set_ylim(0.0, max(data["-logp"]) + 1)

    anno = []
    for ind in data.nlargest(10, ['-logp']).index:
        anno.append(
            plt.text(data.i[ind] + 0.2, data['-logp'][ind] + 0.2, data[genes_col][ind],
                     horizontalalignment='left', size='medium',color='black'))
    adjust_text(anno, only_move={'points': 'y', 'texts': 'y'}, arrowprops=dict(arrowstyle="->", color='b', lw=0.5))
    chrom_df = data.groupby(chr_col)['i'].median()
    plot.ax.set_xlabel(chr_col)
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.axhline(-np.log10(1.00e-05), c='blue', ls='--')
    plot.ax.axhline(-np.log10(5.00e-08), c='red', ls='--')
    plot.savefig(manhattan_output)


def draw_qqplot(*, pvals, qq_output):
    pvals.dropna(inplace=True)
    f, ax = plt.subplots(figsize=(6, 6), facecolor="w", edgecolor="k")
    qqplot(data=pvals,
           marker="o",
           title=qq_output.split('.')[0],
           xlabel=r"Expected $-log_{10}{(P)}$",
           ylabel=r"Observed $-log_{10}{(P)}$",
           dpi=300,
           figname=qq_output,
           ax=ax)

