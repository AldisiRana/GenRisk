# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from qmplot import qqplot
from tqdm import tqdm
from pybiomart import Dataset


def gene_length_normalize(*, genes_info, genes_col='HGNC symbol', length_col='gene_length', scores_df, samples_col):
    """
    Normalize dataset by gene length. if gene lengths file is not provided, info will be retrieved from ensembl.

    Parameters
    ----------
    genes_info : str
        file containing gene lengths. If file is not provided, info will be retrieved from ensembl
    genes_col : str
        column containing gene names.
    length_col : str
        column containing the length of each gene.
    scores_df : pd.DataFrame
        dataframe containing data to normalize.
    samples_col : str
        column containing samples IDs.

    Returns
    -------
    pd.Dataframe
        dataframe containing normalized dataframe.

    """
    unnormalized = []
    if not genes_info:
        dataset = Dataset(name='hsapiens_gene_ensembl',
                          host='http://www.ensembl.org')
        genes_df = dataset.query(attributes=['hgnc_symbol', 'start_position', 'end_position'])
        genes_df['gene_length'] = genes_df['Gene end (bp)'] - genes_df['Gene start (bp)']
    else:
        genes_df = pd.read_csv(genes_info, sep='\t')
    genes_lengths = genes_df.set_index(genes_col).to_dict()[length_col]
    for (name, data) in tqdm(scores_df.drop(columns=[samples_col]).iteritems(), desc="Normalizing genes scores"):
        if name not in genes_lengths.keys():
            unnormalized.append(name)
            continue
        # normalize genes by length
        scores_df[name] = round(scores_df[name] / genes_lengths[name], 5)
    scores_df = scores_df.drop(unnormalized, axis=1)
    return scores_df


def draw_manhattan(*, data, chr_col, pos_col, pvals_col, genes_col, manhattan_output):
    """
    Generate  manhattan plot from a given dataset.

    Parameters
    ----------
    data : pd.DataFrame
        a dataframe with pvalues and gene information.
    chr_col : str
        the column with the chromosomes.
    pos_col : str
        the column containing the position/start.
    pvals_col : str
        the column containing the position/start.
    genes_col : str
        the column containing gene names.
    manhattan_output : str
        the path to output the manhattan plot image.

    Returns
    -------
    Manhattan plot

    """
    data.drop_duplicates(subset=[genes_col], inplace=True)
    data.dropna(subset=[pvals_col], inplace=True)
    data['-logp'] = - np.log10(data[pvals_col])
    data = data.dropna(subset=[genes_col, chr_col, pos_col], axis=0)
    data[chr_col].replace('X', 23, inplace=True)
    data[chr_col] = data[chr_col].astype('int64')
    data = data.sort_values([chr_col, pos_col])
    data.reset_index(inplace=True, drop=True)
    data['ind'] = data.index
    df_grouped = data.groupby((chr_col))
    fig = plt.figure(figsize=(14, 8))  # Set the figure size
    ax = fig.add_subplot(111)
    colors = ['black', 'grey']
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='-logp', color=colors[num % len(colors)], ax=ax)
    ax.set_ylim(0.0, max(data["-logp"]) + 1)

    for i in data[data['-logp'] > 5.6].index:
        ax.annotate(xy=(data.ind[i], data['-logp'][i]),
                              xytext=(data.ind[i] + 0.2, data['-logp'][i] + 0.2),
                              text=data[genes_col][i],
                              arrowprops=dict(arrowstyle="->", color='b', lw=0.5),
                              horizontalalignment='left', size='medium', rotation=20, color='black')
    # adjust_text(anno, arrowprops=dict(arrowstyle="->", color='b', lw=0.5))
    chrom_df = data.groupby(chr_col)['ind'].median()
    data['color group'] = data[chr_col].apply(lambda x: 'A' if x % 2 == 0 else 'B')
    fig.suptitle(manhattan_output.split('.')[0])
    ax.set_ylim(0.0, max(data["-logp"]) + 1)
    ax.set_xlabel(chr_col)
    ax.set_xticks(chrom_df)
    ax.set_xticklabels(chrom_df.index)
    ax.axhline(-np.log10(1.00e-05), c='blue', ls='--')
    ax.axhline(-np.log10(5.00e-08), c='red', ls='--')
    fig.savefig(manhattan_output)
    return ax


def draw_qqplot(*, pvals, qq_output):
    """
    Generate QQ-plot for given data.

    Parameters
    ----------
    pvals : pd.Series
        the list of p_values.
    qq_output : str
        the path to output the QQplot image.

    Returns
    -------
    QQPlot

    """
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
    return ax


def merge_files(*, files_lst, sep, by, cols=None):
    """
    Merge a list of files with the same format into one dataframe.

    Parameters
    ----------
    files_lst
        a list with files to merge.
    sep
        the column seperator.
    by
        the common column between all files.
    cols
        selected columns to remain in dataframe.

    Returns
    -------
    merged dataframe

    """
    df = pd.read_csv(files_lst[0], sep=sep)
    for file in files_lst[1:]:
        df2 = pd.read_csv(file, sep=sep)
        df = pd.merge(df, df2, on=by)
    if cols:
        df = df[[cols]]
    return df

