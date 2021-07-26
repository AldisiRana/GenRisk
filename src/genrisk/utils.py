# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from pybiomart import Dataset
from qmplot import qqplot
from tqdm import tqdm


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
        genes_df = pd.read_csv(genes_lengths_file, sep='\t')
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
    """
    Generate  manhattan plot from a given dataset.

    :param data: a dataframe with pvalues and gene information.
    :param chr_col: the column with the chromosomes.
    :param pos_col: the column containing the position/start.
    :param pvals_col: the column containing the p_values.
    :param genes_col: the column containing gene names.
    :param manhattan_output: the path to output the manhattan plot image.

    :return:
    """
    data.drop_duplicates(subset=[genes_col], inplace=True)
    data.dropna(subset=[pvals_col], inplace=True)
    data['-logp'] = - np.log10(data[pvals_col])
    data = data.dropna(subset=[genes_col, chr_col, pos_col], axis=0)
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
                     horizontalalignment='left', size='medium', rotation=20, color='black'))
    adjust_text(anno, arrowprops=dict(arrowstyle="->", color='b', lw=0.5))
    chrom_df = data.groupby(chr_col)['i'].median()
    plot.ax.set_xlabel(chr_col)
    plot.ax.set_xticks(chrom_df)
    plot.ax.set_xticklabels(chrom_df.index)
    plot.ax.axhline(-np.log10(1.00e-05), c='blue', ls='--')
    plot.ax.axhline(-np.log10(5.00e-08), c='red', ls='--')
    plot.savefig(manhattan_output)


def draw_qqplot(*, pvals, qq_output):
    """
    Generate QQ-plot for given data.

    :param pvals: the list of p_values.
    :param qq_output: the path to output the QQplot image.

    :return:
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

