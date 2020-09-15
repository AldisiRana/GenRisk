# -*- coding: utf-8 -*-
import re

import click
import shutil

import pandas as pd
from scipy.stats import pearsonr
from tqdm import tqdm

from .utils import normalize_gene_len, find_pvalue
from .pipeline import get_gene_info, plink_process, combine_scores


@click.group()
def main():
    """Handle cogenassess functions."""


@main.command()
@click.option('-v', '--vcf', required=True, help='the annotated vcf file')
@click.option('--bed', required=True)
@click.option('--bim', required=True)
@click.option('--fam', required=True)
@click.option('--plink', default='plink')
@click.option('-t', '--temp-dir', required=True)
@click.option('-o', '--output-file', required=True)
@click.option('--beta-param', default=(1.0, 25.0), nargs=2, type=float)
@click.option('--weight-func', default='beta', type=click.Choice(['beta', 'log10']))
@click.option('--remove-temp', is_flag=True)
def score_genes(
    vcf,
    bed,
    bim,
    fam,
    plink,
    beta_param,
    temp_dir,
    output_file,
    weight_func,
    remove_temp,
):
    # check number of processes
    click.echo('getting information from vcf files')
    genes_folder = get_gene_info(vcf=vcf, output_dir=temp_dir, beta_param=beta_param, weight_func=weight_func)
    click.echo('calculating gene scores ...')
    plink_process(genes_folder=genes_folder, plink=plink, bed=bed, bim=bim, fam=fam)
    click.echo('combining score files ...')
    df = combine_scores(input_path=temp_dir, output_path=output_file)
    click.echo(df.info())
    if remove_temp:
        shutil.rmtree(temp_dir)
    click.echo('process is complete.')


@main.command()
@click.option('--bed', required=True)
@click.option('--bim', required=True)
@click.option('--fam', required=True)
@click.option('--plink', default='plink')
@click.option('--genes-folder', required=True)
def run_plink(*, genes_folder, plink, bed, bim, fam):
    click.echo('staring plink processing ...')
    plink_process(genes_folder=genes_folder, plink=plink, bed=bed, bim=bim, fam=fam)
    click.echo('plink processing is complete.')


@main.command()
@click.option('-s', '--scores-file', required=True, help="The scoring file of genes across a population.")
@click.option('-i', '--genotype-file', required=True, help="File containing information about the cohort.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('-g', '--genes',
              help="a list containing the genes to calculate. if not provided all genes will be used.")
@click.option('-t', '--test', required=True, type=click.Choice(['ttest_ind', 'mannwhitneyu', 'logit', 'glm']),
              help='statistical test for calculating P value.')
@click.option('-c', '--cases-column', required=True, help="the name of the column that contains the case/control type.")
@click.option('-m', '--samples-column', required=True, help="the name of the column that contains the samples.")
@click.option('-p', '--pc-file', default=None, help="Principle components values for logistic regression.")
@click.option('--adj-pval', type=click.Choice(
    ['bonferroni', 'sidak', 'holm-sidak', 'holm',
     'simes-hochberg', 'hommel', 'fdr_bh', 'fdr_by', 'fdr_tsbh', 'fdr_tsbky']))
def calculate_pval(
    *,
    scores_file,
    genotype_file,
    output_path,
    genes,
    cases_column,
    samples_column,
    test,
    pc_file,
    adj_pval,
):
    """Calculate the P-value between two given groups."""
    scores_df = pd.read_csv(scores_file, sep=r'\s+')

    click.echo("The process for calculating the p_values will start now.")
    df = find_pvalue(
        scores_df=scores_df,
        output_file=output_path,
        genotype_file=genotype_file,
        genes=genes,
        cases_column=cases_column,
        samples_column=samples_column,
        test=test,
        pc_file=pc_file,
        adj_pval=adj_pval,
    )
    click.echo('Process is complete.')
    click.echo(df.info())


@main.command()
@click.option('-i', '--input-path', required=True, help="The directory that contains the matrices to merge.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('-r', '--remove-input', is_flag=True, help='if flagged will remove input folder')
def merge(
    *,
    output_path,
    input_path,
    remove_input,
):
    """This command merges all matrices in a directory into one big matrix"""
    click.echo("Starting the merging process")
    df = combine_scores(input_path=input_path, output_path=output_path)
    click.echo(df.info())
    if remove_input:
        shutil.rmtree(input_path)
    click.echo("Merging is done.")


@main.command()
@click.option('-m', '--matrix-file', required=True, help="The scoring matrix to normalize.")
@click.option('-g', '--genes-lengths-file',
              help="The file containing the lengths of genes. If not provided it will be produced.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('-s', '--samples-col', default='IID', help='the name of the samples column')
def normalize(
    *,
    matrix_file,
    genes_lengths_file=None,
    output_path=None,
    samples_col
):
    """This command normalizes the scoring matrix by gene length."""
    click.echo("Normalization in process.")
    normalize_gene_len(
        matrix_file=matrix_file,
        genes_lengths_file=genes_lengths_file,
        output_path=output_path,
        samples_col=samples_col
    )


@main.command()
@click.option('--first-file', required=True)
@click.option('--second-file', required=True)
@click.option('--samples-col', default='IID')
@click.option('--output-file', required=True)
def calc_corr(
    *,
    first_file,
    second_file,
    samples_col,
    output_file,
):
    """Calculate the pearson's correlation between same genes in two scoring matices."""
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
    for gene in tqdm(genes, desc='calculating correlation'):
        if gene not in second_df.columns:
            continue
        first_df = pd.read_csv(first_file, sep=r'\s+', index_col=False, usecols=[samples_col, gene])
        second_df = pd.read_csv(second_file, sep=r'\s+', index_col=False, usecols=[samples_col, gene])
        gene_df = pd.merge(first_df, second_df, on='IID')
        corr, pval = pearsonr(gene_df[gene+'_x'], gene_df[gene+'_y'])
        corr_info.append([gene, corr, pval])
    corr_df = pd.DataFrame(corr_info, columns=['genes', 'corr', 'p_value']).sort_values(by=['p_value'])
    corr_df.to_csv(output_file, sep='\t', index=False)
    click.echo('Process is complete.')
    click.echo(corr_df.info())


if __name__ == '__main__':
    main()
