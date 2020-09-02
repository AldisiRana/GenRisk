# -*- coding: utf-8 -*-
import click

import pandas as pd

from .utils import normalize_gene_len, merge_matrices, find_pvalue
from .pipeline import get_gene_info, plink_process, combine_scores


@click.group()
def main():
    """Handle cogenassess functions."""


@main.command()
@click.option('-v', '--vcf', required=True)
@click.option('--bed', required=True)
@click.option('--bim', required=True)
@click.option('--fam', required=True)
@click.option('--plink', default='plink')
@click.option('-t', '--temp-dir', required=True)
@click.option('-o', '--output-file', required=True)
@click.option('--beta-param', default=(1.0, 25.0), nargs=2, type=float)
@click.option('--weight-func', default='beta', type=click.Choice(['beta', 'log10']))
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
):
    # check number of processes
    click.echo('getting information from vcf files')
    genes_folder = get_gene_info(vcf=vcf, output_dir=temp_dir, beta_param=beta_param, weight_func=weight_func)
    click.echo('calculating gene scores ...')
    plink_process(genes_folder=genes_folder, plink=plink, bed=bed, bim=bim, fam=fam)
    click.echo('combining score files ...')
    df = combine_scores(input_path=temp_dir, output_path=output_file)
    click.echo(df.info())
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
@click.option('--scores-file-sep', default='\t', help="the seperator for scores files")
@click.option('-i', '--genotype-file', required=True, help="File containing information about the cohort.")
@click.option('--genotype-file-sep', default='\t', help="The file separator")
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
    scores_file_sep,
    genotype_file_sep,
    adj_pval,
):
    """Calculate the P-value between two given groups."""
    scores_df = pd.read_csv(scores_file, sep=scores_file_sep)

    click.echo("The process for calculating the p_values will start now.")
    df = find_pvalue(
        scores_df=scores_df,
        output_file=output_path,
        genotype_file=genotype_file,
        genotype_file_sep=genotype_file_sep,
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
@click.option('-d', '--directory', required=True, help="The directory that contains the matrices to merge.")
@click.option('-s', '--file-suffix', default='.tsv', help='The suffix of scores files in directory')
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('--samples-col', required=True, multiple=True, help="the name of samples column in matrices")
@click.option('--scores-col', required=True, help="the name of scores column in matrices")
@click.option('--file-sep', default='\t', help="the seperator for scores files")
def merge(
    *,
    directory,
    output_path,
    samples_col,
    scores_col,
    file_sep,
    file_suffix
):
    """This command merges all matrices in a directory into one big matrix"""
    click.echo("Starting the merging process")
    merge_matrices(
        directory=directory,
        file_suffix=file_suffix,
        output_path=output_path,
        scores_col=scores_col,
        file_sep=file_sep,
        samples_col=list(samples_col),
    )
    click.echo("Merging is done.")


@main.command()
@click.option('-m', '--matrix-file', required=True, help="The scoring matrix to normalize.")
@click.option('-g', '--genes-lengths-file',
              help="The file containing the lengths of genes. If not provided it will be produced.")
@click.option('-o', '--output-path', required=True, help='the path for the output file.')
@click.option('--file-sep', default='\t', help="the seperator for scores files")
@click.option('-s', '--samples-col', default='patient_id', help='the name of the samples column')
def normalize(
    *,
    matrix_file,
    genes_lengths_file=None,
    output_path=None,
    file_sep='\t',
    samples_col
):
    """This command normalizes the scoring matrix by gene length."""
    click.echo("Normalization in process.")
    normalize_gene_len(
        matrix_file=matrix_file,
        file_sep=file_sep,
        genes_lengths_file=genes_lengths_file,
        output_path=output_path,
        samples_col=samples_col
    )


if __name__ == '__main__':
    main()
