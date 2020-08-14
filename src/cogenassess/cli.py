# -*- coding: utf-8 -*-
import click

from CoGenAssess.src.cogenassess.pipeline import get_gene_info, plink_process


@click.command()
@click.option('-v', '--vcf', required=True)
@click.option('--bed', required=True)
@click.option('--bim', required=True)
@click.option('--fam', required=True)
@click.option('--plink', default='plink')
@click.option('-o', '--output-dir', required=True)
def score_genes(
    vcf,
    bed,
    bim,
    fam,
    plink,
    output_dir,
):
    click.echo('getting information from vcf files')
    genes_folder = get_gene_info(vcf=vcf)
    plink_process(genes_folder=genes_folder, plink=plink, bed=bed, bim=bim, fam=fam)
    click.echo('process is done.')


if __name__ == '__main__':
    score_genes()
