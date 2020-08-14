# -*- coding: utf-8 -*-
import os
import subprocess

import pandas as pd
import gzip
from scipy.stats import beta
from tqdm import tqdm


def get_gene_info(*, vcf, output_dir):
    line_number = 0
    with gzip.open(vcf, 'r') as infile:
        for line in infile:
            line_number = line_number + 1
            if line.startswith(b'#CHR'):
                break
    df = pd.read_csv(vcf, usecols=['ID', 'ALT', 'INFO'], sep="\t", skiprows=line_number - 1)
    df = df[df.INFO.str.contains('AF=', regex=True, na=False) & df.INFO.str.contains(
        'RawScore=', regex=True, na=False) & df.INFO.str.contains('gene=', regex=True, na=False)]
    df[['PR', 'AF', 'RawScore', 'PHRED', 'gene']] = df.INFO.str.split(";", expand=True, )
    df.replace(to_replace=r'^AF=', value='', regex=True, inplace=True)
    df.replace(to_replace=r'^RawScore=', value='', regex=True, inplace=True)
    df.replace(to_replace=r'^gene=', value='', regex=True, inplace=True)
    df['beta_AF'] = beta.pdf(df.AF.values.astype(float), 1, 25)
    df = df[df['AF'].values.astype(float) < 0.01]
    df['score'] = df['beta_AF'].values.astype(float) * df['RawScore'].values.astype(float)
    genes = list(set(df['gene']))
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    os.chdir(output_dir)
    gene_file = output_dir + '.genes'
    with open(gene_file, 'w') as f:
        for gene in genes:
            f.write("%s\n" % gene)
    [df[df['gene'] == gene][['ID', 'ALT', 'score', 'gene']].to_csv(str(gene) + '.w', index=False, sep='\t') for gene in
     genes]
    [df[df['gene'] == gene][['ID']].to_csv(str(gene) + '.v', index=False, sep='\t') for gene in genes]
    return output_dir


def plink_process(*, genes_folder, plink, bed, bim, fam):
    genes = [line.strip() for line in open(genes_folder + '.genes', 'r')]
    for gene in tqdm(genes, desc='calculating genes scores'):
        try:
            p = subprocess.call(
                plink + " --bed " + bed + " --bim " + bim +
                " --fam " + fam + " --extract " + gene + ".v --score " + gene + ".w 1 2 3 sum --out " + gene)
        except:
            print('Cannot find ' + gene)
