# -*- coding: utf-8 -*-

import os
import re
import subprocess

import pandas as pd
import numpy as np
from scipy.stats import beta
from tqdm import tqdm


def get_gene_info(
    *,
    annotated_file,
    variant_col,
    af_col,
    del_col,
    output_dir,
    genes_col,
    maf_threshold,
    beta_param,
    weight_func
):
    df = pd.read_csv(annotated_file, usecols=[variant_col, af_col, del_col], sep=r'\s+')
    df = df[df[af_col].values.astype(float) < maf_threshold]
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
    [df[df[genes_col] == gene][[variant_col, 'score', genes_col]].to_csv(os.path.join(output_dir, (str(gene) + '.w')),
        index=False, sep='\t') for gene in tqdm(genes, desc="writing w gene files")]
    [df[df[genes_col] == gene][[variant_col]].to_csv(os.path.join(output_dir, (str(gene) + '.v')),
        index=False, sep='\t') for gene in tqdm(genes, desc="writing v gene files")]
    return output_dir


def combine_scores(
    *,
    input_path,
    output_path,
):
    all_files = [os.path.join(path, name) for path, subdirs, files in os.walk(input_path) for name in files]
    profile_files = [f for f in all_files if re.match(r'.+profile$', f)]
    df = pd.read_csv(str(profile_files[0]), usecols=['IID', 'SCORESUM'], sep=r'\s+').astype({'SCORESUM': np.float32})
    r = re.compile(r'\w+/(.*).profile$')
    gene = r.findall(str(profile_files[0]))
    df.rename(columns={'SCORESUM': gene[0]}, inplace=True)
    pf = profile_files
    for i in tqdm(range(1, len(pf)-1), desc='merging in process'):
        df = unisci(df, pf[i])
    df.to_csv(output_path, sep='\t', index=False)
    return df


def unisci(df, f):
    df2 = pd.read_csv(str(f), usecols=['IID', 'SCORESUM'], sep=r'\s+').astype({'SCORESUM': np.float32})
    r = re.compile(r'\w+/(.*).profile$')
    gene2 = r.findall(str(f))
    df2.rename(columns={'SCORESUM': gene2[0]}, inplace=True)
    df = pd.merge(df, df2, on='IID')
    return df


def plink_process(*, genes_folder, plink, bed, bim, fam):
    genes = [line.strip() for line in open(os.path.join(genes_folder, (genes_folder + '.genes')), 'r')]
    for gene in tqdm(genes, desc='calculating genes scores'):
        v_file = os.path.join(genes_folder, (gene + '.v'))
        w_file = os.path.join(genes_folder, (gene + '.w'))
        p = subprocess.call(
            plink + " --bed " + bed + " --bim " + bim +
            " --fam " + fam + " --extract " + v_file + " --score " + w_file + " 1 2 3 sum --out " +
            os.path.join(genes_folder, gene), shell=True
        )
