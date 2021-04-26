# -*- coding: utf-8 -*-

import os
import re
import subprocess

import pandas as pd
import pycaret.classification as cl
import pycaret.regression as pyreg
from pycaret.utils import check_metric
import requests
import numpy as np
from scipy.stats import beta
from tqdm import tqdm
import urllib.request as urllib


def get_gene_info(
    *,
    annotated_file,
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
    :param annotated_file: a file containing the variant, AF, ALT, Gene, and deleterious score.
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
    df = pd.read_csv(annotated_file, usecols=[variant_col, alt_col, af_col, del_col, genes_col], sep=r'\s+')
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
        index=False, sep='\t') for gene in tqdm(genes, desc="writing v gene files")]
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
    r = re.compile(input_path + "(.*).profile$")
    gene = r.findall(str(profile_files[0]))
    df.rename(columns={'SCORESUM': gene[0]}, inplace=True)
    pf = profile_files
    for i in tqdm(range(1, len(pf)-1), desc='merging in process'):
        df = unisci(input_path, df, pf[i])
    df.to_csv(output_path, sep='\t', index=False)
    return df


def unisci(input_path, df, f):
    """
    Merge two dataframes.
    :param input_path: the path to input directory.
    :param df: the main dataframe with all the scores.
    :param f: the file containing the scores of one gene.
    :return: the merged dataframe.
    """
    df2 = pd.read_csv(str(f), usecols=['IID', 'SCORESUM'], sep=r'\s+').astype({'SCORESUM': np.float32})
    r = re.compile(input_path + "(.*).profile$")
    gene2 = r.findall(str(f))
    df2.rename(columns={'SCORESUM': gene2[0]}, inplace=True)
    df = pd.merge(df, df2, on='IID')
    return df


def plink_process(*, genes_folder, plink, bed, bim, fam):
    """
    Use plink to calculate and sum the scores for each gene.
    :param genes_folder: the folder containing the temporary genes files.
    :param plink: the directory of plink (if not default).
    :param bed: the bed file path.
    :param bim: the bim file path.
    :param fam: the fam file path.
    :return:
    """
    genes = [line.strip() for line in open(os.path.join(genes_folder, (genes_folder + '.genes')), 'r')]
    for gene in tqdm(genes, desc='calculating genes scores'):
        v_file = os.path.join(genes_folder, (gene + '.v'))
        w_file = os.path.join(genes_folder, (gene + '.w'))
        p = subprocess.call(
            plink + " --bed " + bed + " --bim " + bim +
            " --fam " + fam + " --extract " + v_file + " --score " + w_file + " 1 2 3 sum --out " +
            os.path.join(genes_folder, gene), shell=True
        )


def get_prs(
    *,
    prs_id,
    vcf=None,
    bed=None,
    bim=None,
    fam=None,
    plink,
    output_file,
):
    """
    Get PGS scores from pgscatalog and calculate PRS score for samples.
    :param prs_id: the PRS ID in the pgscatalog.
    :param vcf: vcf file with sampels info. OR use bed/bim/fam files.
     :param bed: the bed file path.
    :param bim: the bim file path.
    :param fam: the fam file path.
    :param plink:
    :param output_file: the path to the output scores file.
    :return:
    """
    # make sure that the columns are present and matching
    resp = requests.get('https://www.pgscatalog.org/rest/score/%s' % prs_id)
    prs_info = resp.json()
    if resp.status_code != 200 or not prs_info:
        raise Exception('The PRS score might be wrong!')
    url = prs_info['ftp_scoring_file']
    prs_file = prs_id+'.gz'
    urllib.urlretrieve(url, prs_file)
    # unpack file
    if vcf:
        p = subprocess.call(
            plink + " --vcf " + vcf + " --score " + prs_file + " 1 2 3 sum --out " + output_file, shell=True
        )
    else:
        p = subprocess.call(
            plink + " --bed " + bed + " --bim " + bim + " --fam " + fam + " --score " + prs_file + " 1 2 3 sum --out " +
            output_file, shell=True
        )


def create_model(
    *,
    model_name='final_model',
    model_type='reg',
    y_col,
    imbalanced=True,
    normalize=True,
    folds=10,
    training_set,
    testing_set=None,
    test_size=0.25,
    metric=None,
):
    """
    Create a prediction model (classifier or regressor) using the provided dataset.
    :param output_folder: the path to the folder to save outputs.
    :param model_name: the name of the prediction model.
    :param model_type: type of model (reg or classifier).
    :param y_col: the column containing the target (qualitative or quantitative).
    :param imbalanced: True means data is imbalanced.
    :param normalize: True if data needs normalization.
    :param folds: how many folds for cross-validation.
    :param training_set: the training set for the model.
    :param testing_set: if exists an extra evaluation step will be done using the testing set.
    :param test_size: the size to split the training/testing set.
    :param metric: the metric to evaluate the best model.
    :return: the metrics.
    """
    metrics = None
    if model_type == 'reg':
        if not metric:
            metric = 'RMSE'
        reg = pyreg.setup(target=y_col, data=training_set, normalize=normalize, train_size=1-test_size, fold=folds,
                          silent=True)
        best_model = pyreg.compare_models(sort=metric)
        reg_model = pyreg.create_model(best_model)
        reg_tuned_model = pyreg.tune_model(reg_model)
        pyreg.plot_model(reg_tuned_model, save=True)
        pyreg.plot_model(reg_tuned_model, plot='feature', save=True)
        pyreg.plot_model(reg_tuned_model, plot='error', save=True)
        final_model = pyreg.finalize_model(reg_tuned_model)
        if testing_set.bool():
            unseen_predictions = pyreg.predict_model(final_model, data=testing_set)
            r2 = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'R2')
            rmse = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'RMSE')
            metrics = ['R2: ' + r2, 'RMSE: ' + rmse]
        pyreg.save_model(final_model, model_name)
        pyreg.save_experiment('model_session')
        with open('model.log', 'w') as f:
            f.writelines("%s\n" % output for output in reg)
            f.writelines("%s\n" % output for output in metrics)
    elif model_type == 'classifier':
        if not metric:
            metric = 'AUC'
        classifier = cl.setup(target=y_col, fix_imbalance=imbalanced, data=training_set, train_size=1-test_size,
                              silent=True)
        best_model = cl.compare_models(sort=metric)
        cl_model = cl.create_model(best_model)
        cl_tuned_model = pyreg.tune_model(cl_model)
        cl.plot_model(cl_tuned_model, plot='pr', save=True)
        cl.plot_model(cl_tuned_model, plot='confusion_matrix', save=True)
        cl.plot_model(cl_tuned_model, plot='feature', save=True)
        final_model = cl.finalize_model(cl_tuned_model)
        if testing_set.bool():
            unseen_predictions = cl.predict_model(final_model, data=testing_set)
            auc = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'AUC')
            accuracy = check_metric(unseen_predictions[y_col], unseen_predictions.Label, 'accuracy')
            metrics = ['AUC: ' + auc, 'Accuracy: ' + accuracy]
        cl.save_model(final_model, model_name)
        cl.save_experiment('model_session')
        with open('model.log', 'w') as f:
            f.writelines("%s\n" % output for output in classifier)
            f.writelines("%s\n" % output for output in metrics)
    else:
        return Exception('Model requested is not available. Please choose reg or classifier.')
    return final_model
