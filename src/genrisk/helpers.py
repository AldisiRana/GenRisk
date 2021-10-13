# -*- coding: utf-8 -*-
from datetime import datetime
import re

import logging
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn.metrics as metrics


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


def generate_scatterplot(*, x_axis, y_axis, output):
    """

    Generate a scatterplot.
    :param x_axis: the x axis column.
    :param y_axis: the y axis column.
    :param output: the output image path.
    :return: a scatterplot
    """
    plt.scatter(x_axis, y_axis, alpha=0.5)
    m, b = np.polyfit(x_axis, y_axis, 1)
    plt.plot(x_axis, m * x_axis + b, 'r')
    plt.title('Actual vs predicted scatterplot')
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.savefig(output + '_regressor_scatterplot.png')
    return plt.show()


def generate_confusion_matrix(x_set, y_set, output):
    """
    Generate a confusion matrix for a dataset.

    :param x_set: the features set.
    :param y_set: the target set.
    :param output: the output image path
    :return: confusion matrix plot
    """
    confusion = metrics.plot_confusion_matrix(x_set, y_set)
    confusion.ax_.set_title('Classifier confusion matrix')
    plt.savefig(output + '_classifier_confusion_matrix.png')
    return plt.show()


def write_output(*, input_list, output):
    """
    Write list info a txt file.

    :param input_list: a list of lines
    :param output: the output file
    """
    textfile = open(output, "w")
    for line in input_list:
        textfile.write(line)
    textfile.close()


def create_logger():
    """
    Create a logger to output command information.

    :return: logger
    """
    logger = logging.getLogger('genrisk')
    logger.setLevel(logging.DEBUG)
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    now = datetime.now()
    log_file = 'genrisk_' + now.strftime("%d%m%Y_%H%M") + '.log'
    # create file handler which logs even debug messages
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger
