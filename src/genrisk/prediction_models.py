import joblib
import numpy as np
import pycaret.classification as pycl
import pycaret.regression as pyreg
import sklearn.metrics as metrics

from .helpers import write_output, generate_confusion_matrix, generate_scatterplot


def regression_model(
    *,
    y_col,
    training_set,
    normalize,
    folds,
    metric,
    model_name,
    testing_set,
    imbalanced,
    seed,
    include_models,
    normalize_method,
    feature_selection
):
    """
    Build a regression model for prediction.

    Parameters
    ----------
    y_col : str
        the name of the target column.
    training_set : pd.DataFrame
        DataFrame containing the training data.
    normalize : bool
        if True the dataset will be normalized before training.
    test_size : float
        Between [0.0-1.0]. The size of the split for test within the training set.
    folds : int
        number of folds for cross validation.
    metric : str
        the metric used for evaluating the best model.
    model_name : str
        the name to save the model.
    testing_set : pd.DataFrame
        the external dataset for evaluating the best model.
    imbalanced
    seed : int
        random number to initilize the process.
    include_models : List
        a list of models to be included in the process.
    normalize_method : str
        The method used for normalizing the data.

    Returns
    -------
    Final regression model

    """
    if not metric:
        metric = 'RMSE'
    setup = pyreg.setup(target=y_col, data=training_set, normalize=normalize, normalize_method=normalize_method,
                        test_data=testing_set, fold=folds, silent=True, session_id=seed,
                        feature_selection=feature_selection)
    best_model = pyreg.compare_models(sort=metric, include=include_models)
    pyreg.pull().to_csv(model_name + '_compare_models.tsv', sep='\t', index=False)
    reg_model = pyreg.create_model(best_model)
    reg_tuned_model = pyreg.tune_model(reg_model, optimize=metric)
    pyreg.pull().to_csv(model_name + '_tuned_model.tsv', sep='\t', index=False)
    final_model = pyreg.finalize_model(reg_tuned_model)
    pyreg.plot_model(final_model, save=True)
    pyreg.plot_model(final_model, plot='feature', save=True)
    pyreg.plot_model(final_model, plot='error', save=True)
    pyreg.save_model(final_model, model_name)
    if len(testing_set.index) != 0:
        unseen_predictions = test_regressor(
            model_path=model_name+'.pkl', x_set=testing_set.drop(columns=[y_col]), y_col=testing_set[y_col], output=model_name
        )
        unseen_predictions.to_csv(model_name + '_external_testing_results.tsv', sep='\t', index=True)
    return final_model


def classification_model(
    *,
    y_col,
    training_set,
    normalize,
    folds,
    metric,
    model_name,
    testing_set,
    imbalanced,
    seed,
    include_models,
    normalize_method,
    feature_selection
):
    """
    Build a classification model for prediction.

    Parameters
    ----------
    y_col : str
        the name of the target column.
    training_set : pd.DataFrame
        DataFrame containing the training data.
    normalize : bool
        if True the dataset will be normalized before training.
    test_size : float
        Between [0.0-1.0]. The size of the split for test within the training set.
    folds : int
        number of folds for cross validation.
    metric : str
        the metric used for evaluating the best model.
    model_name : str
        the name to save the model.
    testing_set : pd.DataFrame
        the external dataset for evaluating the best model.
    imbalanced : bool
        if True the imbalance will be fixed before the training.
    seed : int
        random number to initilize the process.
    include_models : List
        a list of models to be included in the process.
    normalize_method : str
        The method used for normalizing the data.

    Returns
    -------
    Final classification model

    """
    if not metric:
        metric = 'AUC'
    if len(training_set.groupby(y_col).groups) == 2:
        training_set[y_col] = np.interp(
            training_set[y_col], (training_set[y_col].min(), training_set[y_col].max()), (0, 1))
        testing_set[y_col] = np.interp(
            testing_set[y_col], (testing_set[y_col].min(), testing_set[y_col].max()), (0, 1))
    setup = pycl.setup(target=y_col, fix_imbalance=imbalanced, normalize=normalize, normalize_method=normalize_method,
                    data=training_set, test_data=testing_set, silent=True, fold=folds, session_id=seed,
                    feature_selection=feature_selection)
    best_model = pycl.compare_models(sort=metric, include=include_models)
    pycl.pull().to_csv(model_name + '_compare_models.tsv', sep='\t', index=False)
    cl_model = pycl.create_model(best_model)
    cl_tuned_model = pycl.tune_model(cl_model, optimize=metric)
    pycl.pull().to_csv(model_name + '_tuned_model.tsv', sep='\t', index=False)
    final_model = pycl.finalize_model(cl_tuned_model)
    pycl.plot_model(final_model, plot='pr', save=True)
    pycl.plot_model(final_model, plot='feature', save=True)
    pycl.save_model(final_model, model_name)
    unseen_predictions = test_classifier(
        model_path=model_name+'.pkl', x_set=testing_set.drop(columns=[y_col]), y_col=testing_set[y_col], output=model_name
    )
    unseen_predictions.to_csv(model_name + '_external_testing_results.tsv', sep='\t', index=True)
    return final_model


def test_classifier(
    *,
    model_path,
    x_set,
    y_col,
    output
):
    """
    Test a model with an independent dataset.

    Parameters
    ----------
    model_path
        the model used for prediction
    x_set : pd.DataFrame
        the independent testing set (without the target)
    y_col : pd.Series
        the target true values.
    output : str
        the name for the evaluation output

    Returns
    -------
    pd.DataFrame
        the testing dataset with the true and predicted values.
    """
    model = joblib.load(model_path)
    x_set['Label'] = model.predict(x_set)
    x_set['True value'] = y_col
    if len(x_set.groupby('True value').groups) == 2:
        x_set['True value'] = np.interp(
            x_set['True value'], (x_set['True value'].min(), x_set['True value'].max()), (0, 1))
    report = metrics.classification_report(x_set['True value'], x_set['Label'])
    acc = metrics.accuracy_score(x_set['True value'], x_set['Label'])
    auc = metrics.roc_auc_score(x_set['True value'], x_set['Label'])
    input_list = [
        output, '\nTesting model report: \n', report + '\n', 'AUC = ' + str(auc) + '\n', 'Accuracy = ' + str(acc) + '\n'
    ]
    write_output(input_list=input_list, output=output + "_report.txt")
    generate_confusion_matrix(y_true=x_set['True value'], y_pred=x_set['Label'], output=output)
    return x_set


def test_regressor(
    *,
    model_path,
    x_set,
    y_col,
    output
):
    """
    Test a model with an independent dataset.

    Parameters
    ----------
    model_path
        the model used for prediction
    x_set : pd.DataFrame
        the independent testing set (without the target)
    y_col : pd.Series
        the target true values.
    output : str
        the name for the evaluation output

    Returns
    -------
    pd.DataFrame
        the testing dataset with the true and predicted values.

    """
    model = joblib.load(model_path)
    x_set['Label'] = model.predict(x_set)
    x_set['True value'] = y_col
    r2 = metrics.r2_score(y_col, x_set['Label'])
    rmse = metrics.mean_squared_error(y_col, x_set['Label'], squared=False)
    plot = generate_scatterplot(
        x_axis=x_set['Label'], y_axis=y_col, output=output)
    input_list = [output + '\nTesting model report: \n', 'R^2 = ' + str(r2) + '\n', 'RMSE = ' + str(rmse) + '\n']
    write_output(input_list=input_list, output=output + "_report.txt")
    return x_set
