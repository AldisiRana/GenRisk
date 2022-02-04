
import pycaret.classification as pycl
import pycaret.regression as pyreg
import sklearn.metrics as metrics

from genrisk.helpers import write_output, generate_confusion_matrix, generate_scatterplot


def regression_model(
    *,
    y_col,
    training_set,
    normalize,
    test_size,
    folds,
    metric,
    model_name,
    testing_set,
    imbalanced,
    seed,
    include_models,
    normalize_method
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
                        train_size=1 - test_size, fold=folds, silent=True, session_id=seed)
    best_model = pyreg.compare_models(sort=metric, include=include_models)
    pyreg.pull().to_csv(model_name + '_compare_models.tsv', sep='\t', index=False)
    reg_model = pyreg.create_model(best_model)
    reg_tuned_model = pyreg.tune_model(reg_model, optimize=metric)
    pyreg.pull().to_csv(model_name + '_tuned_model.tsv', sep='\t', index=False)
    final_model = pyreg.finalize_model(reg_tuned_model)
    pyreg.plot_model(final_model, save=True)
    pyreg.plot_model(final_model, plot='feature', save=True)
    pyreg.plot_model(final_model, plot='error', save=True)
    if len(testing_set.index) != 0:
        unseen_predictions = test_regressor(
            model=final_model, x_set=testing_set.drop(columns=[y_col]), y_col=testing_set[y_col], output=model_name
        )
        unseen_predictions.to_csv(model_name + '_external_testing_results.tsv', sep='\t', index=True)
    pyreg.save_model(final_model, model_name)
    return final_model


def classification_model(
    *,
    y_col,
    training_set,
    normalize,
    test_size,
    folds,
    metric,
    model_name,
    testing_set,
    imbalanced,
    seed,
    include_models,
    normalize_method,
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
    setup = pycl.setup(target=y_col, fix_imbalance=imbalanced, normalize=normalize, normalize_method=normalize_method,
                       data=training_set, train_size=1 - test_size, silent=True, fold=folds, session_id=seed)
    best_model = pycl.compare_models(sort=metric, include=include_models)
    pycl.pull().to_csv(model_name + '_compare_models.tsv', sep='\t', index=False)
    cl_model = pycl.create_model(best_model)
    cl_tuned_model = pycl.tune_model(cl_model, optimize=metric)
    pycl.pull().to_csv(model_name + '_tuned_model.tsv', sep='\t', index=False)
    final_model = pycl.finalize_model(cl_tuned_model)
    pycl.plot_model(final_model, plot='pr', save=True)
    pycl.plot_model(final_model, plot='confusion_matrix', save=True)
    pycl.plot_model(final_model, plot='feature', save=True)
    if len(testing_set.index) != 0:
        unseen_predictions = test_classifier(
            model=final_model, x_set=testing_set.drop(columns=[y_col]), y_col=testing_set[y_col], output=model_name
        )
        unseen_predictions.to_csv(model_name + '_external_testing_results.tsv', sep='\t', index=True)
    pycl.save_model(final_model, model_name)
    return final_model


def test_classifier(
    *,
    model,
    x_set,
    y_col,
    output
):
    """
        Test a model with an independent dataset.

        Parameters
        ----------
        model
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
    x_set['Label'] = model.predict(x_set)
    x_set['True value'] = y_col
    report = metrics.classification_report(y_col, x_set.Label)
    acc = metrics.accuracy_score(y_col, x_set.Label)
    auc = metrics.auc(y_col, x_set.Label)
    plot = generate_confusion_matrix(x_set=x_set, y_set=y_col, output=output)
    input_list = [
        output, '\nTesting model report: \n', report + '\n', 'AUC = ' + str(auc) + '\n', 'Accuracy = ' + str(acc) + '\n'
    ]
    write_output(input_list=input_list, output=output + "_report.txt")
    return x_set


def test_regressor(
    *,
    model,
    x_set,
    y_col,
    output
):
    """
    Test a model with an independent dataset.

    Parameters
    ----------
    model
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
    x_set['Label'] = model.predict(x_set)
    x_set['True value'] = y_col
    r2 = metrics.r2_score(y_col, x_set.Label)
    rmse = metrics.mean_squared_error(y_col, x_set.Label, squared=False)
    plot = generate_scatterplot(
        x_axis=x_set.Label, y_axis=y_col, output=output)
    input_list = [output + '\nTesting model report: \n', 'R^2 = ' + str(r2) + '\n', 'RMSE = ' + str(rmse) + '\n']
    write_output(input_list=input_list, output=output + "_report.txt")
    return x_set
