
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

    :param y_col:
    :param training_set:
    :param normalize:
    :param test_size:
    :param folds:
    :param metric:
    :param model_name:
    :param testing_set:
    :param imbalanced:
    :param seed:
    :param include_models:
    :param normalize_method:
    :return:
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
            model=final_model, x_set=testing_set, y_col=testing_set[y_col], output=model_name
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

    :param y_col:
    :param training_set:
    :param normalize:
    :param test_size:
    :param folds:
    :param metric:
    :param model_name:
    :param testing_set:
    :param imbalanced:
    :param seed:
    :param include_models:
    :param normalize_method:
    :return:
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
            model=final_model, x_set=testing_set, y_col=testing_set[y_col], output=model_name
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
    unseen_predictions = pycl.predict_model(model, data=x_set)
    report = metrics.classification_report(y_col, unseen_predictions.Label)
    acc = metrics.accuracy_score(y_col, unseen_predictions.Label)
    auc = metrics.auc(y_col, unseen_predictions.Label)
    plot = generate_confusion_matrix(x_set=x_set, y_set=y_col, output=output)
    input_list = [
        output, '\nTesting model report: \n', report + '\n', 'AUC = ' + str(auc) + '\n', 'Accuracy = ' + str(acc) + '\n'
    ]
    write_output(input_list=input_list, output=output + "_report.txt")
    return unseen_predictions


def test_regressor(
    *,
    model,
    x_set,
    y_col,
    output
):
    x_set['Label'] = model.predict(x_set)
    r2 = metrics.r2_score(y_col, x_set.Label)
    rmse = metrics.mean_squared_error(y_col, x_set.Label, squared=False)
    plot = generate_scatterplot(
        x_axis=x_set.Label, y_axis=y_col, output=output)
    input_list = [output + '\nTesting model report: \n', 'R^2 = ' + str(r2) + '\n', 'RMSE = ' + str(rmse) + '\n']
    write_output(input_list=input_list, output=output + "_report.txt")
    return x_set
