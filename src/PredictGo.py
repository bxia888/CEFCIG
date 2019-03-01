import pandas as pd
from copy import deepcopy
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
# matplotlib.use('Agg')
from collections import defaultdict
from scipy import interp
# from scipy.stats import cumfreq
from scipy.stats import norm
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from utils import get_training_testing, label_label, predict_decision, predict_proba, score, \
    preprocessing_table, center_normalization, predict_logisticregression, roc
from plot import roc_plot
stats = importr('stats')


def preknown_cost(meta_df, train_df, real_table, columns, preknown_list, verbose=False):
    cur_train_df = train_df[columns].copy()
    cur_train_df = label_label(cur_train_df, meta_df)
    scaler = center_normalization(cur_train_df.iloc[:, :-1])
    cur_train_df.iloc[:, :-1] = preprocessing_table(scaler, cur_train_df.iloc[:, :-1])
    predictor = predict_logisticregression(cur_train_df.iloc[:, :-1], cur_train_df.iloc[:, -1], c=.2)

    for i in range(len(predictor.coef_[0])):
        if columns[i].find('h3k27me3') != -1 and columns[i].find('kurtosis') == -1:
            predictor.coef_[0][i] = -1
        elif columns[i].find('h3k27me3') != -1 and columns[i].find('kurtosis') != -1:
            predictor.coef_[0][i] = 1

    cur_real_table = real_table[columns].copy()

    cur_real_table = preprocessing_table(scaler, cur_real_table)

    result = predict_decision(predictor, cur_real_table)
    result2 = predict_proba(predictor, cur_real_table)
    result = np.concatenate([result, result2], axis=1)
    result_df = pd.DataFrame(result, index=cur_real_table.index)
    result_df.columns = ['distance', 'non-CIG_prob', 'CIG_prob']
    del result_df['non-CIG_prob']
    result_df = result_df.sort_values(by=['distance'], ascending=False)
    result_df['rank'] = range(1, result_df.shape[0] + 1)

    cur_hit_df = result_df[(result_df.index.isin(preknown_list))]
    if verbose:
        print cur_hit_df

    cur_hit = cur_hit_df[cur_hit_df['rank'] <= 10].shape[0] * 15
    cur_hit += cur_hit_df[(cur_hit_df['rank'] > 10) & (cur_hit_df['rank'] < 20)].shape[0] * 10
    cur_hit += cur_hit_df[(cur_hit_df['rank'] > 20) & (cur_hit_df['rank'] < 50)].shape[0] * 8
    cur_hit += cur_hit_df[(cur_hit_df['rank'] > 50) & (cur_hit_df['rank'] < 100)].shape[0] * 3
    cur_hit += cur_hit_df[(cur_hit_df['rank'] > 100) & (cur_hit_df['rank'] < 200)].shape[0] * 1
    cur_hit += cur_hit_df[(cur_hit_df['rank'] > 200) & (cur_hit_df['rank'] < 500)].shape[0] * 0.1

    if verbose:
        print cur_hit, columns
    return cur_hit


def optimization(meta_df, train_df, real_table, active_columns, suppress_columns, preknown_list):
    pool = active_columns
    pool2 = suppress_columns
    existing_pool = []
    paths = {}
    while len(pool) > 0:
        if len(existing_pool) == 0:
            pass
        else:
            cur_hit = preknown_cost(meta_df, train_df.copy(), real_table.copy(), existing_pool, preknown_list)
            paths[len(existing_pool)] = (cur_hit, deepcopy(existing_pool))

        if len(pool) == 1:
            existing_pool.append(pool[0])
            cur_hit = preknown_cost(meta_df, train_df.copy(), real_table.copy(), existing_pool, preknown_list)
            paths[len(existing_pool)] = (cur_hit, deepcopy(existing_pool))
            break
        to_be_add = None
        best = None
        for i in range(len(pool)):
            new_pool = existing_pool + [pool[i]]
            cur_hit = preknown_cost(meta_df, train_df.copy(), real_table.copy(), new_pool, preknown_list)
            if best is None:
                to_be_add = i
                best = cur_hit
            elif best < cur_hit:
                to_be_add = i
                best = cur_hit
            else:
                pass
        if to_be_add is None:
            to_be_add = 0
        existing_pool.append(pool[to_be_add])
        pool = pool[:to_be_add] + pool[to_be_add + 1:]

    while len(pool2) > 0:
        if len(existing_pool) == 0:
            pass
        else:
            cur_hit = preknown_cost(meta_df, train_df.copy(), real_table.copy(), existing_pool, preknown_list)
            paths[len(existing_pool)] = (cur_hit, deepcopy(existing_pool))
        best = None
        if len(pool2) == 1:
            existing_pool.append(pool2[0])
            cur_hit = preknown_cost(meta_df, train_df.copy(), real_table.copy(), existing_pool, preknown_list)
            paths[len(existing_pool)] = (cur_hit, deepcopy(existing_pool))
            break
        to_be_add = None
        for i in range(len(pool2)):
            new_pool = existing_pool + [pool2[i]]
            cur_hit = preknown_cost(meta_df, train_df.copy(), real_table.copy(), new_pool, preknown_list)
            if best is None:
                to_be_add = i
                best = cur_hit
            elif best < cur_hit:
                to_be_add = i
                best = cur_hit
            else:
                pass
        if to_be_add is None:
            to_be_add = 0
        existing_pool.append(pool2[to_be_add])
        pool2 = pool2[:to_be_add] + pool2[to_be_add + 1:]
    return paths



class PredictGo:
    """
    Optimize the genomic features of a set of genes, store current best parameters
    """
    def __init__(self, genes, gtf, training_table, parameters, preknown, gene_meta_df):
        """

        :param genes: gene objects, all genes in gtf
        :param gtf: gtf
        :param training_table: training table for the predictor
        :param parameters: parameters from GridGo
        :param preknown: preknown knowledge, empty list if nothing known or None
        :param gene_meta_df: curated genes meta information
        """
        self.genes = genes
        self.gtf = gtf
        self.training_table = training_table
        self.parameters = parameters
        self.preknown = preknown
        self.gene_meta_df = gene_meta_df
        self.prediction = None
        self.gene_ids = [x.gene_id for x in self.genes]


        index = []
        labels = []
        gene_ids = []
        expressions = []
        markers = set()
        features = ['total_width', 'total_signal', 'height', 'coverage', 'skewness', 'kurtosis']
        genebody_types = ['TSS', 'genebody']

        for gene in self.genes:
            gene_ids.append(gene.gene_id)
            labels.append(gene.label)
            index.append(gene.gene_id + '_' + gene.celltype)
            expressions.append(gene.exp)
            for type_signal in gene.signal.keys():
                markers.add(type_signal)
        markers = list(markers)
        columns = ['_'.join([marker, feature, genebody_type]) for marker in markers
                   for feature in features
                   for genebody_type in genebody_types]

        self.real_table = pd.DataFrame(index=index, columns=columns)
        self.real_table.index.name = 'gene_id'
        self.real_table['label'] = labels
        self.real_table['gene_id'] = gene_ids
        self.real_table['RNA_exp'] = expressions

        self.update_real_table()
        self.real_table = self.real_table.fillna(0)

    def get_prediction_result(self):
        return self.prediction

    def get_real_table(self):
        return self.real_table

    def get_parameters(self):
        return self.parameters

    def update_stats(self, marker, feature, genebody, up_stream_distance, down_stream_distance, height):
        for gene in self.genes:
            if genebody:
                cur_start = gene.start + up_stream_distance
                cur_end = gene.end + down_stream_distance
                genebody_type = 'genebody'
            else:
                cur_start = gene.start + up_stream_distance
                cur_end = gene.start + down_stream_distance
                genebody_type = 'TSS'
            if feature == "total_signal":
                cur_stats = gene.get_total_signal(marker, cur_start, cur_end, height, step=10)
            elif feature == "total_width":
                cur_stats = gene.get_total_width(marker, cur_start, cur_end, height, step=10)
            elif feature == "height":
                cur_stats = gene.get_height(marker, cur_start, cur_end, height, step=10)
            elif feature == "kurtosis":
                cur_stats = gene.get_kurtosis(marker, cur_start, cur_end, height, step=10)
            elif feature == "skewness":
                cur_stats = gene.get_skewness(marker, cur_start, cur_end, height, step=10)
            elif feature == "coverage":
                cur_stats = gene.get_coverage(marker, cur_start, cur_end, height, step=10)
            else:
                cur_stats = None
            cur_column = '_'.join([marker, feature, genebody_type])
            self.real_table.ix[gene.gene_id+'_'+gene.celltype, cur_column] = cur_stats
        return

    def update_real_table(self):
        for marker in self.parameters.keys():
            for genebody_type in self.parameters[marker].keys():
                for feature in self.parameters[marker][genebody_type].keys():
                    if genebody_type == 'TSS':
                        genebody = False
                    else:
                        genebody = True
                    _, best_comb = self.parameters[marker][genebody_type][feature]
                    up_stream_distance, down_stream_distance, height = best_comb
                    self.update_stats(marker, feature, genebody, up_stream_distance, down_stream_distance, height)
        return

    def predict(self):
        active_columns = [x for x in self.real_table.columns if (x.lower().find('h3k4me3') != -1 or
                                                                 x.lower().find('h3k27ac') != -1 or
                                                                 x.lower().find('h3k4me1') != -1 or
                                                                 x.lower().find('h3k79me2') != -1 or
                                                                 x.lower().find('rna') != -1)
                                                                and x.find('single') == -1]

        suppressive_columns = [x for x in self.real_table.columns if (x.lower().find('h3k27me3') != -1 or
                                                                      x.lower().find('h3k9me3') != -1)
                                                                and x.find('single') == -1]

        paths = optimization(self.gene_meta_df, self.training_table, self.real_table,
                             active_columns, suppressive_columns, self.preknown)

        best_path = None
        best_score = None
        for key, value in paths.items():
            cur_score, cur_path = value
            if cur_score > best_score:
                best_score = cur_score
                best_path = cur_path
        train_df = self.training_table[best_path].copy()
        train_df = label_label(train_df, self.gene_meta_df)
        scaler = center_normalization(train_df.iloc[:, :-1])
        train_df.iloc[:, :-1] = preprocessing_table(scaler, train_df.iloc[:, :-1])
        predictor = predict_logisticregression(train_df.iloc[:, :-1], train_df.iloc[:, -1], c=.2)

        real_table = self.real_table[best_path].copy()
        real_table = preprocessing_table(scaler, real_table)

        result = predict_decision(predictor, real_table)
        result2 = predict_proba(predictor, real_table)

        result = np.concatenate([result, result2], axis=1)
        df = pd.DataFrame(result, index=real_table.index)
        df.columns = ['distance', 'non-CIG_prob', 'CIG_prob']
        del df['non-CIG_prob']

        for l in np.arange(-20, 1, 0.25):
            df['p_value'] = 1 - norm.cdf(df['distance'], loc=l)
            df['FDR'] = stats.p_adjust(FloatVector(df["p_value"].tolist()),
                                       method='BH')
            if df[df['FDR'] < 0.05].shape[0] < 600:
                break

        df = df.sort_values(by=['distance'], ascending=False)
        df['rank'] = range(1, df.shape[0] + 1)
        self.prediction = df


    def prediction_roc(self, groups, output_file_path, labels, test_size=0.2, c=.2, verbose=False, save=True,
                       iterations=100, plot=True, kind='test', method='AUC',
                       band=False, fig=None, ax=None, color=None,
                       test_function=get_training_testing, p='', penalty='l1'):
        """

        :param groups: list containing groups of ROC features, for example [[H3K4me3 columns], [H3K27ac columns]]
        :param output_file_path: output file path and prefix
        :param labels: labels of each group
        :param test_size:
        :param c:
        :param verbose:
        :param save:
        :param iterations:
        :param plot:
        :param kind:
        :param method:
        :param band:
        :param fig:
        :param ax:
        :param color:
        :param test_function:
        :param p:
        :param penalty:
        :return:
        """
        train_df = label_label(self.training_table, self.gene_meta_df)

        all_auc_train = {}
        all_tpr_train = {}

        all_auc_test = {}
        all_tpr_test = {}

        scores = defaultdict(float)

        mean_fpr = np.linspace(0, 1, 101)

        for r in range(iterations):
            x_train, x_test, y_train, y_test = test_function(train_df, random_state=r, test_size=test_size)

            y_trues_train = None
            y_scores_train = None

            y_trues_test = None
            y_scores_test = None

            first = True

            for i in range(len(groups)):
                cur_column = groups[i]
                cur_x_train = x_train.ix[:, cur_column]
                cur_x_test = x_test.ix[:, cur_column]

                if len(cur_x_train.shape) == 1:
                    cur_x_train = cur_x_train.to_frame()
                    cur_x_test = cur_x_test.to_frame()

                cur_predictor = predict_logisticregression(cur_x_train, y_train, penalty=penalty, c=c)

                if method == 'score':
                    cur_score = score(cur_x_test, y_test, cur_predictor)
                    scores[labels[i]] += cur_score

                cur_y_score_train = predict_decision(cur_predictor, cur_x_train, False)
                cur_y_score_test = predict_decision(cur_predictor, cur_x_test, False)
                if first:
                    y_trues_train = y_train
                    y_trues_test = y_test
                    y_scores_train = cur_y_score_train.values
                    y_scores_test = cur_y_score_test.values
                    first = False
                else:
                    y_trues_train = np.concatenate((y_trues_train, y_train), axis=1)
                    y_trues_test = np.concatenate((y_trues_test, y_test), axis=1)
                    y_scores_train = np.concatenate((y_scores_train, cur_y_score_train.values), axis=1)
                    y_scores_test = np.concatenate((y_scores_test, cur_y_score_test.values), axis=1)

            auc_train, fpr_train, tpr_train = roc(y_trues_train, y_scores_train, labels)
            auc_test, fpr_test, tpr_test = roc(y_trues_test, y_scores_test, labels)

            for label in labels:
                if kind == 'train' and band:
                    plt.plot(fpr_train[label], tpr_train[label], lw=0.4, alpha=0.1, color='grey')
                elif kind == 'test' and band:
                    plt.plot(fpr_test[label], tpr_test[label], lw=0.4, alpha=0.1, color='grey')
                tpr_train[label] = interp(mean_fpr, fpr_train[label], tpr_train[label])
                tpr_test[label] = interp(mean_fpr, fpr_test[label], tpr_test[label])

            for l in range(len(labels)):
                label = labels[l]
                if label not in all_auc_train:
                    all_auc_train[label] = auc_train[label]
                    all_auc_test[label] = auc_test[label]
                    all_tpr_train[label] = [tpr_train[label]]
                    all_tpr_test[label] = [tpr_test[label]]
                else:
                    all_auc_train[label] += auc_train[label]
                    all_auc_test[label] += auc_test[label]
                    all_tpr_train[label].append(tpr_train[label])
                    all_tpr_test[label].append(tpr_test[label])
        for label in labels:
            all_auc_train[label] /= iterations
            all_auc_test[label] /= iterations

            if method == 'score':
                scores[label] /= iterations

        if plot:
            if kind == 'train':
                roc_plot(all_auc_train, mean_fpr, all_tpr_train, len(labels), output_file_path,
                         labels, verbose=verbose, save=save, band=band, fig=fig, ax=ax, color=color, p=p)
            elif kind == 'test':
                roc_plot(all_auc_test, mean_fpr, all_tpr_test, len(labels), output_file_path,
                         labels, verbose=verbose, save=save, band=band, fig=fig, ax=ax, color=color, p=p)

        if method == 'AUC':
            results_train_df = pd.DataFrame.from_dict(all_auc_train, orient='index')
            results_train_df.columns = ['train']
            results_test_df = pd.DataFrame.from_dict(all_auc_test, orient='index')
            results_test_df.columns = ['test']

            result_df = results_train_df.join(results_test_df)
            return result_df
        elif method == 'score':
            result_df = pd.DataFrame.from_dict(scores, orient='index')
            result_df.columns = ['accuracy']
        else:
            result_df = None

        if kind == 'test':
            for key in all_tpr_test.keys():
                all_tpr_test[key] = np.mean(all_tpr_test[key], axis=0)
                # TPRs_test_df = pd.DataFrame.from_dict(all_tpr_test)
                # TPRs_test_df.to_csv(title+'TPR.csv')
        return result_df
