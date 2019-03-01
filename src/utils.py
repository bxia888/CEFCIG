# from scipy import stats
import numpy as np
import pickle
import pandas as pd
import traceback
from sklearn.linear_model import LogisticRegression
from scipy import stats
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn import preprocessing
# from scipy.stats import norm, lognorm
from sklearn.metrics import roc_curve, auc


def peak_skewness(peak_signals):
    """
    :param peak_signals: has two rows, first rows is genome index, 2nd is the signals
    :return:
    """
    values = []
    for i in range(peak_signals.shape[1]):
        index = peak_signals[0, i]
        count = int(peak_signals[1, i] * 10)
        values += [index] * count
    return stats.skew(np.asarray(values))


def peak_kurtosis(peak_signals):
    """
    :param peak_signals: has two rows, first rows is genome index, 2nd is the signals
    :return:
    """
    values = []
    for i in range(peak_signals.shape[1]):
        index = peak_signals[0, i]
        count = int(peak_signals[1, i] * 10)
        values += [index] * count
    return stats.kurtosis(np.asarray(values))


def save_obj(obj, name):
    with open(name, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)


def load_obj(name):
    f = open(name, 'rb')
    result = pickle.load(f)
    f.close()
    return result


def genome_size_chrom(path="./hg19_chr_sizes.txt"):
    """
    :param path: the location where is the genome size file
    :return: a dictionary in which key is the chromosome name (str), value is the chromosome size (int)
    """
    genome = {}
    genome_size_file = open(path, "r")
    for line in genome_size_file.readlines():
        chr_name, chr_size = line.split()
        genome[chr_name] = int(chr_size)
    genome_size_file.close()
    return genome


def label_label(table, positive_table):
    result = [0] * table.shape[0]
    table['label'] = result
    table.loc[table.index.isin(positive_table.index), 'label'] = 1
    return table


def predict_proba(predictor, real_table, label=False):
    if label:
        table = real_table.ix[:, :-1].copy()
    else:
        table = real_table.copy()
    result = predictor.predict_proba(table)
    return pd.DataFrame(result, index=table.index)


def predict_decision(predictor, real_table, label=False):
    if label:
        table = real_table.ix[:, :-1].copy()
    else:
        table = real_table.copy()
    result = predictor.decision_function(table)
    return pd.DataFrame(result, index=table.index)


def center_normalization(table):
    return preprocessing.StandardScaler().fit(table)


def preprocessing_table(scaler, table):
    new_table_data = scaler.transform(table)
    new_table = pd.DataFrame(new_table_data, index=table.index, columns=table.columns)
    return new_table


def score(x, y, predictor):
    return predictor.score(x, y)


def get_training_testing(table, test_size=.2, random_state=0, n_class=2):
    y = table[table.columns[-1]].tolist()
    y = label_binarize(y, classes=range(n_class))
    x = table.iloc[:, :-1]
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=test_size, random_state=random_state)
    return x_train, x_test, y_train, y_test


def predict_logisticregression(x_train, y_train, c=1., penalty='l1'):
    """
    :param x_train: training table, value part
    :param y_train: training table, label part
    :param c: regularization
    :param penalty: l1 or l2 regularizaion
    :return: the sklearn logistic regression object
    """
    predictor = LogisticRegression(penalty=penalty, C=c).fit(x_train, y_train)
    return predictor


def roc(y_true, y_score, labels):
    """
    Compute ROC curve and ROC area for each class
    :param y_true:
    :param y_score:
    :param labels:
    :return:
    """
    try:
        n_classes = y_true.shape[1]
    except:
        n_classes = 1
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(n_classes):
        fpr[labels[i]], tpr[labels[i]], _ = roc_curve(y_true[:, i], y_score[:, i])

        roc_auc[labels[i]] = auc(fpr[labels[i]], tpr[labels[i]])
    return roc_auc, fpr, tpr


def logp_wilcoxon(groupA, groupB, bags=1):
    """
    return the negative log P value for two groups
    :param groupA:
    :param groupB:
    :return:
    """
    p_values = []

    alloc = 0.75

    if bags == 1:
        try:
            rank_diff, p1 = stats.mannwhitneyu(groupA, groupB, alternative='less')
            rank_diff, p2 = stats.mannwhitneyu(groupB, groupA, alternative='less')
            if p1 < p2:
                return np.log10(p1)
            else:
                return np.log10(p2)
        except:
            print groupA, groupB
            traceback.print_exc()
            return np.log10(1)
    else:
        for i in range(bags):
            cur_groupA = np.random.choice(groupA, int(len(groupA) * alloc), replace=True)
            cur_groupB = np.random.choice(groupB, int(len(groupB) * alloc), replace=True)
            try:
                rank_diff, p1 = stats.mannwhitneyu(cur_groupA, cur_groupB, alternative='less')
                rank_diff, p2 = stats.mannwhitneyu(cur_groupB, cur_groupA, alternative='less')
                if p1 < p2:
                    p = p1
                else:
                    p = p2
                p_values.append(p)
            except:
                p_values.append(1)
        final_p = np.mean(p_values)
        return np.log10(final_p)


def logp_fisher(gene_df, positive_results, negative_results, top_enrich=500, ascending=False):
    total_genes = len(positive_results) + len(negative_results)
    target_genes = gene_df.shape[0]
    sort_result = sorted(list(positive_results)+list(negative_results))
    cutoff_value = sort_result[top_enrich]
    positive_results = np.asarray(positive_results)
    overlap = len(positive_results[positive_results>cutoff_value])
    not_overlap = top_enrich - overlap
    # print overlap, top_enrich, not_overlap, total_genes
    p = stats.fisher_exact([[overlap, target_genes - overlap], [not_overlap, total_genes - not_overlap - target_genes]])
    return np.log10(p)
