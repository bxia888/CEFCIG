import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(style="white")


def roc_plot(roc_auc, fpr, tpr, n_classes, output_file_path, labels, verbose=False,
             save=False, band=False, plot=False, fig=None, ax=None, color=None, p='', legend=True):
    # plt.figure()
    if not plot:
        if fig is None and ax is None:
            fig, ax = plt.subplots(figsize=(6, 5))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
    color_idx = np.linspace(0, 1, n_classes)
    colors = ['red', 'black', 'cyan', 'pink', 'orange', 'skyblue', 'yellow', 'pink', 'purple']
    # colors = ['black','cyan', 'red', 'black', 'grey', 'green', 'skyblue', 'gray', 'cyan', 'purple',]
    # colors = ['red', 'blue', 'grey', 'grey', 'green', 'skyblue', 'gray', 'cyan', 'purple', ]
    c = 0
    font = {'fontname': 'Helvetica', 'fontsize': 15}
    for i, j in zip(range(n_classes), color_idx):
        # print i, j
        y_info = np.mean(tpr[labels[i]], axis=0)
        y_info[0] = 0
        if color is None:
            plt.plot(fpr[0:], y_info, color=colors[c], linestyle='-', alpha=1,
                     label='{0} (area = {1:0.2f})'
                           ''.format(labels[i], roc_auc[labels[i]]))
            c += 1
        else:
            plt.plot(fpr[0:], y_info[0:], color=color, linestyle='--',alpha=1,
                     label='{0} (area = {1:0.2f})'
                           ''.format(labels[i], roc_auc[labels[i]]))
            # label = '{0} (area = {1:0.2f})'
            # ''.format(labels[i], roc_auc[labels[i]])

        if band:
            mean_tpr = np.mean(tpr[labels[i]], axis=0)
            mean_tpr[-1] = 1.0
            std_tpr = np.std(tpr[labels[i]], axis=0)
            tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
            tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
            plt.fill_between(fpr, tprs_lower, tprs_upper, color=plt.cm.plasma(j), alpha=.3)

    if p != '':
        ax.text(0.95, 0.5, 'p value: '+p, verticalalignment='bottom', horizontalalignment='right',
                transform=ax.transAxes, **font)

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate', **font)
    plt.ylabel('True Positive Rate', **font)
    plt.title(output_file_path)
    if legend:
        plt.legend(loc="lower right")
    if verbose:
        plt.show()
    # plt.axis('off')
    if save:
        if output_file_path == '' or output_file_path is None:
            plt.savefig('result' + '.pdf', transparent=True)
        else:
            plt.savefig(output_file_path+'.pdf', transparent=True)
        plt.close('all')
    return roc_auc, fpr, tpr
