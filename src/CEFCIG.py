#!/usr/bin/env python

import sys
import argparse

from gene import *
from GridGo import *
from utils import *
from wig import Wig
from PredictGo import PredictGo
from CigNet import CigNet


def printHelp():
    print '\nCEFCIG'
    print 'For help information for each function, try:\nCEFCIG <function> -h'
    print '\nFunctions:'
    print '\tGridGo:\n\t\t Optimize the parameter of chipseq peak calling.'
    print '\tPredictGo:\n\t\t Predict cell identity genes or gene of interests based on different signals provided (Chipseq and RNAseq).\n'
    print '\tCigNet:\n\t\t Predict master transcription factors based CellNet database and predictgo prediction result'
    print '\nKaifu Chen, et al. chenkaifu@gmail.com, Chen lab, Cardiovascular department, Houston methodist.'
    print ''


def gridgo():
    if (len(sys.argv) < 7) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\nCEFCIG GridGo <gene meta table> <gene_expression> <gene GTF> <wig meta table> " \
              "\n\nfor more help, please try: CEFCIG grigo -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\nCEFCIG GridGo <gene meta table> <gene GTF> <wig meta table> <gene expression>>\n\n",
                                     description='', epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'GridGo' to perform parameter optimization")

    parser.add_argument('gene_meta_table', default=None,
                        help="gene meta table, containing at least of the positive list of genes of interests.")
    parser.add_argument('gene_expression', default=None,
                        help="gene expression.")
    parser.add_argument('gene_GTF', default=None,
                        help="the file path of gene gtf following the format of UCSC genome browser.")
    parser.add_argument('wig_meta_table', default=None,
                        help="table containing wig file meta data information")
    parser.add_argument('output_dir', default=None,
                        help="output directory")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0

    elif len(sys.argv) >= 7:  # at least two parameter need to be specified
        args = parser.parse_args()  # all paramter values are now saved in args
    else:
        print "\nfor help, please try: CEFCIG GridGo -h\n"
        return 0
    print ''

    if args is not None:
        gene_meta_table_path = args.gene_meta_table
        if gene_meta_table_path.endswith('.xlsx'):
            gene_meta_df = pd.read_excel(gene_meta_table_path)
        elif gene_meta_table_path.endswith('.csv'):
            gene_meta_df = pd.read_csv(gene_meta_table_path)
        elif gene_meta_table_path.endswith('.xls') or gene_meta_table_path.endswith('.txt'):
            gene_meta_df = pd.read_csv(gene_meta_table_path, sep='\t')
        else:
            gene_meta_df = None

        gtf_path = args.gene_GTF
        if gtf_path.endswith('.xlsx'):
            gtf = pd.read_excel(gtf_path)
        elif gtf_path.endswith('.csv'):
            gtf = pd.read_csv(gtf_path)
        elif gtf_path.endswith('.xls') or gtf_path.endswith('.txt'):
            gtf = pd.read_csv(gtf_path, sep='\t')
        else:
            gtf = None

        gtf = gtf.iloc[:, 0:5]
        gtf.columns = ['gene_id', 'chr', 'strand', 'txStart', 'txEnd']

        if len(gene_meta_df['label'].unique()) == 1:
            # only one kind of label found in gene table
            control_gtf = gtf[~gtf['gene_id'].isin(gene_meta_df['gene_id'].unique())]
            control = control_gtf.sample(gene_meta_df.shape[0], random_state=0)
            control_genes = list(control['gene_id'].values)
            control_df = pd.DataFrame(index=range(gene_meta_df.shape[0], 2 * gene_meta_df.shape[0]))
            control_df['gene_id'] = control_genes
            control_df['cell_type'] = list(gene_meta_df['cell_type'])
            control_df['label'] = [0] * gene_meta_df.shape[0]
            gene_meta_df = gene_meta_df.append(control_df)

        wigs_metadata_path = args.wig_meta_table
        if wigs_metadata_path.endswith('.xlsx'):
            wigs_meta_df = pd.read_excel(wigs_metadata_path)
        elif wigs_metadata_path.endswith('.csv'):
            wigs_meta_df = pd.read_csv(wigs_metadata_path)
        elif wigs_metadata_path.endswith('.xls') or wigs_metadata_path.endswith('.txt'):
            wigs_meta_df = pd.read_csv(wigs_metadata_path, sep='\t')
        else:
            wigs_meta_df = None

        gtf = gtf.set_index(['gene_id'])

        gene_objs = []

        for i in range(gene_meta_df.shape[0]):
            gene_id, cell_type, label = gene_meta_df.iloc[i, 0], gene_meta_df.iloc[i, 1], gene_meta_df.iloc[i, 2]
            if gene_id not in gtf.index:
                # print gene_id, label
                continue
            chr, start, end = gtf.ix[gene_id, 'chr'], gtf.ix[gene_id, 'txStart'], gtf.ix[gene_id, 'txEnd']
            if isinstance(chr, pd.Series):
                continue
            gene_obj = Gene(gene_id, cell_type, label, chr, start, end)
            gene_objs.append(gene_obj)

        expression_path = args.gene_expression
        if expression_path.endswith('.xlsx'):
            expression = pd.read_excel(expression_path, index_col=0)
            expression.index = expression.index.str.upper()
        elif expression_path.endswith('.csv'):
            expression = pd.read_csv(expression_path, index_col=0)
            expression.index = expression.index.str.upper()
        elif expression_path.endswith('.xls') or expression_path.endswith('.txt'):
            expression = pd.read_csv(expression_path, sep='\t', index_col=0)
            expression.index = expression.index.str.upper()
        else:
            expression = None

        # print expression
        expression = expression[~expression.index.duplicated(keep='first')]

        for gene_obj in gene_objs:
            if gene_obj.gene_id in expression.index:
                cur_cols = [col for col in expression.columns if col.find(gene_obj.celltype) != -1]
                if len(cur_cols) == 0:
                    print gene_obj.celltype
                else:
                    cur_col = cur_cols[0]
                gene_obj.exp = expression.ix[gene_obj.gene_id, cur_col]
            # else:
            #     print gene_obj.gene_id

        wigs = {}
        for i in range(wigs_meta_df.shape[0]):
            cell_type, marker, wig_path = wigs_meta_df.iloc[i, 0], wigs_meta_df.iloc[i, 1], wigs_meta_df.iloc[i, 2]
            if cell_type not in wigs:
                wigs[cell_type] = {}
            wigs[cell_type][marker] = Wig(wig_path)

        for gene_obj in gene_objs:
            cell_type = gene_obj.celltype
            cur_signals = {}
            start = gene_obj.start - 10000 if gene_obj.start - 10000 > 0 else 0
            end = gene_obj.end + 10000
            chr = gene_obj.chr
            for marker in wigs[cell_type].keys():
                cur_wig_obj = wigs[cell_type][marker]
                cur_signals[marker] = cur_wig_obj.genome[chr].get_signals(start, end)
            gene_obj.add_signal(cur_signals)

        output_dir = args.output_dir
        output_dir = output_dir if output_dir.endswith('/') else output_dir + '/'
        save_obj(gene_objs, output_dir+'gene_objs.pkl')

        gridgo = GridGo(gene_objs, gtf)

        up_stream_distance_range = range(-10000, 0, 1000)
        down_stream_size_range = range(0, 10000, 1000)
        cutoff_range = range(1, 51, 5)

        gridgo.grid(up_stream_distance_range, down_stream_size_range, cutoff_range,
                    up_stream_distance_grid=1000, down_stream_size_grid=1000, cutoff_grid=5,
                    up_stream_distance_range_step=1000, down_stream_size_range_step=1000, cutoff_range_step=1,
                    up_stream_distance_step=2, down_stream_size_step=2, cutoff_step=2,
                    up_stream_distance_limit=1000, down_stream_size_limit=1000, cutoff_limit=1,
                    process=8, wigs=None, fisher_c=500)
        save_obj(gridgo, output_dir+'GridGo_step1.pkl')

        gridgo.update_training_table()
        save_obj(gridgo, output_dir+'GridGo_step2.pkl')

        train_df = gridgo.get_training_table()
        parameter_df = gridgo.get_parameters()
        train_df.to_csv(output_dir+'training_table.xls', sep='\t')
        parameter_df.to_csv(output_dir + 'parameter.xls', sep='\t')


def predictgo():
    if (len(sys.argv) < 8) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\nCEFCIG PredictGo <gene meta table> <gene GTF> <target wig meta table> <gene expression>" \
              "\n\nfor more help, please try: CEFCIG PredictGo -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\nCEFCIG PredictGo <gene meta table> <gridgo_obj> <gene GTF> <target wig meta table> <gene expression>>\n\n",
                                     description='', epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'PredictGo' to perform prediction of cell identity genes or genes of interests")

    parser.add_argument('gene_meta_table', default=None,
                        help="gene meta table, containing at least of the positive list of genes of interests.")
    parser.add_argument('gridgo_obj', default=None,
                        help="gridgo object, gridgo optimized result.")
    parser.add_argument('gene_GTF', default=None,
                        help="the file path of gene gtf following the format of UCSC genome browser.")
    parser.add_argument('target_wig_meta_table', default=None,
                        help="table containing target wig file meta data information")
    parser.add_argument('gene_expression', default=None,
                        help="gene expression")
    parser.add_argument('output_dir', default=None,
                        help="output directory")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0

    elif len(sys.argv) >= 6:  # at least two parameter need to be specified
        args = parser.parse_args()  # all paramter values are now saved in args
    else:
        print "\nfor help, please try: CEFCIG PredictGo -h\n"
        return 0
    print ''

    if args is not None:
    # if True:
        gtf_path = args.gene_GTF

        if gtf_path.endswith('.xlsx'):
            gtf = pd.read_excel(gtf_path)
        elif gtf_path.endswith('.csv'):
            gtf = pd.read_csv(gtf_path)
        elif gtf_path.endswith('.xls') or gtf_path.endswith('.txt'):
            gtf = pd.read_csv(gtf_path, sep='\t')
        else:
            gtf = None

        gtf = gtf.iloc[:, 0:5]
        gtf.columns = ['gene_id', 'chr', 'strand', 'txStart', 'txEnd']

        real_wigs_metadata_path = args.target_wig_meta_table
        if real_wigs_metadata_path.endswith('.xlsx'):
            real_wigs_meta_df = pd.read_excel(real_wigs_metadata_path)
        elif real_wigs_metadata_path.endswith('.csv'):
            real_wigs_meta_df = pd.read_csv(real_wigs_metadata_path)
        elif real_wigs_metadata_path.endswith('.xls') or real_wigs_metadata_path.endswith('.txt'):
            real_wigs_meta_df = pd.read_csv(real_wigs_metadata_path, sep='\t')
        else:
            real_wigs_meta_df = None

        gtf = gtf.set_index(['gene_id'])

        gene_objs = []

        target_celltype = real_wigs_meta_df.iloc[:, 0].unique()[0]

        for gene_id in gtf.index:
            chr, start, end = gtf.ix[gene_id, 'chr'], gtf.ix[gene_id, 'txStart'], gtf.ix[gene_id, 'txEnd']
            if isinstance(chr, pd.Series):
                continue
            if chr == 'chrM':
                continue
            gene_obj = Gene(gene_id, target_celltype, None, chr, start, end)
            gene_objs.append(gene_obj)

        target_wigs = {}
        for i in range(real_wigs_meta_df.shape[0]):
            cell_type, marker, wig_path = real_wigs_meta_df.iloc[i, 0], real_wigs_meta_df.iloc[i, 1], \
                                          real_wigs_meta_df.iloc[i, 2]
            if cell_type not in target_wigs:
                target_wigs[cell_type] = {}
            target_wigs[cell_type][marker] = Wig(wig_path)

        for gene_obj in gene_objs:
            cell_type = gene_obj.celltype
            cur_signals = {}
            start = gene_obj.start - 10000 if gene_obj.start - 10000 > 0 else 0
            end = gene_obj.end + 10000
            chr = gene_obj.chr

            for marker in target_wigs[cell_type].keys():
                cur_wig_obj = target_wigs[cell_type][marker]
                cur_signals[marker] = cur_wig_obj.genome[chr].get_signals(start, end)
            gene_obj.add_signal(cur_signals)

        expression_path = args.gene_expression
        if expression_path.endswith('.xlsx'):
            expression = pd.read_excel(expression_path, index_col=0)
            expression.index = expression.index.str.upper()
        elif expression_path.endswith('.csv'):
            expression = pd.read_csv(expression_path, index_col=0)
            expression.index = expression.index.str.upper()
        elif expression_path.endswith('.xls') or expression_path.endswith('.txt'):
            expression = pd.read_csv(expression_path, sep='\t', index_col=0)
            expression.index = expression.index.str.upper()
        else:
            expression = None

        expression = expression[~expression.index.duplicated(keep='first')]

        expression.index = expression.index.str.upper()
        expression.columns = ['RNA_exp']
        for gene_obj in gene_objs:
            if gene_obj.gene_id in expression.index:
                gene_obj.exp = expression.ix[gene_obj.gene_id, 'RNA_exp']


        gridgo_obj_path = args.gridgo_obj
        gridgo_obj = load_obj(gridgo_obj_path)
        training_table = gridgo_obj.get_training_table()
        parameters = gridgo_obj.parameters

        preknown = []

        metadata_path = args.gene_meta_table
        if metadata_path.endswith('.xlsx'):
            gene_meta_df = pd.read_excel(metadata_path)
        elif metadata_path.endswith('.csv'):
            gene_meta_df = pd.read_csv(metadata_path)
        elif metadata_path.endswith('.xls') or metadata_path.endswith('.txt'):
            gene_meta_df = pd.read_csv(metadata_path, sep='\t')
        else:
            gene_meta_df = None

        gene_meta_df = gene_meta_df.set_index(['gene_id'])
        predict_obj = PredictGo(gene_objs, gtf, training_table, parameters, preknown, gene_meta_df)

        output_dir = args.output_dir
        output_dir = output_dir if output_dir.endswith('/') else output_dir + '/'

        # save_obj(predict_obj, output_dir + '/predictGo.pkl')

        predict_obj.predict()

        final_df = predict_obj.prediction
        final_df.to_csv(output_dir + 'prediction_results.csv')
        # save_obj(predict_obj, output_dir + '/predictGo.pkl')


def cignet():
    if (len(sys.argv) < 5) and ('-h' not in sys.argv) and ('--help' not in sys.argv):
        # at least one parameter need to be specified, will print help message if no parameter is specified
        print "\nusage:\n\nCEFCIG CigNet <predictgo prediction result table> <gene expression> <output directory>" \
              "\n\nfor more help, please try: CEFCIG CigNet -h\n"
        return 1

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     usage="\n\nCEFCIG CigNet <predictgo prediction result table> <gene expression> <output directory>\n\n",
                                     description='', epilog="Chen lab, Houston Methodist")
    parser.add_argument('command', default=None, help="set as 'CigNet' to perform prediction of master transcription factors")
    parser.add_argument('predictgo_prediction_result_table', default=None,
                        help="prediction result table from predictgo")
    parser.add_argument('gene_expression', default=None,
                        help="gene expression table")
    parser.add_argument('output_dir', default=None,
                        help="output directory")

    if '-h' in sys.argv or '--help' in sys.argv:  # print help information once required by user
        print '\n'
        parser.print_help()
        print '\n'
        return 0

    elif len(sys.argv) >= 5:  # at least two parameter need to be specified
        args = parser.parse_args()  # all paramter values are now saved in args
    else:
        print "\nfor help, please try: CEFCIG CigNet -h\n"
        return 0
    print ''

    if args is not None:
        cignet_obj = load_obj('../data/cignet_obj.pkl')
        predictor, scaler, CellNet, features = cignet_obj

        CIG_prediction_result_path = args.predictgo_prediction_result_table
        CIG_prediction_result = pd.read_csv(CIG_prediction_result_path, index_col=0)

        expression_path = args.gene_expression
        if expression_path.endswith('.xlsx'):
            expression = pd.read_excel(expression_path, index_col=0)
            expression.index = expression.index.str.upper()
        elif expression_path.endswith('.csv'):
            expression = pd.read_csv(expression_path, index_col=0)
            expression.index = expression.index.str.upper()
        elif expression_path.endswith('.xls') or expression_path.endswith('.txt'):
            expression = pd.read_csv(expression_path, sep='\t', index_col=0)
            expression.index = expression.index.str.upper()
        else:
            expression = None

        expression.columns = ['RNA_exp']
        expression = expression[~expression.index.duplicated(keep='first')]

        CigNet_obj = CigNet(predictor, scaler, CellNet, features, CIG_prediction_result, expression)
        prediction_result = CigNet_obj.prediction_result

        output_dir = args.output_dir
        output_dir = output_dir if output_dir.endswith('/') else output_dir + '/'
        prediction_result.to_csv(output_dir + 'CigNet_prediction.csv')

        # save_obj(CigNet_obj, output_dir + 'CigNet_obj.pkl')


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == 'GridGo':
            gridgo()
        elif sys.argv[1] == 'PredictGo':
            predictgo()
        elif sys.argv[1] == 'CigNet':
            cignet()
        else:
            printHelp()
    else:
        print 'Computational Epigenetic Framework for Cell Identity Gene Discovery (CEFCIG)\n'
        print 'For a list of functions in CEFCIG, please try:\nCEFCIG -h'
        print 'Kaifu Chen, et al. chenkaifu@gmail.com, Chen lab, Cardiovascular department, Houston methodist.'

