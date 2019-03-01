import pandas as pd
import numpy as np
from scipy.stats import norm
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
from utils import predict_decision, predict_proba


class CigNet:
    """
    CigNet Class
    """
    def __init__(self, predictor, scaler, CellNet, features, CIG_prediction_result, expression):
        """
        :param predictor: predictor of CIGNet, pretrained and save as pkl object
        :param scaler: scaler object, used to normalize real table
        :param CIG_prediction_result: PredictGo results
        :param CellNet: CellNet database table
        """
        self.predictor = predictor
        self.scaler = scaler
        self.CellNet = CellNet
        self.features = features
        self.expression = expression

        if CIG_prediction_result is not None:
            self.CIG_prediction_result = CIG_prediction_result
            self.candidates = self.CIG_prediction_result[self.CIG_prediction_result['FDR'] < 0.05].index.tolist()
            self.net_df = self.GRN_table()
            self.real_table = self.GRN_real_table()
            self.prediction_result = self.CigNet_prediction()

    def load_result(self, CIG_prediction_result):
        self.CIG_prediction_result = CIG_prediction_result
        self.candidates = self.CIG_prediction_result[self.CIG_prediction_result['FDR'] < 0.05].index.tolist()
        self.net_df = self.GRN_table()
        self.real_table = self.GRN_real_table()
        self.prediction_result = self.CigNet_prediction()
        return

    def GRN_table(self):
        return self.CellNet[(self.CellNet['to'].isin(self.candidates)) |
                            (self.CellNet['from'].isin(self.candidates))].copy()

    def GRN_real_table(self):
        results = []
        for gene_id in self.candidates:
            edge_results = set()
            for v in self.net_df[['to', 'from']].values:
                edge_results.add(tuple(v)[0])
                edge_results.add(tuple(v)[1])

            edge_index = edge_results

            edge_stat_df = pd.DataFrame(index=edge_index)
            edge_stat_df['in'] = self.net_df.groupby('to').count()['from']
            edge_stat_df['out'] = self.net_df.groupby('from').count()['to']
            edge_stat_df.fillna(value=0, inplace=True)
            edge_stat_df['total'] = edge_stat_df['in'] + edge_stat_df['out']

            if gene_id not in edge_stat_df.index:
                continue
            parent_edge = edge_stat_df.ix[gene_id, 'in'] + 1
            children_edge = edge_stat_df.ix[gene_id, 'out'] + 1

            children = self.net_df[self.net_df['from'] == gene_id]
            parent = self.net_df[self.net_df['to'] == gene_id]

            children_coef = children['coef'].abs().mean()
            parent_coef = parent['coef'].abs().mean()

            cur_children_genes = children['to'].tolist()

            cur_parent_genes = parent['from'].tolist()

            ## get CIG_distance
            predict_df = self.CIG_prediction_result
            CIG_distance = predict_df.ix[gene_id, 'distance']
            CIG_children_distance = predict_df.ix[cur_children_genes, 'distance'].mean()
            CIG_parent_distance = predict_df.ix[cur_parent_genes, 'distance'].mean()

            if self.expression is not None:
                if gene_id not in self.expression.index:
                    self_exp = 0
                else:
                    self_exp = self.expression.ix[gene_id, 'RNA_exp']

                cur_children_genes = children['to'].tolist()
                cur_children_genes = [g for g in cur_children_genes if g in self.expression.index]
                if len(cur_children_genes) == 0:
                    children_exp = 0
                else:
                    cur_children_exp = self.expression.ix[cur_children_genes, 'RNA_exp']
                    cur_children_exp.fillna(value=0, inplace=True)
                    children_exp = np.mean(cur_children_exp)

                cur_parent_genes = parent['from'].tolist()
                cur_parent_genes = [g for g in cur_parent_genes if g in self.expression.index]
                if len(cur_parent_genes) == 0:
                    parent_exp = 0
                else:
                    cur_parent_exp = self.expression.ix[cur_parent_genes, 'RNA_exp']
                    cur_parent_exp.fillna(value=0, inplace=True)
                    parent_exp = np.mean(cur_parent_exp)
            else:
                self_exp = 0
                children_exp = 0
                parent_exp = 0

            results.append([gene_id, parent_edge, children_edge, parent_coef, children_coef, CIG_distance,
                            CIG_parent_distance, CIG_children_distance, self_exp, children_exp, parent_exp])

        result_df = pd.DataFrame(results)
        result_df.columns = ['gene_id',
                             'parent_edge', 'children_edge',
                             'parent_coef', 'children_coef',
                             'CIG_distance', 'CIG_parent_distance', 'CIG_children_distance',
                             'RNA_exp', 'children_exp', 'parent_exp']
        result_df = result_df.fillna(0)
        result_df = result_df.set_index(['gene_id'])
        if self.expression is not None:
            for gene in result_df.index:
                if gene not in self.expression.index:
                    continue
                result_df.ix[gene, 'RNA_exp'] = self.expression.ix[gene, 'RNA_exp']

        result_df = result_df[self.features]
        return result_df

    def CigNet_prediction(self):
        new_table_data = self.scaler.transform(self.real_table)
        real_table = pd.DataFrame(new_table_data, index=self.real_table.index, columns=self.real_table.columns)

        result = predict_decision(self.predictor, real_table)
        result2 = predict_proba(self.predictor, real_table)
        result = np.concatenate([result, result2], axis=1)
        result_df = pd.DataFrame(result, index=real_table.index)
        result_df.columns = ['distance', 'non-driver_prob', 'driver_prob']

        stats = importr('stats')
        for l in [-2, -1.75, -1.5, -1.25 - 1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 1]:
            result_df['p_value'] = 1 - norm.cdf(result_df['distance'], loc=l)
            result_df['q_value'] = stats.p_adjust(FloatVector(result_df["p_value"].tolist()), method='BH')
            if result_df[result_df['q_value'] < 0.05].shape[0] * 1. / result_df.shape[0] < 0.65:
                break
        candidate_list = self.CellNet['from'].unique()
        result_df = result_df[result_df.index.isin(candidate_list)]
        result_df = result_df.sort_values(by='distance', ascending=False)
        result_df['Rank'] = result_df['distance'].rank(axis=0, ascending=False)
        return result_df
