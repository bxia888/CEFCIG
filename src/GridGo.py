import pandas as pd
from copy import deepcopy
import numpy as np
from multiprocessing import Process, Queue
from utils import logp_wilcoxon, logp_fisher


class GridGo:
    """
    Optimize the genomic features of a set of genes, store current best parameters
    """
    def __init__(self, genes, gtf):
        """
        :param genes: list of (gene object)
        """
        self.genes = genes
        self.gtf = gtf
        self.parameters = {}
        self.gene_ids = [x.gene_id for x in self.genes]

        index = []
        labels = []
        gene_ids = []
        celltypes = []
        expressions = []
        markers = set()
        features = ['total_width', 'total_signal', 'height', 'coverage', 'skewness', 'kurtosis']
        genebody_types = ['TSS', 'genebody']

        for gene in self.genes:
            gene_ids.append(gene.gene_id)
            labels.append(gene.label)
            index.append(gene.gene_id + '_' + gene.celltype)
            celltypes.append(gene.celltype)
            expressions.append(gene.exp)
            for type_signal in gene.signal.keys():
                markers.add(type_signal)
        markers = list(markers)
        columns = ['_'.join([marker, feature, genebody_type]) for marker in markers
                   for feature in features
                   for genebody_type in genebody_types]
        for marker in markers:
            self.parameters[marker] = {}
            for genebody_type in genebody_types:
                self.parameters[marker][genebody_type] = {}
                for feature in features:
                    self.parameters[marker][genebody_type][feature] = (None, None)

        self.table = pd.DataFrame(index=index, columns=columns)
        self.table['label'] = labels
        self.table['gene_id'] = gene_ids
        self.table['cell_type'] = celltypes
        self.table['RNA_exp'] = expressions
        # print self.table

    @staticmethod
    def next_grid(all_range, range_step, cur_step, cur_center, reduction_step, step_limit, search):
        """
        :param all_range:
        :param range_step:
        :param cur_step:
        :param cur_center:
        :param reduction_step:
        :param step_limit:
        :param search: boolean, whether there is next iteration.
        :return:
        """
        new_step = cur_step / reduction_step
        if new_step < step_limit:
            new_step = step_limit
        center_index = all_range.index(cur_center)

        if search:
            if range_step == cur_step:
                start_index = center_index - 2 if center_index - 2 >= 0 else 0
                end_index = center_index + 3 if center_index + 3 <= len(all_range) else len(all_range)
                return all_range[start_index:end_index], new_step
            else:
                start_index = center_index - int(cur_step * 2 / range_step) if center_index - int(
                    cur_step * 2 / range_step) >= 0 else 0
                end_index = center_index + int(cur_step * 2 / range_step) + 1 if center_index + int(
                    cur_step * 2 / range_step) + 1 <= len(all_range) else len(all_range)
                return all_range[start_index:end_index][::int(new_step / range_step)], new_step
        else:
            start_index = center_index - 4 if center_index - 4 >= 0 else 0
            end_index = center_index + 5 if center_index + 5 <= len(all_range) else len(all_range)
            return all_range[start_index:end_index][::int(new_step / range_step)], new_step

    def get_training_table(self):
        train_df = self.table.fillna(0)
        train_df.index.name = 'index'
        train_df = train_df.set_index(['gene_id'], drop=True)
        return train_df

    def get_parameters(self):
        return self.parameters

    def get_stats(self, marker, feature, genebody, up_stream_distance, down_stream_distance, height):
        positive_results = []
        negative_results = []

        for gene in self.genes:
            if genebody:
                cur_start = gene.start + up_stream_distance
                cur_end = gene.end + down_stream_distance
            else:
                cur_start = gene.start + up_stream_distance
                cur_end = gene.start + down_stream_distance
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
            if gene.label == 1:
                positive_results.append(cur_stats)
            else:
                negative_results.append(cur_stats)
        return positive_results, negative_results

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
            self.table.ix[gene.gene_id+'_'+gene.celltype, cur_column] = cur_stats
        return

    def go(self, gtf, marker, feature, cost_function_type, genebody,
           up_stream_distance_range, down_stream_size_range, cutoff_range,
           up_stream_distance_grid=10000, down_stream_size_grid=10000, cutoff_grid=10,
           up_stream_distance_range_step=1000, down_stream_size_range_step=1000, cutoff_range_step=1,
           up_stream_distance_step=2, down_stream_size_step=2, cutoff_step=2,
           up_stream_distance_limit=1000, down_stream_size_limit=1000, cutoff_limit=1,
           process=8, wigs=None, fisher_c=500):
        """

        :param gtf: gene gtf, id need to match the id in meta_df
        :param marker: ChIP-seq marker
        :param feature: total_width, total_signal, kurtosis, skewness, coverage, height
        :param cost_function_type: wilcoxon or fisher exact
        :param genebody: genebody or TSS, boolean
        :param up_stream_distance_range: range of search start site relative to TSS
        :param down_stream_size_range: search range if genebody, then it determines how far to search away from TTS, and
        TTS range is the up_stream_distance_range to TSS
        :param cutoff_range: range of height
        :param up_stream_distance_grid: initial grid size
        :param down_stream_size_grid: initial grid size
        :param cutoff_grid: initial grid size
        :param up_stream_distance_range_step: smallest grid
        :param down_stream_size_range_step: smallest grid
        :param cutoff_range_step: smallest grid
        :param up_stream_distance_step: degree of grid to be shrink
        :param down_stream_size_step: degree of grid to be shrink
        :param cutoff_step: degree of grid to be shrink
        :param up_stream_distance_limit: the grid size limit
        :param down_stream_size_limit: the grid size limit
        :param cutoff_limit: the grid size limit
        :param process: number of process, multiple processing
        :param wigs:
        :param fisher_c: fisher test number cutoff
        :return:
        """

        path = []

        new_up_stream_distance_range = deepcopy(up_stream_distance_range)
        new_down_stream_size_range = deepcopy(down_stream_size_range)
        new_cutoff_range = deepcopy(cutoff_range)

        search = True

        grid_up_stream_distance_range = new_up_stream_distance_range[
                                        up_stream_distance_grid / up_stream_distance_range_step / 2::
                                        (up_stream_distance_grid / up_stream_distance_range_step)]
        grid_down_stream_size_range = new_down_stream_size_range[
                                      down_stream_size_grid / down_stream_size_range_step / 2::
                                      (down_stream_size_grid / down_stream_size_range_step)]
        grid_cutoff_range = new_cutoff_range[int(cutoff_grid / cutoff_range_step / 2)::
                            int(cutoff_grid / cutoff_range_step)]

        past_path = {}

        best_log_P = None
        best_comb = None

        while True:
            combinations = np.stack(np.meshgrid(grid_up_stream_distance_range,
                                                grid_down_stream_size_range,
                                                grid_cutoff_range), -1).reshape(-1, 3)

            if not genebody:
                new_combinations = []
                for comb in combinations:
                    if comb[1] > comb[0]:
                        new_combinations.append(comb)
                combinations = np.asarray(new_combinations)

            chunks = []
            cur_index = 0
            reminder = len(combinations) % process
            chunk_size = len(combinations) / process
            for i in range(process):
                if reminder > 0:
                    chunks.append(combinations[cur_index + i * chunk_size:cur_index + (i + 1) * chunk_size + 1])
                    cur_index += 1
                    reminder -= 1
                else:
                    chunks.append(combinations[cur_index + i * chunk_size: cur_index + (i + 1) * chunk_size])

            total_chunk_size = 0
            for chunk in chunks:
                total_chunk_size += len(chunk)
            if total_chunk_size != len(combinations):
                print 'multiple processes chunk size is not correct'
                return None

            queue = Queue()
            processes = []

            cur_past_path = deepcopy(past_path)
            for i in range(process):
                cur_chunk = chunks[i]
                p = Process(target=self.go_process,
                            args=(queue, cur_chunk, gtf,
                                  feature, cur_past_path, cost_function_type, marker, genebody,
                                  wigs, fisher_c))
                processes.append(p)
                p.start()

            cur_path = []

            for i in range(process):
                cur_best_logP, cur_best_comb, cur_process_path = queue.get()
                cur_path += cur_process_path
                if (best_log_P is None or cur_best_logP < best_log_P) and cur_best_logP is not None:
                    best_log_P = cur_best_logP
                    best_comb = cur_best_comb

            for p in processes:
                p.join()

            for trace in cur_path:
                if tuple(trace[0]) not in past_path:
                    past_path[tuple(trace[0])] = trace[1]

            path += cur_path

            if not search:
                break

            if up_stream_distance_grid == up_stream_distance_limit and \
               down_stream_size_grid == down_stream_size_limit and \
               cutoff_grid == cutoff_limit:
                search = False

            up_stream_center, down_stream_size_center, cutoff_center = best_comb
            grid_up_stream_distance_range, up_stream_distance_grid = self.next_grid(up_stream_distance_range,
                                                                                    up_stream_distance_limit,
                                                                                    up_stream_distance_grid,
                                                                                    up_stream_center,
                                                                                    up_stream_distance_step,
                                                                                    up_stream_distance_limit,
                                                                                    search)
            grid_down_stream_size_range, down_stream_size_grid = self.next_grid(down_stream_size_range,
                                                                                down_stream_size_limit,
                                                                                down_stream_size_grid,
                                                                                down_stream_size_center,
                                                                                down_stream_size_step,
                                                                                down_stream_size_limit,
                                                                                search)
            grid_cutoff_range, cutoff_grid = self.next_grid(cutoff_range,
                                                            cutoff_limit,
                                                            cutoff_grid,
                                                            cutoff_center,
                                                            cutoff_step,
                                                            cutoff_limit,
                                                            search)
        return path, best_comb

    def go_process(self, queue, combinations, gtf,
                        feature, cur_past_path, cost_function_type, marker, genebody,
                        wigs, fisher_c):
        """

        :param queue: queue from multiple process
        :param combinations: grid combinations need to check
        :param gtf: gene gtf
        :param feature: genomic feature
        :param cur_past_path: current grid search path
        :param cost_function_type: type of cost function, wilcoxon or fisher exact test
        :param marker: genomic marker
        :param genebody: boolean, genebody or TSS (True is genebody)
        :param wigs:
        :param fisher_c: fisher exact test
        :return:
        """
        path = []
        best_log_P = None
        best_comb = None

        for comb in combinations:
            cur_up_stream_distance, cur_down_stream_distance, cur_cutoff = comb
            if (cur_up_stream_distance, cur_down_stream_distance, cur_cutoff) in cur_past_path:
                cur_logP = cur_past_path[(cur_up_stream_distance, cur_down_stream_distance, cur_cutoff)]
            else:
                if cost_function_type == "wilcoxon":
                    cur_logP = self.wilcoxon_cost_function(cur_up_stream_distance, cur_down_stream_distance, cur_cutoff,
                                                           feature, marker, genebody, wigs, fisher_c)
                elif cost_function_type == "fisher":
                    cur_logP = self.fisher_cost_function(cur_up_stream_distance, cur_down_stream_distance, cur_cutoff,
                                                         feature, marker, genebody, wigs, fisher_c)
            path.append((comb, cur_logP))
            if (best_log_P is None or best_log_P > cur_logP) and cur_logP is not None:
                best_log_P = cur_logP
                best_comb = comb
        queue.put((best_log_P, best_comb, path))
        return

    def wilcoxon_cost_function(self, cur_up_stream_distance, cur_down_stream_distance, cur_cutoff,
                               feature, marker, genebody, wigs, top_enrich):
        positive_results, negative_results = \
            self.get_stats(marker, feature, genebody,
                           cur_up_stream_distance, cur_down_stream_distance, cur_cutoff)

        cur_logP = logp_wilcoxon(positive_results, negative_results)
        return cur_logP

    def fisher_cost_function(self, cur_up_stream_distance, cur_down_stream_distance, cur_cutoff,
                             feature, marker, genebody, wigs, top_enrich):
        positive_results, negative_results = \
            self.get_stats(marker, feature, genebody,
                           cur_up_stream_distance, cur_down_stream_distance, cur_cutoff)
        cur_logP = logp_fisher(self.gene_ids, positive_results, negative_results, top_enrich=top_enrich)
        return cur_logP

    def grid(self, up_stream_distance_range, down_stream_distance_range, cutoff_range,
             up_stream_distance_grid=10000, down_stream_size_grid=10000, cutoff_grid=10,
             up_stream_distance_range_step=1000, down_stream_size_range_step=1000, cutoff_range_step=1,
             up_stream_distance_step=2, down_stream_size_step=2, cutoff_step=2,
             up_stream_distance_limit=1000, down_stream_size_limit=1000, cutoff_limit=1,
             process=8, wigs=None, fisher_c=500):
        for marker in self.parameters.keys():
            for genebody_type in self.parameters[marker].keys():
                for feature in self.parameters[marker][genebody_type].keys():
                    if genebody_type == 'TSS':
                        genebody = False
                    else:
                        genebody = True
                    path, best_comb = self.go(self.gtf, marker, feature, "wilcoxon", genebody,
                                              up_stream_distance_range, down_stream_distance_range, cutoff_range,
                                              up_stream_distance_grid, down_stream_size_grid, cutoff_grid,
                                              up_stream_distance_range_step, down_stream_size_range_step, cutoff_range_step,
                                              up_stream_distance_step, down_stream_size_step, cutoff_step,
                                              up_stream_distance_limit, down_stream_size_limit, cutoff_limit,
                                              process=process, wigs=wigs, fisher_c=fisher_c)
                    # print marker, feature, genebody_type, best_comb
                    self.parameters[marker][genebody_type][feature] = (path, best_comb)

    def update_training_table(self):
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
