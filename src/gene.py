import numpy as np
from utils import peak_skewness, peak_kurtosis


class Gene:
    """
    Gene object, save a gene's information
    """
    def __init__(self, gene_id, celltype, label, chr, start, end, step=10, signal=None, exp=None, cur_signal=None):
        """
        :param gene_id: unique gene identifier
        :param celltype: biosample term name
        :param label: gene label, positive or control, 1 or 0
        :param signal: dict, save the tracks of genomic markers
        :param exp: gene expression
        :param cur_signal: dict, save the tracks of genomic markers based on the current parameters
        """
        self.gene_id = gene_id
        self.celltype = celltype
        self.label = label
        self.chr = chr
        self.start = start
        self.end = end
        self.step = step
        if signal is not None:
            self.signal = signal
        else:
            self.signal = {}
        self.exp = exp
        if cur_signal is not None:
            self.cur_signal = cur_signal
        else:
            self.cur_signal = {}

    @classmethod
    def load_signal(cls, gene_id, celltype, label, chr, start, end, **kwargs):
        """
        :param gene_id: unique gene identifier
        :param celltype: biosample term name
        :param label: positive or control
        :param kwargs: key: genomic marker type, eg. H3K4me3.., value: numpy array of genomic tracks,
        first row is genome index, 2nd row is genome value, step need to 10.
        :return: gene object
        """
        gene_obj = cls(gene_id, celltype, label, chr, start, end)
        for key, value in kwargs.items():
            gene_obj.signal[key] = value
        return gene_obj

    @classmethod
    def load_exp(cls, gene_id, celltype, label, chr, start, end, exp):
        """
        :param gene_id: unique gene identifier
        :param celltype: biosample term name
        :param label: positive or control
        :param exp: RNA expression
        :return: gene object
        """
        gene_obj = cls(gene_id, celltype, label, chr, start, end)
        gene_obj.exp = exp
        return gene_obj

    def add_exp(self, exp):
        self.exp = exp

    def add_signal(self, signals):
        for key, value in signals.items():
            self.signal[key] = value

    def get_exp(self):
        return self.exp

    def get_signal(self, marker, start, end, step):
        """
        :param marker: type of genomic marker, eg. H3K4me3
        :param start: start pos, genome index
        :param end: end pos, genome index
        :param step: genome step
        :return: total signal
        """
        if end < start:
            print "start position need to be smaller than end position"
            return None

        marker_tracks = self.signal[marker]
        if end % step == 0:
            start = start / step
            end = end / step
        else:
            start = start / step
            end = end / step + 1

        start_index = np.searchsorted(marker_tracks[0, :], start)
        end_index = np.searchsorted(marker_tracks[0, :], end)

        return marker_tracks[:, start_index:end_index]

    def update_cur_signal(self, marker, start, end, step):
        cur_signal = self.get_signal(marker, start, end, step)
        # print cur_signal
        self.cur_signal[marker] = cur_signal
        return

    def get_total_signal(self, marker, start, end, height, step):
        self.update_cur_signal(marker, start, end, step)
        cur_values = self.cur_signal[marker]
        if cur_values is None or len(cur_values) == 0 or len(cur_values[1, :]) == 0:
            total_signal = 0
        else:
            total_signal = cur_values[1, cur_values[1, :] > height].sum() * step
        return total_signal

    def get_total_width(self, marker, start, end, height, step):
        self.update_cur_signal(marker, start, end, step)
        cur_values = self.cur_signal[marker]
        if cur_values is None or len(cur_values) == 0 or len(cur_values[1, :]) == 0:
            total_width = 0
        else:
            total_width = cur_values[1, cur_values[1, :] > height].shape[0] * step
        return total_width

    def get_height(self, marker, start, end, height, step):
        self.update_cur_signal(marker, start, end, step)
        cur_values = self.cur_signal[marker]
        if cur_values is None or len(cur_values) == 0 or len(cur_values[1, :]) == 0:
            max_height = 0
        else:
            max_height = cur_values[1, :].max()
        if max_height >= height:
            return max_height
        else:
            return 0

    def get_kurtosis(self, marker, start, end, height, step):
        self.update_cur_signal(marker, start, end, step)
        cur_values = self.cur_signal[marker]
        if cur_values is None or len(cur_values) == 0 or len(cur_values[1, :]) == 0:
            return 0
        else:
            return peak_kurtosis(cur_values[:, cur_values[1, :] > height])

    def get_skewness(self, marker, start, end, height, step):
        self.update_cur_signal(marker, start, end, step)
        cur_values = self.cur_signal[marker]
        if cur_values is None or len(cur_values) == 0 or len(cur_values[1, :]) == 0:
            return 0
        else:
            return peak_skewness(cur_values[:, cur_values[1, :] > height])

    def get_coverage(self, marker, start, end, height, step):
        self.update_cur_signal(marker, start, end, step)
        cur_values = self.cur_signal[marker]
        if cur_values is None or len(cur_values) == 0 or len(cur_values[1, :]) == 0:
            return 0
        else:
            total_width = cur_values[1, cur_values[1, :] > height].shape[0] * step
            coverage = total_width*1./(end-start)
            return coverage
