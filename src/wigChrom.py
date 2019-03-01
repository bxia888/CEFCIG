import numpy as np
from utils import peak_skewness, peak_kurtosis


class WigChrom:
    """
    This is a class for storing the chromosome information from a wig file
    """
    def __init__(self, chr_name, start, size, step, span, fixed=True):
        self.chr_name = chr_name
        self.start = start
        self.size = size
        self.step = step
        self.span = span
        self.fixed = fixed

        vector_size = self.size / step if self.size % step == 0 else self.size/step + 1
        self.signals = np.zeros(vector_size)
        self.hasReload = False

    def reload(self):
        if not self.hasReload:
            values = self.signals[self.signals > 0]
            indexes = np.where(self.signals > 0)[0]

            signals = [indexes, values]

            del self.signals
            self.hasReload = True
            self.signals = np.asarray(signals)
        else:
            return

    def get_signals(self, start, end):
        if start < 0:
            print "start position need to be bigger than 0"
            return
        # if end > self.size:
        #     print "end position need to be smaller than genome size"
        #     return
        if end % self.step == 0:
            start = start/self.step
            end = end/self.step
        else:
            start = start / self.step
            end = end / self.step + 1

        start_index = np.searchsorted(self.signals[0, :], start)
        end_index = np.searchsorted(self.signals[0, :], end)

        return self.signals[:, start_index:end_index]

    def get_peaks(self, start, end, cutoff, min_width):
        """
        call peaks
        :param start:
        :param end:
        :return:
        """
        if start is None and end is None:
            signals = self.signals.copy()
        else:
            signals = self.get_signals(start, end)
        candidates = np.where(signals[1, :] >= cutoff)[0]
        signals = signals[:, candidates]
        indexes = signals[0, :]

        peaks = []
        start = 0
        for i in range(len(indexes) - 1):
            if indexes[i + 1] - indexes[i] == 1:
                continue
            else:
                peak_signal = signals[:, start:i+1]
                peak_start = int(signals[0, start] * self.step)
                peak_end = int((signals[0, i] + 1) * self.step)
                peak_height = np.max(peak_signal[1, :])
                peak_width = peak_end - peak_start
                peak_total_signal = np.sum(peak_signal[1, :]) * self.step
                cur_peak_skewness, cur_peak_kutosis = peak_skewness(peak_signal), peak_kurtosis(peak_signal)

                start = i + 1

                if min_width is None or peak_width >= min_width:
                    peaks.append([self.chr_name, peak_start, peak_end, (peak_start+peak_end)/2, peak_width,
                                  peak_total_signal, peak_height, 0,
                                  cur_peak_skewness, cur_peak_kutosis])

        if start < signals.shape[1]:
            peak_signal = signals[:, start:len(indexes)]
            peak_start = int(signals[0, start] * self.step)
            peak_end = int((signals[0, -1] + 1) * self.step)
            peak_height = np.max(peak_signal[1, :])
            peak_width = peak_end - peak_start
            peak_total_signal = np.sum(peak_signal[1, :]) * self.step
            cur_peak_skewness, cur_peak_kutosis = peak_skewness(peak_signal), peak_kurtosis(peak_signal)
            peaks.append([self.chr_name, peak_start, peak_end, (peak_start+peak_end)/2, peak_width,
                          peak_total_signal, peak_height, 0,
                          cur_peak_skewness, cur_peak_kutosis])
        return peaks

    def get_sk(self, start, end):
        """
        calculate the skewness and kurtosis, based on reference peak information
        :param start:
        :param end:
        :return:
        """
        signals = self.get_signals(start, end)
        skewness_value = peak_skewness(signals)
        kurtosis_value = peak_kurtosis(signals)
        return skewness_value, kurtosis_value
