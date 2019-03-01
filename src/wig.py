import pandas as pd
import wigChrom
from utils import genome_size_chrom


class Wig:
    """
    This is the class for input a list Wig files and generate the partitions for the wig files
    """
    def __init__(self, wig_file, genome_size_path="../ref_data/hg19_chr_sizes.txt"):
        self.genome_size = genome_size_chrom(genome_size_path)
        self.genome = {}
        self.initiate(wig_file)
        self.splitted_chroms = {}
        wig_file_name = wig_file.split("/")[-1]
        self.file_name = wig_file_name[:-4]

    def initiate(self, wig_file):
        wig = open(wig_file, "r")
        chr_name = None
        start = None
        step = None
        span = None

        cur_position = 0

        for line in wig.readlines():
            cur_line = line.rstrip().split()
            if len(cur_line) > 1:
                if cur_line[1].startswith("chrom="):
                    chr_name = cur_line[1][cur_line[1].find("=")+1:]
                if cur_line[2].startswith("start="):
                    start = int(cur_line[2][cur_line[2].find("=")+1:])
                if cur_line[3].startswith("step="):
                    step = int(cur_line[3][cur_line[3].find("=") + 1:])
                if cur_line[4].startswith("span="):
                    span = int(cur_line[4][cur_line[4].find("=") + 1:])
                if chr_name in self.genome_size:
                    size = self.genome_size[chr_name]
                    self.genome[chr_name] = wigChrom.WigChrom(chr_name, start, size, step, span)
                else:
                    chr_name = "unknown"
                cur_position = start/step
            else:
                if chr_name == "unknown":
                    continue
                if cur_position < self.genome[chr_name].signals.shape[0]:
                    self.genome[chr_name].signals[cur_position] = float(cur_line[0])
                cur_position += 1
        wig.close()
        self.reload()

    def reload(self):
        """
        reload the chromosome to only keep the none zero value and the index
        :return:
        """
        for chr_name in self.genome.keys():
            self.genome[chr_name].reload()

    def callpeak(self, cutoff, min_width=40, output=None):
        peaks = []
        for chr_name in self.genome.keys():
            cur_peaks = self.genome[chr_name].get_peaks(None, None, cutoff, min_width)
            peaks += cur_peaks

        df = pd.DataFrame(peaks)
        df.columns = ['chr', 'start', 'end', 'center', 'width_above_cutoff', 'total_signal', 'height', 'height_logP',
                      'skewness', 'kurtosis']
        if output is None:
            df.to_csv(self.file_name+'_'+str(cutoff)+'.csv', index=None)
        else:
            df.to_csv(output+'_'+str(cutoff)+'.csv', index=None)
        return df
