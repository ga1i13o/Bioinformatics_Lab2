import numpy as np


class SequenceGenerator(object):
    def __init__(self, filename=".fa", n_reads=100, probs=None, seq_len=50):
        self.filename = filename
        self.reads = n_reads
        if probs is None:
            self.probs = [0.25, 0.25, 0.25, 0.25]
        else:
            self.probs = probs
        self.seq_len = seq_len
        self.bases = ["A", "T", "C", "G"]

    def generate_file(self):
        if self.filetype == "fa":
            self.__gen_fasta()
        elif self.filetype == "fq":
            self.__gen_fastq()

    def __gen_fasta(self):
        lines = []
        for i in range(self.reads):
            lines.append(self.__gen_fasta_id(i))
            sequence = np.random.choice(self.bases, size=self.seq_len, p=self.probs_list)
            sequence = "".join(sequence)
            lines.append(sequence+"\n")

        with open(self.filename, 'w') as output:
            for line in lines:
                output.write(line)

    def __gen_fastq(self):
        lines = []
        ascii_quality = np.arange(33, 127)
        prob_quality = 1/(127-33)*np.ones(127-33)
        for i in range(self.reads):
            lines.append(self.__gen_fastq_id(i, pos="begin", length=self.seq_len))
            basis = np.random.choice(self.bases, size=self.seq_len, p=self.probs_list)
            basis = "".join(basis)
            lines.append(basis+"\n")
            lines.append(self.__gen_fastq_id(i, pos="end", length=self.seq_len))
            ascii_qualities = np.random.choice(ascii_quality, size=self.seq_len, p=prob_quality)
            str_qualities = map(lambda x: chr(x), ascii_qualities)
            qualities = "".join(str_qualities)
            lines.append(qualities+"\n")

        with open(self.filename, 'w') as output:
            for line in lines:
                output.write(line)

    @staticmethod
    def __gen_fasta_id(n_read):
        pref = ">read_"
        return pref+str(n_read)+"\n"

    @staticmethod
    def __gen_fastq_id(n_read, pos, length):
        if pos == "begin":
            pref = "@"
        elif pos == "end":
            pref = "+"
        id = "read_" + str(n_read) + " length="+str(length)
        return pref+id+"\n"

    @property
    def probs(self):
        return self.__probs

    @probs.setter
    def probs(self, values):
        self.__probs = {"A":values[0], "T":values[1], "C":values[2], "G":values[3]}
        self.probs_list = values

    @property
    def filename(self):
        return self.__filename

    @filename.setter
    def filename(self, value):
        self.__filename = value
        if "." in value:
            self.filetype = value.split(".")[1]