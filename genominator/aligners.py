from .read_objects import FastaObj
import numpy as np
from statistics import mode


class LocalAligner(object):
    def __init__(self, match=1, mismatch=0, gap=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.score_matrix = None
        self.predecessor_matrix = None
        self.max_score = 0
        self.chunks_aligned = []

    def align(self, seq_a, seq_b):
        print("[!!] Aligning... [!!]")
        n_rows = len(seq_a) + 1
        n_cols = len(seq_b) + 1
        self.score_matrix = np.zeros(shape=(n_rows, n_cols), dtype=int)
        self.predecessor_matrix = np.zeros(shape=(n_rows, n_cols), dtype=object)

        for i in range(n_rows):
            self.score_matrix[i, 0] = 0
        for j in range(n_cols):
            self.score_matrix[0, j] = 0

        for i in range(1, n_rows):
            for j in range(1, n_cols):
                match = self.score_matrix[i-1, j-1] + self.__sim(seq_a[i-1], seq_b[j-1])
                delete = self.score_matrix[i-1, j] + self.gap
                insert = self.score_matrix[i, j-1] + self.gap
                self.score_matrix[i, j] = self.__choose_move(match, insert, delete, i, j)
        self.max_score = self.score_matrix.max()
        self.__rewind_moves(seq_a, seq_b)
        self.__write_output()
        print("[!!] All done [!!]")

    def __write_output(self):
        print(f"Local alignment score of each chunk:\t{self.max_score}")
        print(self.score_matrix)
        print(f"\nList of found local alignments")
        for chunk in self.chunks_aligned:
            links = "".join(['|' if chunk[0][i] == chunk[1][i] else ' ' for i in range(len(chunk[0]))])
            print(f"\t{chunk[0]}\n\t{links}\n\t{chunk[1]}")

    def __rewind_moves(self, seq1, seq2):
        print("[!!] Backtracking... [!!]")
        local_chunks = np.argwhere( self.score_matrix == self.max_score)

        for chunk in local_chunks:
            seq1_align = ""
            seq2_align = ""
            i = chunk[0]
            j = chunk[1]
            while self.score_matrix[i, j] != 0:
                move = self.predecessor_matrix[i, j]
                if move == "match":
                    seq1_align += seq1[i-1]
                    seq2_align += seq2[j-1]
                    i, j = i-1, j-1
                elif move == "insert":
                    seq1_align += seq1[i-1]
                    seq2_align += "_"
                    i -= 1
                elif move == "delete":
                    seq1_align += "_"
                    seq2_align += seq2[j-1]
                    j -= 1
            self.chunks_aligned.append((seq1_align[::-1], seq2_align[::-1]))

    def __choose_move(self, match, delete, insert, i, j):
        max_score = max(match, delete, insert)
        if 0 >= max_score:
            self.predecessor_matrix[i, j] = "skip"
            return 0
        if match == max_score:
            self.predecessor_matrix[i, j] = "match"
            return match
        if delete == max_score:
            self.predecessor_matrix[i, j] = "delete"
            return delete
        self.predecessor_matrix[i, j] = "insert"
        return insert

    def __sim(self, x, y):
        if x == y:
            return self.match
        else:
            return self.mismatch


class GlobalAligner(object):
    def __init__(self, match=1, mismatch=0, gap=-1):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.score_matrix = None
        self.predecessor_matrix = None

    def align(self, seq_a, seq_b):
        print("[!!] Aligning... [!!]")
        n_rows = len(seq_a) + 1
        n_cols = len(seq_b) + 1
        self.score_matrix = np.zeros(shape=(n_rows, n_cols), dtype=int)
        self.predecessor_matrix = np.zeros(shape=(n_rows, n_cols), dtype=object)

        for i in range(n_rows):
            self.score_matrix[i, 0] = self.gap * i
        for j in range(n_cols):
            self.score_matrix[0, j] = self.gap * j

        for i in range(1, n_rows):
            for j in range(1, n_cols):
                match = self.score_matrix[i-1, j-1] + self.__sim(seq_a[i-1], seq_b[j-1])
                delete = self.score_matrix[i-1, j] + self.gap
                insert = self.score_matrix[i, j-1] + self.gap
                self.score_matrix[i, j] = self.__choose_move(match, insert, delete, i, j)
        result = self.__rewind_moves(seq_a, seq_b)
        self.__write_output(result)
        print("[!!] All done [!!]")

    def __write_output(self, result):
        links = "".join(['|' if result[0][i] == result[1][i] else ' ' for i in range(len(result[0]))])
        print(f"Global alignment score:\t{self.score_matrix[-1, -1]}")
        print(self.score_matrix)
        print(f"\nFound global alignment:\n\t{result[0]}\n\t{links}\n\t{result[1]}")

    def __rewind_moves(self, seq1, seq2):
        print("[!!] Backtracking... [!!]")
        i = len(seq1)
        j = len(seq2)
        seq1_align = ""
        seq2_align = ""

        while i > 0 or j > 0:
            move = self.predecessor_matrix[i, j]
            if move == "match":
                seq1_align += seq1[i-1]
                seq2_align += seq2[j-1]
                i, j = i-1, j-1
            elif move == "insert":
                seq1_align += seq1[i-1]
                seq2_align += "_"
                i -= 1
            elif move == "delete":
                seq1_align += "_"
                seq2_align += seq2[j-1]
                j -= 1


        return seq1_align[::-1], seq2_align[::-1]

    def __choose_move(self, match, delete, insert, i, j):
        if match == max(match, delete, insert):
            self.predecessor_matrix[i, j] = "match"
            return match
        if delete == max(match, delete, insert):
            self.predecessor_matrix[i, j] = "delete"
            return delete
        self.predecessor_matrix[i, j] = "insert"
        return insert

    def __sim(self, x, y):
        if x == y:
            return self.match
        else:
            return self.mismatch


class ConsensusFinder(object):
    def __init__(self, in_file=""):
        self.input_file = in_file
        self.fragments = []
        self.consensus = []

    def find_consensus(self):
        self.__load_fragments()
        self.__group_fragments()
        print("[!!] All done [!!]")

    def write_output(self, output_file="out_consensus.txt"):
        print(f"[!!] Output dumped to {output_file} [!!]")
        with open(output_file, "w") as out_file:
            for read in self.consensus:
                out = read + '\n'
                out_file.write(out)

    def __group_fragments(self):
        print("[!!] Constructing consensus... [!!]")
        begin_region = 0
        end_region = 1
        for i in range(1, len(self.fragments)):
            if self.fragments[i][1] >= self.fragments[i-1][1] + self.fragments[i-1][0].length:
                end_region = i
                self.consensus.append(self.__extract_consensus(begin_region, end_region))
                begin_region = i
        if end_region != len(self.fragments) - 1:
            end_region = len(self.fragments)
            self.consensus.append(self.__extract_consensus(3, end_region))

    def __extract_consensus(self, begin, end):
        n_fragments = end-begin
        consensus_start_position = self.fragments[begin][1]
        consensus_length = self.fragments[end-1][1] +\
                           self.fragments[end-1][0].length -\
                           self.fragments[begin][1]
        consensus_result = ""
        consensus_matrix = np.zeros(shape=(n_fragments,
                                            consensus_length), dtype=int)

        for i in range(n_fragments):
            consensus_matrix[i][self.fragments[begin+i][1] - consensus_start_position:\
                                self.fragments[begin+i][1]+self.fragments[begin+i][0].length - consensus_start_position] =\
                        [ord(base) for base in self.fragments[begin+i][0].read]
        for i in range(consensus_length):
            if consensus_matrix[:,i].sum() % consensus_matrix[:,i].max() != 0:
                consensus_result += chr(mode(consensus_matrix[:,i]))
            else:
                consensus_result += chr(consensus_matrix[:,i].max())
        return consensus_result

    def __load_fragments(self):
        print("[!!] Loading fragments... [!!]")
        with open(self.input_file, "r") as in_file:
            lines = in_file.readlines()
        for line in lines:
            fields = line.split(" ")
            fasta_read = FastaObj(id_line=fields[0], read=fields[1])
            self.fragments.append((fasta_read, int(fields[2].strip())))
        self.fragments.sort(key=lambda x: x[1])

    @property
    def input_file(self):
        return self.__input_file

    @input_file.setter
    def input_file(self, value):
        self.__input_file = value
        if "." in value:
            self.filetype = value.split(".")[1]