import os
from .read_objects import FastaObj
from .utils import binary_search


class Comparator(object):
    def __init__(self, input_1="", input_2="", output="out.txt"):
        self.input_1 = input_1
        self.input_2 = input_2
        self.output_file = output
        self.is_input_loaded = 0
        self.loaded_file = []
        self.output = []

    def load_file(self):
        if self.is_input_loaded:
            print("The smaller file of the 2 is already loaded")
            return
        print("[!!] Loading file... [!!]")
        size1 = os.path.getsize(self.input_1)
        size2 = os.path.getsize(self.input_2)
        if size1 < size2:
            self.__load_smaller_file(self.input_1)
            self.is_input_loaded = 1
            print(f"[!!] File {self.input_1}, being the smallest, was loaded [!!]")
        else:
            self.__load_smaller_file(self.input_2)
            self.is_input_loaded = 2
            print(f"[!!] File {self.input_2}, being the smallest, was loaded [!!]")

    def compare_files(self):
        print("[!!] Comparing... [!!]")
        if not self.is_input_loaded:
            self.load_file()
        if self.is_input_loaded == 1:
            file_to_read = self.input_2
        else:
            file_to_read = self.input_1

        with open(file_to_read, "r") as in_file:
            for line_id in in_file:
                seq = next(in_file)
                read = FastaObj(line_id, seq)
                search = binary_search(self.loaded_file, read)
                if search > -1:
                    self.output.append( FastaObj(read.id+self.loaded_file[search].id,
                                         read.read) )

    def write_output(self):
        print(f"[!!] Output dumped to {self.output_file} [!!]")
        with open(self.output_file, "w") as out_file:
            for read in self.output:
                out = ">" + read.id + '\n' + read.read + '\n'
                out_file.write(out)

    def __load_smaller_file(self, file):
        with open(file, "r") as in_file:
            lines = in_file.readlines()
        for i in range(0, len(lines), 2):

            read = FastaObj(lines[i], lines[i+1])
            self.loaded_file.append(read)
        self.loaded_file.sort(key=lambda x: x.read)

    @property
    def input_1(self):
        return self.__input_1

    @input_1.setter
    def input_1(self, value):
        self.__input_1 = value
        if "." in value:
            self.filetype = value.split(".")[1]

    @property
    def input_2(self):
        return self.__input_2

    @input_2.setter
    def input_2(self, value):
        self.__input_2 = value
        if "." in value:
            self.filetype = value.split(".")[1]