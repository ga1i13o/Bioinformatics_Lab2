class Statista(object):
    def __init__(self, input_file="", gc_thresh=0):
        self.input_file = input_file
        self.gc_thresh = gc_thresh
        self.basis = ["A", "T", "C", "G"]
        self.line_filterer = {"fa": self.__fa_valid_line, "fq": self.__fq_valid_line}
        self.available_operations = {"basis_number":self.__count_basis, "low_complex_seqs": self.__count_low_complex,
                                     "gc_content": self.__gc_content}
        self.ops_to_perform = []
        self.output_data = {k:None for k in self.available_operations.keys()}

    def set_operations(self, ops):
        for op in ops:
            if op in self.available_operations.keys():
                self.ops_to_perform.append(op)
            else:
                print(f"Operation '{op}' not supported")

    def compute_stats(self):
        with open(self.input_file) as in_file:
            prev_line = ""
            for line in in_file:
                if self.line_filterer[self.filetype](line, prev_line):
                    for op in self.ops_to_perform:
                        self.available_operations[op](line, prev_line)
                prev_line = line

    def write_output(self, out_file):
        with open(out_file, "w") as out_file:
            for k,v in self.output_data.items():
                out = k + " : " +str(v) +"\n"
                out_file.write(out)

    def __fa_valid_line(self, line, _):
        if line.startswith(">"):
            return False
        return True

    def __fq_valid_line(self, _, prev_line):
        if prev_line.startswith("@"):
            return True
        return False

    def __count_basis(self, line, _):
        if self.output_data["basis_number"] is None:
            self.output_data["basis_number"] = {"A":0, "T":0, "C":0, "G":0}

        for base in self.basis:
            self.output_data["basis_number"][base] += line.count(base)

    def __count_low_complex(self, line, _):
        if self.output_data["low_complex_seqs"] is None:
            self.output_data["low_complex_seqs"] = 0

        reads_count = 0
        for base in self.basis:
            reads_count += line.count(5*base)
        if reads_count > 0:
            self.output_data["low_complex_seqs"] += 1

    def __gc_content(self, line, prev_line):
        if self.output_data["gc_content"] is None:
            self.output_data["gc_content"] = []
        gc_count = line.count("GC")
        if gc_count > self.gc_thresh:
            self.output_data["gc_content"].append((self.__extract_id(prev_line), gc_count))

    def __extract_id(self, line):
        return line[1:].strip()

    @property
    def gc_thresh(self):
        return self.__gc_thresh

    @gc_thresh.setter
    def gc_thresh(self, val):
        if val is None:
            val = 0
        self.__gc_thresh = val

    @property
    def input_file(self):
        return self.__input_file

    @input_file.setter
    def input_file(self, value):
        self.__input_file = value
        if "." in value:
            self.filetype = value.split(".")[1]

