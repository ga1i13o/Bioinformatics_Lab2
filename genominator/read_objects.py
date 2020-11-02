from functools import total_ordering


@total_ordering
class FastaObj(object):
    def __init__(self, id_line="", read=""):
        begin = 0
        if id_line.startswith(">"):
            begin = 1
        self.id = id_line[begin:].strip()
        self.read = read.strip()
        self.length = len(self.read)

    def __eq__(self, other):
        if not isinstance(other, FastaObj):
            # don't attempt to compare against unrelated types
            return NotImplemented
        if self.length == other.length:
            if self.read == other.read:
                return True
        return False

    def __lt__(self, other):
        if not isinstance(other, FastaObj):
            return NotImplemented
        return self.read < other.read

    def __str__(self):
        return self.id + " " + self.read

    def __repr__(self):
        return self.id + " " + self.read