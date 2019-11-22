
# Python object models to be manipulated


class Path:
    def __init__(self):
        self.bins = []  # Bin
        self.links = []  # LinkEntry

    class Bin:
        def __init__(self, bin_id, coverage, inversion_rate, mean_pos):
            self.bin_id = bin_id
            self.coverage = coverage
            self.inversion_rate = inversion_rate
            self.mean_pos = mean_pos

    class LinkEntry:
        def __init__(self, left, right):
            self.left = left
            self.right = right


class Component:
    """Block of co-linear variation within a Graph Matrix"""
    def __init__(self):
        self.matrix = [[]]  # square grid
        self.upstream = []  # Links
        self.downstream = []  # Links

    def split(self, abs_index):
        """Splits one Component into two"""
        pass


def test_construction():
    p = Path()
    p.bins = [Path.Bin(1,.9, 0, 0.11), Path.Bin(2,.8, 0, 0.15),]
    p.links = [Path.LinkEntry(1,2)]

