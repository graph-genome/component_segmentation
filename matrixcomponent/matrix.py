"""Python object models to be manipulated"""

from dataclasses import dataclass
from typing import List
import numpy
import recordclass
from sortedcontainers import SortedDict


class Bin(recordclass.dataobject):
    bin_id: int
    coverage: float
    inversion: float
    nucleotide_ranges: 'numpy.array' # List[List[int]] is encoded as a Numpy flat array - this saves memory
    path_id: int

## Path is all for input files


class Path:
    name: str
    bins: 'SortedDict[int,Bin]'
    path_dividers: 'numpy.array'
    self_loops: 'numpy.array'
    num_links: int
    max_bin_id: int

    def __init__(self, name=''):
        self.name = name
        self.bins = SortedDict()  # Bin


## For Output to RDF  ###########
@dataclass
class LinkColumn(recordclass.dataobject):
    upstream: int
    downstream: int
    participants: 'numpy.array' # ids of participated path_names

    # todo: not more than 2^32 bins are supported - refactor the ontology code processing the Link items
    def __hash__(self):
        return (self.upstream << 32) + self.downstream


class Component:
    """Block of co-linear variation within a Graph Matrix
        # bin_id and seq are global to column and could be reduced to save memory,
        # careful construction can reuse Bin.sequence memory pointer"""
    first_bin: int
    last_bin: int
    x = 0
    compressedX = 0
    occupants: set  # pure ids
    matrix: List
    arrivals: List[LinkColumn]
    departures: List[LinkColumn]

    def __init__(self, first_bin: int, last_bin: int):
        self.first_bin = first_bin
        self.last_bin = last_bin
        if first_bin == 0:
            # don't confuse the special case bin 0 with regular content
            assert last_bin == 0
        self.occupants = set()
        self.matrix = []
        self.arrivals = []  # reverse ordered Links
        self.departures = []  # ordered Links

    def width(self, includes_connectors, compressed=False):
        depart_n = len(self.departures)
        if includes_connectors:
            depart_n -= 1
        matrix_width = (self.last_bin - self.first_bin + 1)
        if compressed:
            matrix_width = 3  # binScalingFactor
        return len(self.arrivals) + depart_n + matrix_width
    # (len(self.departures)-1) because last departure is adjacent connectors

    def next_x_coord(self, includes_connectors):
        # 1 column of padding
        return self.x + self.width(includes_connectors) + (1 if includes_connectors else 0)

    def next_compressedX(self, includes_connectors):
        # 1 column of padding
        return self.compressedX + self.width(includes_connectors, compressed=True) + (1 if includes_connectors else 0)
