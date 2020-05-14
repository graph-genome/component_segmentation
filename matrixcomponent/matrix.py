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


class Component:
    """Block of co-linear variation within a Graph Matrix
        # bin_id and seq are global to column and could be reduced to save memory,
        # careful construction can reuse Bin.sequence memory pointer"""
    first_bin: int
    last_bin: int
    occupants: set  # pure ids
    matrix: List
    arrivals: List[LinkColumn]
    departures: List[LinkColumn]
    x = 0

    def __init__(self, first_bin: int, last_bin: int):
        self.first_bin = first_bin
        self.last_bin = last_bin
        self.occupants = set()
        self.matrix = []
        self.arrivals = []  # reverse ordered Links
        self.departures = []  # ordered Links

    def width(self):
        return len(self.arrivals) + (len(self.departures) - 1) + self.last_bin - self.first_bin
    # (len(self.departures)-1) because last departure is adjacent connectors

    def next_x_coord(self):
        return self.x + self.width() + 1 # 1 column of padding
