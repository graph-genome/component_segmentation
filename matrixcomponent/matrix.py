"""Python object models to be manipulated"""

from dataclasses import dataclass
from typing import List
import numpy
from sortedcontainers import SortedDict


@dataclass
class Bin:
    bin_id: int
    coverage: float
    inversion: float
    nucleotide_ranges: List[List[int]]
    sequence: str = ''

## Path is all for input files


class Path:
    name: str
    bins: 'SortedDict[int,Bin]'
    links: 'numpy.array'

    def __init__(self, name=''):
        self.name = name
        self.bins = SortedDict()  # Bin


## For Output to RDF  ###########
@dataclass
class LinkColumn:
    upstream: int
    downstream: int
    num_paths: int
    participants: 'numpy.array' # ids of participated path_names


class Component:
    """Block of co-linear variation within a Graph Matrix
        # bin_id and seq are global to column and could be reduced to save memory,
        # careful construction can reuse Bin.sequence memory pointer"""
    first_bin: int
    last_bin: int
    occupants: List[bool]
    matrix: List[List[Bin]]
    arrivals: List[LinkColumn]
    departures: List[LinkColumn]
    x = 0

    def __init__(self, first_bin: int, last_bin: int):
        self.first_bin = first_bin
        self.last_bin = last_bin
        self.occupants = []
        self.matrix = []
        self.arrivals = []  # reverse ordered Links
        self.departures = []  # ordered Links

    def width(self):
        return len(self.arrivals) + (len(self.departures)-1) + self.last_bin - self.first_bin
    # (len(self.departures)-1) because last departure is adjacent connectors

    def next_x_coord(self):
        return self.x + self.width() + 1 # 1 column of padding
