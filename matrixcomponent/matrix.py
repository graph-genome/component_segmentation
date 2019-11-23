"""Python object models to be manipulated"""

from dataclasses import dataclass
from typing import List, Any, Set


## Path is all for input files

class Path:
    name: str
    links: 'List[Path.LinkEntry]'
    bins: 'List[Path.Bin]'

    def __init__(self, name=''):
        self.name = name
        self.bins = []  # Bin
        self.links = []  # LinkEntry
        self.__bin_set = set()

    class Bin:
        def __init__(self, bin_id, coverage, inversion_rate, mean_pos, sequence=''):
            self.bin_id = bin_id
            self.coverage = coverage
            self.inversion_rate = inversion_rate
            self.mean_pos = mean_pos
            self.sequence = sequence

    class LinkEntry:
        def __init__(self, upstream, downstream):
            self.upstream = upstream
            self.downstream = downstream
            # TODO: self.insert_size will require a topology search to find this

    def __contains__(self, item):  # used by " x in Path "
        return item in self.__bin_set

    def finalize_bins(self):
        self.__bin_set = {x.bin_id for x in self.bins}  # build and cache a set



## For Output to RDF
@dataclass
class LinkColumn:
    upstream: int
    downstream: int
    participants: Set[str]


class Component:
    """Block of co-linear variation within a Graph Matrix"""
    arrivals: List[LinkColumn]
    departures: List[LinkColumn]

    def __init__(self, first_bin: int, last_bin: int):
        self.first_bin = first_bin
        self.last_bin = last_bin
        # bin_id and seq are global to column and could be reduced to save memory,
        # careful construction can reuse Bin.sequence memory pointer
        self.arrivals = []  # reverse ordered Links
        self.departures = []  # ordered Links

    def split(self, abs_index):
        """Splits one Component into two
        :param abs_index the first position of the new component"""
        pass


@dataclass
class PangenomeSchematic:
    components: List[Component]


