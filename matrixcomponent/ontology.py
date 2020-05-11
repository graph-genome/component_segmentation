from typing import List
from rdflib import Namespace, Graph, Literal, URIRef, RDF


class Region:
    ns: URIRef
    begin: int
    end: int

    def __str__(self):
        suffix = str(self.begin)  # region/2
        if self.begin != self.end:
            suffix = suffix + "-" + str(self.end)  # region/2-10
        return "region/" + suffix

    def ns_term(self):
        return self.ns + str(self) # path1/...

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        region = self.ns_term()  # str representation

        if self.begin == self.end:
            # add the object itself
            graph.add((region, RDF.type, faldo.ExactPosition))
            graph.add((region, faldo.position, Literal(self.begin)))
        else:
            # add the object itself
            graph.add( (region, RDF.type, faldo.Region) )

            # add its properties, recursively if needed
            graph.add( (region, faldo.begin, self.ns + "position/" + str(self.begin)) )
            graph.add( (region, faldo.end, self.ns + "position/" + str(self.end)) )


class Cell:
    id: int
    path_id: str
    ns: URIRef
    position_percent: float
    inversion_percent: float
    cell_region: List[Region] # [[1,4], [7,11]]

    def __init__(self):
        self.cell_region = []

    def __str__(self):
        return "cell" + self.path_id # cell1

    def ns_term(self):
        return self.ns + str(self)

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        cell = self.ns_term()
        inner_ns = URIRef(self.path_id + "/")

        # add the object itself
        graph.add((cell, RDF.type, vg.Cell))

        # add its properties, recursively if needed
        graph.add((cell, vg.positionPercent, Literal(self.position_percent)))
        graph.add((cell, vg.inversionPercent, Literal(self.inversion_percent)))
        for region in self.cell_region:
            region.ns = inner_ns
            graph.add((cell, vg.cellRegion, region.ns_term()))  # can have multiple bins
            region.add_to_graph(graph, vg, faldo)


class Bin:
    ns: URIRef
    bin_rank: int
    forward_bin_edge: str
    reverse_bin_edge: str
    cells: List[Region]

    def __init__(self):
        self.cells = []
        self.forward_bin_edge = ''
        self.reverse_bin_edge = ''

    def __str__(self):
        return "bin" + str(self.bin_rank)  # bin1

    def ns_term(self):
        return self.ns + str(self)

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        bin = self.ns_term()
        inner_ns = bin + "/"

        # add the object itself
        graph.add((bin, RDF.type, vg.Bin))

        # add its properties, recursively if needed
        graph.add((bin, vg.binRank, Literal(self.bin_rank)))
        for cell in self.cells:
            cell.ns = inner_ns
            graph.add((bin, vg.cells, cell.ns_term()))  # can have multiple bins
            cell.add_to_graph(graph, vg, faldo)

        # add the reference to another object if needed
        if self.forward_bin_edge:
            graph.add((bin, vg.forwardBinEdge, URIRef(self.forward_bin_edge)))

        if self.reverse_bin_edge:
            graph.add((bin, vg.reverseBinEdge, URIRef(self.reverse_bin_edge)))


class Link:
    id: int
    ns: URIRef
    arrival: str
    departure: str
    paths: List[str]
    forward_link_edge: int
    reverse_link_edge: int

    def __init__(self):
        self.paths = []
        self.forward_link_edge = -1
        self.reverse_link_edge = -1
        self.arrival = ''
        self.departure = ''

    def __str__(self):
        return "link" + str(self.id)  # bin1


    def ns_term(self):
        return self.ns + str(self)


    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        link = self.ns_term()
        inner_ns = link + "/"

        # add the object itself
        graph.add((link, RDF.type, vg.Link))

        # add its properties, recursively if needed
        if self.arrival:
            graph.add((link, vg.arrival, URIRef(self.arrival)))

        if self.departure:
            graph.add((link, vg.departure, URIRef(self.departure)))

        for path in self.paths:
            graph.add((link, vg.linkPaths, URIRef(path)))  # can have multiple bins

        # add the reference to another object if needed
        if self.forward_link_edge > -1:
            graph.add((link, vg.forwardLinkEdge, self.ns + "link" + str(self.forward_link_edge)))

        if self.reverse_link_edge > -1:
            graph.add((link, vg.reverseLinkEdge, self.ns + "link" + str(self.reverse_link_edge)))


class Component:
    id: int
    ns: URIRef
    forward_component_edge: str  # id of the next Component
    reverse_component_edge: str  # id of the previous Component
    component_rank: int
    bins: List[Bin]

    def __init__(self, id):
        self.id = id
        self.bins = []
        self.forward_component_edge = ''
        self.reverse_component_edge = ''

    def __str__(self):
        return "component" + str(self.id) # component1

    def ns_term(self):
        return self.ns + self.__str__()

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        component = self.ns_term()
        inner_ns = component + "/"

        # add the object itself
        graph.add((component, RDF.type, vg.Component))

        # add its properties, recursively if needed
        graph.add((component, vg.componentRank, Literal(self.id)))
        for bin in self.bins:
            bin.ns = inner_ns
            graph.add((component, vg.bins, bin.ns_term()))  # can have multiple bins
            bin.add_to_graph(graph, vg, faldo)  # add the inner content of each

        # add the reference to another object if needed
        if self.forward_component_edge:
            graph.add((component, vg.forwardComponentEdge, URIRef(self.forward_component_edge)))

        if self.reverse_component_edge:
            graph.add((component, vg.reverseComponentEdge, URIRef(self.reverse_component_edge)))


class ZoomLevel:
    id: str
    ns: URIRef
    zoom_factor: int
    components: List[Component]
    links: List[Link]

    def __init__(self):
        self.components = []
        self.links = []

    def __str__(self):
        return "zoom" +  str(self.zoom_factor) #zoom1000

    def ns_term(self):
        return self.ns + self.__str__()

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        zoomfactor = self.ns_term()
        inner_ns = zoomfactor + "/"

        # add the object itself
        graph.add((zoomfactor, RDF.type, vg.ZoomLevel))

        # add its properties, recursively if needed
        graph.add((zoomfactor, vg.zoomFactor, Literal(self.zoom_factor)))
        for i,comp in enumerate(self.components):
            comp.ns = inner_ns
            graph.add((zoomfactor, vg.components, comp.ns_term()))
            comp.add_to_graph(graph, vg, faldo)

        for link in self.links:
            link.ns = inner_ns
            link.add_to_graph(graph, vg, faldo)
