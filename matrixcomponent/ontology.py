from typing import List
from rdflib import Namespace, Graph, Literal, URIRef, RDF


class Path:
    id: str
    ns: URIRef

    def __init__(self):
        ns = Namespace("") # should be empty!

    def __str__(self):
        return "path" + str(self.id)  # path1

    def ns_term(self):
        return self.ns + "/" + str(self)

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        path = self.ns + str(self)  # str representation

        # add the object itself
#        graph.add((path, RDF.type, faldo.Path))

        # add its properties, recursively if needed


class Region:
    ns: URIRef
    begin: int
    end: int

    def __str__(self):
        return "region/" + str(self.begin) + "-" + str(self.end)  # region/2-10

    def ns_term(self):
        return self.ns + str(self) # path1/...

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        region = self.ns_term()  # str representation

        # add the object itself
        graph.add( (region, RDF.type, faldo.Region) )

        # add its properties, recursively if needed
        graph.add( (region, faldo.begin, self.ns + "position/" + str(self.begin)) )
        graph.add( (region, faldo.end, self.ns + "position/" + str(self.end)) )


class Cell:
    id: int
    path_id: int
    ns: URIRef
    position_percent: float
    inversion_percent: float
    cell_region: List[Region] # [[1,4], [7,11]]

    def __init__(self):
        self.cell_region = []

    def __str__(self):
        return "cell" + str(self.id)  # cell1

    def ns_term(self):
        return self.ns + str(self)

    def add_to_graph(self, graph: Graph, vg, faldo: Namespace) -> None:
        cell = self.ns_term()
        inner_ns = URIRef("path" + str(self.path_id) + "/")

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
    path_id: int
    bin_edge: object # pointer to the next Bin
    cells: List[Region]

    def __init__(self):
        self.bin_edge = None
        self.cells = []

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
            cell.path_id = self.path_id
            graph.add((bin, vg.cells, cell.ns_term()))  # can have multiple bins
            cell.add_to_graph(graph, vg, faldo)

        # add the reference to another object if needed
        if self.bin_edge is not None:
            graph.add((inner_ns, vg.binEdge, Literal(self.bin_edge))) # __str__ will be called


class Link:
    id: int
    ns: URIRef
    arrival: Bin
    departure: Bin
    paths: List[Path]

    def __init__(self):
        self.paths = []
        self.arrival = None
        self.departure = None

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
            self.arrival.ns = self.ns
            graph.add((link, vg.arrival, self.arrival.ns_term()))

        if self.departure:
            self.departure.ns = self.ns
            graph.add((link, vg.departure, self.departure.ns_term()))

        for path in self.paths:
            graph.add((link, vg.linkPaths, URIRef("path" + str(path))))  # can have multiple bins


class Component:
    id: int
    ns: URIRef
    component_edge: object # Component
    component_rank: int
    bins: List[Bin]
    links: List[Link]

    def __init__(self, id):
        self.id = id
        self.component_edge = None
        self.bins = []
        self.links = []

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
            graph.add((component, vg.bins, bin.ns_term())) # can have multiple bins
            bin.add_to_graph(graph, vg, faldo) # add the inner content of each

        for link in self.links:
            link.ns = inner_ns
            link.add_to_graph(graph, vg, faldo) # add the inner content of each


        # add the reference to another object if needed
        if self.component_edge is not None:
            graph.add((component, vg.componentEdge, Literal(self.component_edge))) # __str__ will be called


class ZoomLevel:
    id: str
    ns: URIRef
    zoom_factor: int
    components: List[Component]

    def __init__(self):
        self.components = []

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
