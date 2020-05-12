import json
from collections import OrderedDict
from statistics import mean

import math
from bisect import bisect

from math import ceil
from typing import List

from dataclasses import dataclass

from rdflib import URIRef, Graph, Namespace

from matrixcomponent import JSON_VERSION, ontology
from matrixcomponent.matrix import Component, Bin, LinkColumn

from DNASkittleUtils.Contigs import Contig, write_contigs_to_file

@dataclass
class PangenomeSchematic:
    json_version: int
    bin_width: int
    first_bin: int
    last_bin: int
    components: List[Component]
    path_names: List[str]
    total_nr_files: int
    pangenome_length: int
    file_dict = {}  # dict for bin2file.json

    def json_dump(self):
        def dumper(obj):
            if isinstance(obj, Bin):  # should be in Bin class def
                return [obj.coverage, obj.inversion, obj.nucleotide_ranges]
            if isinstance(obj, LinkColumn): # todo: get rid of booleans once the JS side can digest path_names ids
                bools = [False] * obj.num_paths
                for i in obj.participants:
                    bools[i] = True
                return {'upstream':obj.upstream, 'downstream':obj.downstream, 'participants':bools}

            if isinstance(obj, set):
                return list(obj)
            try:
                return obj.__dict__
            except:
                return obj
        # whitespace and [] takes up the majority of the file size
        return json.dumps(self, default=dumper, indent=None, separators=(',', ':\n'))

    def n_links(self):
        return sum([len(x.arrivals) + len(x.departures) for x in self.components])

    def n_components(self):
        return len(self.components)

    def update_first_last_bin(self):
        self.first_bin = 1  # these have not been properly initialized
        self.last_bin = self.components[-1].last_bin

    def split_and_write(self, cells_per_file, folder, fasta : Contig, ontology_folder):
        """Splits one Schematic into multiple files with their own
        unique first and last_bin based on the volume of data desired per
        file specified by cells_per_file.  """
        avg_paths = self.lazy_average_occupants()
        columns_per_file = ceil(cells_per_file / avg_paths)
        column_counts = [c.x for c in self.components]
        cut_points = self.find_cut_points_in_file_split(columns_per_file, column_counts)

        # variables cut and end_cut are componentIDs
        # binIDs are in components.{first,last}_bin
        bin2file_mapping = []
        for i, cut in enumerate(cut_points[:-1]):
            end_cut = cut_points[i + 1]
            these_comp = self.components[cut:end_cut]
            if these_comp:  # when we're out of data because the last component is 1 wide
                schematic = PangenomeSchematic(JSON_VERSION, self.bin_width, these_comp[0].first_bin,
                                               these_comp[-1].last_bin, these_comp, self.path_names,
                                               self.total_nr_files, self.pangenome_length)
                schematic.filename = self.filename(i)  # save for consistency IMPORTANT

                chunk_summary = {"file": schematic.filename,
                                 "first_bin": schematic.first_bin,
                                 "last_bin": schematic.last_bin,
                                 # used to calculate amount of padding for x values
                                 "component_count": schematic.n_components(),
                                 # extra columns beyond last_bin - first_bin
                                 "link_count": schematic.n_links()}
                if schematic.bin_width == 1:
                    chunk_summary["fasta"] = self.fasta_filename(i)
                bin2file_mapping.append(chunk_summary)

                p = folder.joinpath(schematic.filename)
                with p.open('w') as fpgh9:
                    fpgh9.write(schematic.json_dump())

                if fasta is not None:
                    x = schematic.bin_width
                    fa_first, fa_last = (schematic.first_bin * x), ((schematic.last_bin + 1) * x)
                    header = f"first_bin: {schematic.first_bin} " + f"last_bin: {schematic.last_bin}"
                    chunk = [Contig(header, fasta.seq[fa_first:fa_last])]
                    c = folder.joinpath(schematic.fasta_filename(i))
                    write_contigs_to_file(c, chunk)

                if ontology_folder:
                    zoom_level = ontology.ZoomLevel()
                    zoom_level.zoom_factor = schematic.bin_width
                    zoom_level.ns = URIRef('pg/')

                    prev_comp_id = -1
                    cell_counter = 0
                    ocomp_dict = {}
                    obin_dict = {}
                    oposition_dict = {}
                    for ic, component in enumerate(schematic.components):
                        ocomp = ontology.Component(ic+1)
                        ocomp.ns = zoom_level.ns_term() + '/'
                        zoom_level.components.append(ocomp)

                        # save the sequence 1-2-3-..-n as a bi-directed list
                        if prev_comp_id in ocomp_dict:
                            prev_comp = ocomp_dict[prev_comp_id]
                            ocomp.reverse_component_edge = prev_comp.ns_term()
                            prev_comp.forward_component_edge = ocomp.ns_term()

                        ocomp_dict[ic] = ocomp
                        prev_comp_id = ic

                        # bins
                        for bins in component.matrix:
                            prev_bin_id = -1
                            for bin in bins:
                                if bin:
                                    cur_bin_id = bin.bin_id
                                    obin = ontology.Bin()
                                    obin.ns = ocomp.ns_term() + '/'
                                    obin.bin_rank = cur_bin_id
                                    obin_dict[cur_bin_id] = obin

                                    # save the sequence 1-2-3-..-n as a bi-directed list
                                    if prev_bin_id in obin_dict:
                                        prev_bin = obin_dict[prev_bin_id]
                                        prev_bin.forward_bin_edge = obin.ns_term()
                                        obin.reverse_bin_edge = prev_bin.ns_term()

                                    prev_bin_id = cur_bin_id
                                    ocomp.bins.append(obin)

                                    cell_counter = cell_counter + 1
                                    ocell = ontology.Cell()
                                    ocell.id = cell_counter
                                    ocell.path_id = self.path_names[bin.path_id]  # saved in the populate_component_matrix
                                    ocell.inversion_percent = bin.inversion
                                    ocell.position_percent = bin.coverage

                                    # todo: are begin,end the real bin_ids or the compressed ones? a sparse list sense
                                    cell_ns = URIRef(ocell.path_id + "/")
                                    for [begin, end] in bin.nucleotide_ranges:
                                        oregion = ontology.Region()
                                        oregion.begin = begin
                                        oregion.end = end
                                        ocell.cell_region.append(oregion)

                                        oposition_begin = ontology.Position(begin, cell_ns)
                                        oposition_end = ontology.Position(end, cell_ns)
                                        oposition_dict[oposition_begin.ns_term()] = oposition_begin
                                        oposition_dict[oposition_end.ns_term()] = oposition_end

                                    obin.cells.append(ocell)

                    # links between components and their bins

                    olink_dict = {}
                    link_counter = 0
                    for component in schematic.components:
                        # search in all arrivals component <-> component links; departures are iterated automatically
                        # every link from departures will be in some other arrival
                        for link in component.departures:
                            if len(link.participants):
                                link_counter = link_counter + 1
                                olink = ontology.Link()
                                olink.id = link_counter
                                olink_dict[link_counter] = olink


                    link_counter = 0
                    for component in schematic.components:
                        # search in all arrivals component <-> component links; departures are iterated automatically
                        # every link from departures will be in some other arrival
                        for link in component.departures:
                            if len(link.participants):
                                link_counter = link_counter + 1
                                olink = olink_dict[link_counter]

                                if link.upstream in obin_dict:
                                    from_bin = obin_dict[link.upstream]
                                    olink.departure = from_bin.ns_term()

                                if link.downstream in obin_dict:
                                    to_bin = obin_dict[link.downstream]
                                    olink.arrival = to_bin.ns_term()

                                olink.paths = [self.path_names[k] for k in link.participants]
                                zoom_level.links.append(olink)

                    g = Graph()
                    vg = Namespace('http://biohackathon.org/resource/vg#')
                    faldo = Namespace('http://biohackathon.org/resource/faldo#')
                    g.bind('vg', vg)
                    g.bind('faldo', faldo)

                    zoom_level.add_to_graph(g, vg, faldo)  # here the magic happens
                    for oposition in oposition_dict.values():
                        oposition.add_to_graph(g, vg, faldo)

                    p = ontology_folder.joinpath(schematic.ttl_filename(i))
                    g.serialize(destination=str(p), format='turtle', encoding='utf-8')


        return bin2file_mapping

    def find_cut_points_in_file_split(self, columns_per_file, column_counts):
        """Use binary search bisect to find the component boundaries closest to
        target breakpoints for splitting files."""
        cut_points, prev_point = [0], 0
        for start_bin in range(0, column_counts[-1], columns_per_file):
            cut = bisect(column_counts, start_bin + columns_per_file)
            cut_points.append(max(prev_point + 1, cut))
            prev_point = cut
        cut_points.append(len(self.components))  # don't chop of dangling end
        self.total_nr_files = len(cut_points) - 1
        return cut_points

    def lazy_average_occupants(self):
        """grab four random components and check how many occupants they have"""
        samples = [self.components[int(len(self.components) * (perc/100))] for perc in range(1, 99)]
        avg_paths = mean([sum(x.occupants) for x in samples])
        return avg_paths

    def pad_file_nr(self, file_nr):
        return str(file_nr).zfill(int(math.log10(self.total_nr_files) + 1))

    def filename(self, nth_file):
        return f'chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.schematic.json'

    def fasta_filename(self, nth_file):
        return f'seq_chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.fa'

    def ttl_filename(self, nth_file):
        return f'seq_chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.ttl'

    def write_index_file(self, folder, bin2file_mapping):

        file_contents = {'bin_width': self.bin_width,
                         'last_bin': self.last_bin,
                         'files': bin2file_mapping}
        self.file_dict[f'{self.bin_width}'] = file_contents

        master_index_file = folder.parent.joinpath(f'bin2file.json')
        master_contents = {'json_version': JSON_VERSION,
                           'pangenome_length': self.pangenome_length,
                           'zoom_levels': OrderedDict(sorted(self.file_dict.items(), reverse=True))}
        with master_index_file.open('w') as f:
            f.write(json.dumps(master_contents, indent=4))
            print("Saved file2bin mapping to", master_index_file)

    def prerender(self):
        """Calculates X coordinates and summary statistics for all Components"""
        x = 0
        for component in self.components:
            component.x = x
            # component.first_bin=0 does not take up rendering space, the next component is 0
            x = component.next_x_coord() if component.first_bin else 0
        self.update_first_last_bin()

