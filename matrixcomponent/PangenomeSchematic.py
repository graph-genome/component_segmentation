import gzip
import json
import logging
import os
import shutil

from collections import OrderedDict
from statistics import mean

import math
from bisect import bisect

from math import ceil
from typing import List

from dataclasses import dataclass

from joblib import delayed
from rdflib import URIRef, Graph, Namespace

from matrixcomponent import JSON_VERSION, ontology
from matrixcomponent.matrix import Component, Bin, LinkColumn

from DNASkittleUtils.Contigs import Contig, write_contigs_to_file

LOGGER = logging.getLogger(__name__)
"""logging.Logger: The logger for this module"""

# move out of the class to be able to use with joblib
def write_rdf(schematic, ontology_output_path):
    zoom_level = ontology.ZoomLevel()
    zoom_level.zoom_factor = schematic.bin_width
    zoom_level.ns = URIRef('http://example.org/vg/')

    prev_comp_id = -1
    cell_counter = 0
    ocomp_dict = {}
    obin_dict = {}
    oposition_dict = {}

    min_bin_id = schematic.first_bin
    max_bin_id = schematic.last_bin

    for ic, component in enumerate(schematic.components):
        ocomp = ontology.Component(ic + 1)
        ocomp.ns = zoom_level.ns_term() + '/'
        zoom_level.components.append(ocomp)

        # save the sequence 1-2-3-..-n as a bi-directed list
        if prev_comp_id in ocomp_dict:
            prev_comp = ocomp_dict[prev_comp_id]
            ocomp.reverse_component_edge = prev_comp.ns_term()
            prev_comp.forward_component_edge = ocomp.ns_term()

        ocomp_dict[ic] = ocomp
        prev_comp_id = ic

        obin_ns = ocomp.ns_term() + '/'
        obin_tmp = ontology.Bin()
        obin_tmp.ns = obin_ns

        # bins
        for bins in component.matrix:
            for bin in bins[1][1]:  # follow the compressed format
                if bin:
                    cur_bin_id = bin.bin_id
                    obin = ontology.Bin()
                    obin.ns = obin_ns
                    obin.bin_rank = cur_bin_id
                    obin_dict[cur_bin_id] = obin

                    if (cur_bin_id > min_bin_id):
                        obin_tmp.bin_rank = cur_bin_id - 1
                        obin.reverse_bin_edge = obin_tmp.ns_term()  # string value
                    if (cur_bin_id < max_bin_id):
                        obin_tmp.bin_rank = cur_bin_id + 1
                        obin.forward_bin_edge = obin_tmp.ns_term()  # string value

                    ocomp.bins.append(obin)

                    cell_counter = cell_counter + 1
                    ocell = ontology.Cell()
                    ocell.id = cell_counter
                    ocell.path_id = schematic.path_names[bin.path_id]  # saved in the populate_component_matrix
                    ocell.inversion_percent = bin.inversion
                    ocell.position_percent = bin.position

                    # todo: are begin,end the real bin_ids or the compressed ones? a sparse list sense
                    cell_ns = URIRef("{0}/".format(ocell.path_id))
                    for be in range(0, len(bin.nucleotide_ranges), 2):
                        begin, end = bin.nucleotide_ranges[be], bin.nucleotide_ranges[be + 1]
                        real_begin = begin if begin else end
                        real_end = end if end else begin

                        oregion = ontology.Region()
                        oregion.begin = real_begin
                        oregion.end = real_end
                        ocell.cell_region.append(oregion)

                        path = schematic.path_names[bin.path_id]
                        oposition_begin = ontology.Position(real_begin, begin < end, path, cell_ns)
                        oposition_end = ontology.Position(real_end, begin < end, path, cell_ns)
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
                else:
                    LOGGER.info(f"No upstream {link.upstream}")

                if link.downstream in obin_dict:
                    to_bin = obin_dict[link.downstream]
                    olink.arrival = to_bin.ns_term()
                else:
                    LOGGER.info(f"No downstream {link.downstream}")

                olink.paths = [schematic.path_names[k] for k in link.participants]
                olink.linkZoomLevel = zoom_level.ns_term()
                zoom_level.links.append(olink)

    g = Graph()
    vg = Namespace('http://biohackathon.org/resource/vg#')
    faldo = Namespace('http://biohackathon.org/resource/faldo#')
    g.bind('vg', vg)
    g.bind('faldo', faldo)

    # here the magic happens
    zoom_level.add_to_graph(g, vg, faldo)
    for oposition in oposition_dict.values():
        oposition.add_to_graph(g, vg, faldo)
    for path in schematic.path_names:
        ontology.Path(path).add_to_graph(g, vg, faldo)

    # format='nt' works 10x faster than 'turtle'; the produced files are 3x bigger
    # do not forget to compress them afterwards - gives 15x diskspace reduction
    g.serialize(destination=ontology_output_path, format='nt')
    with open(ontology_output_path, 'rb') as fin:
        with gzip.open(ontology_output_path + '.gz', 'wb') as fout:
            shutil.copyfileobj(fin, fout)
            fout.close()
    os.remove(ontology_output_path)


@dataclass
class PangenomeSchematic:
    json_version: int
    bin_width: int
    first_bin: int
    last_bin: int
    includes_connectors: bool
    components: List[Component]
    path_names: List[str]
    total_nr_files: int
    pangenome_length: int
    file_dict = {}  # dict for bin2file.json

    def json_dump(self):
        def dumper(obj):
            if isinstance(obj, Bin):  # should be in Bin class def
                flat_ranges = obj.nucleotide_ranges
                ranges = []
                for i in range(0, len(flat_ranges), 2):
                    ranges.append([flat_ranges[i], flat_ranges[i+1]])
                return [obj.coverage, obj.inversion, ranges]
            if isinstance(obj, LinkColumn):
                return {'upstream': obj.upstream, 'downstream': obj.downstream, 'participants': obj.participants.tolist()}
            if isinstance(obj, set):
                return list(obj)
            try:
                return obj.__dict__
            except:
                return obj
        # whitespace and [] takes up the majority of the file size
        return json.dumps(self, default=dumper, indent=None, separators=(',', ':\n'))

    def n_links(self, no_adjacent_links):
        columns = sum([len(x.arrivals) + len(x.departures) for x in self.components])
        if not no_adjacent_links:
            columns -= len(self.components)  # -1 for each component
        return columns
    # (len(self.departures)-1) because last departure is adjacent connectors

    def n_components(self):
        return len(self.components)

    def update_first_last_bin(self):
        self.first_bin = 1  # these have not been properly initialized
        self.last_bin = self.components[-1].last_bin

    def split_and_write(self, cells_per_file, folder, fasta : Contig, no_adjacent_links, ontology_folder, parallel):
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
                                               these_comp[-1].last_bin, self.includes_connectors,
                                               these_comp, self.path_names,
                                               self.total_nr_files, self.pangenome_length)
                schematic.filename = self.filename(i)  # save for consistency IMPORTANT

                chunk_summary = {"file": schematic.filename,
                                 "first_bin": schematic.first_bin,
                                 "last_bin": schematic.last_bin,
                                 # used to calculate amount of padding for x values
                                 "component_count": schematic.n_components(),
                                 # extra columns beyond last_bin - first_bin
                                 "link_count": schematic.n_links(no_adjacent_links)}
                if schematic.bin_width == 1:
                    chunk_summary["fasta"] = self.fasta_filename(i)
                bin2file_mapping.append(chunk_summary)

                p = folder.joinpath(schematic.filename)
                with p.open('w') as fpgh9:
                    fpgh9.write(schematic.json_dump())

                if fasta is not None and schematic.bin_width == 1:
                    x = schematic.bin_width
                    if schematic.first_bin == 0:
                        schematic.first_bin = 1  # don't include bin 0 special
                    fa_first = (schematic.first_bin-1) * x  # 1 indexed
                    fa_last = schematic.last_bin * x # inclusive end
                    header = f"first_bin: {schematic.first_bin} " + f"last_bin: {schematic.last_bin}"
                    chunk = [Contig(header, fasta.seq[fa_first:fa_last])]
                    c = folder.joinpath(schematic.fasta_filename(i))
                    write_contigs_to_file(c, chunk)

        if ontology_folder:
            # generator for the pairs (Schematics, path) - will be instantiated on the item access!
            prepared_schematics = ( (PangenomeSchematic(JSON_VERSION, self.bin_width, self.components[cut:cut_points[i + 1]][0].first_bin,
                                                   self.components[cut:cut_points[i + 1]][-1].last_bin, self.includes_connectors,
                                                   self.components[cut:cut_points[i + 1]], self.path_names,
                                                   self.total_nr_files, self.pangenome_length),
                                     str(ontology_folder.joinpath(self.ttl_filename(i))) )
                                   for i, cut in enumerate(cut_points[:-1]) if self.components[cut:cut_points[i + 1]] )

            if parallel:
                results = parallel(delayed(write_rdf)(sch, path) for (sch, path) in prepared_schematics)
            else:
                for (sch, path) in prepared_schematics:
                    write_rdf(sch, path)

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
        avg_paths = mean([len(x.occupants) for x in samples])

        return avg_paths

    def pad_file_nr(self, file_nr):
        return str(file_nr).zfill(int(math.log10(self.total_nr_files) + 1))

    def filename(self, nth_file):
        return f'chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.schematic.json'

    def fasta_filename(self, nth_file):
        return f'seq_chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.fa'

    def ttl_filename(self, nth_file):
        return f'seq_chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.nt'

    def write_index_file(self, folder, bin2file_mapping):

        file_contents = {'bin_width': self.bin_width,
                         'last_bin': self.last_bin,
                         'files': bin2file_mapping}
        self.file_dict[f'{self.bin_width}'] = file_contents

        master_index_file = folder.parent.joinpath(f'bin2file.json')
        master_contents = {'json_version': JSON_VERSION,
                           'pangenome_length': self.pangenome_length,
                           'includes_connectors': self.includes_connectors,
                           'zoom_levels': OrderedDict(sorted(self.file_dict.items(), reverse=True))}
        with master_index_file.open('w') as f:
            f.write(json.dumps(master_contents, indent=4))
            print("Saved file2bin mapping to", master_index_file)

    def prerender(self):
        """Calculates X coordinates and summary statistics for all Components"""
        x = 0
        compressedX = 0
        for component in self.components:
            component.x = x
            component.compressedX = compressedX
            # component.first_bin=0 does not take up rendering space, the next component is 0
            x = component.next_x_coord(self.includes_connectors) if component.first_bin else 0
            compressedX = component.next_compressedX(self.includes_connectors) if component.first_bin else 0
        self.update_first_last_bin()

