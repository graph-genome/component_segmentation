import json
from collections import OrderedDict
from statistics import mean

import math
from bisect import bisect

from math import ceil
from typing import List

from dataclasses import dataclass

from matrixcomponent import JSON_VERSION, ontology
from matrixcomponent.matrix import Component, Bin, LinkColumn

from rdflib import Graph, Namespace, URIRef

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

    def update_first_last_bin(self):
        self.first_bin = 1  # these have not been properly initialized
        self.last_bin = self.components[-1].last_bin

    def split_and_write(self, cells_per_file, folder, fasta: Contig, ontology_folder):
        """Splits one Schematic into multiple files with their own
        unique first and last_bin based on the volume of data desired per
        file specified by cells_per_file.  """
        avg_paths = self.lazy_average_occupants()
        bins_per_file = ceil(cells_per_file / avg_paths)
        column_counts = self.rolling_sum_column_count()
        cut_points = self.find_cut_points_in_file_split(bins_per_file, column_counts)

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

                if fasta is not None:
                    schematic.fasta_filename = self.fasta_filename(i)

                p = folder.joinpath(schematic.filename)
                with p.open('w') as fpgh9:
                    fpgh9.write(schematic.json_dump())

                if fasta is not None:
                    x = schematic.bin_width
                    fa_first, fa_last = (schematic.first_bin * x), ((schematic.last_bin + 1) * x)
                    header = f"first_bin: {schematic.first_bin} " + f"last_bin: {schematic.last_bin}"
                    chunk = [Contig(header, fasta.seq[fa_first:fa_last])]
                    c = folder.joinpath(schematic.fasta_filename)
                    write_contigs_to_file(c, chunk)

                if schematic.bin_width == 1:
                    bin2file_mapping.append({"file": schematic.filename,
                                             "fasta": schematic.fasta_filename,
                                             "first_bin": schematic.first_bin,
                                             "last_bin": schematic.last_bin})
                else:
                    bin2file_mapping.append({"file": schematic.filename,
                                         "first_bin": schematic.first_bin,
                                         "last_bin": schematic.last_bin})

                if ontology_folder:
                    # going from top to bottom
                    zoom_level = ontology.ZoomLevel()
                    zoom_level.zoom_factor = schematic.bin_width
                    zoom_level.ns = URIRef('pg/')

                    cell_counter = 0
                    obin_dict = {}
                    for ic, component in enumerate(schematic.components):
                        ocomp = ontology.Component(ic + 1)
                        zoom_level.components.append(ocomp)

                        # bins
                        for bins in component.matrix:
                            for bin in bins:
                                if bin:
                                    obin = ontology.Bin()
                                    obin.bin_rank = bin.bin_id
                                    obin.path_id  = bin.path_id # saved in the populate_component_matrix
                                    obin_dict[bin.bin_id] = obin
                                    ocomp.bins.append(obin)

                                    cell_counter = cell_counter + 1
                                    ocell = ontology.Cell()
                                    ocell.id = cell_counter
                                    ocell.inversion_percent = bin.inversion
                                    ocell.position_percent = bin.coverage

                                    for [begin, end] in bin.nucleotide_ranges:
                                        oregion = ontology.Region()
                                        oregion.begin = begin
                                        oregion.end = end
                                        ocell.cell_region.append(oregion)

                                    obin.cells.append(ocell)

                    # links between components and their bins
                    for component, ocomp in zip(schematic.components, zoom_level.components):
                        # links: arrivals
                        link_counter = 0
                        for link in component.arrivals:
                            link_counter = link_counter + 1
                            olink = ontology.Link()
                            olink.id = link_counter

                            if link.upstream in obin_dict:
                                olink.departure = obin_dict[link.upstream]

                            if link.downstream in obin_dict:
                                olink.arrival   = obin_dict[link.downstream]

                            olink.paths = (link.participants + 1).tolist()
                            ocomp.links.append(olink)

                        for link in component.departures:
                            link_counter = link_counter + 1
                            olink = ontology.Link()
                            olink.id = link_counter

                            if link.upstream in obin_dict:
                                olink.departure = obin_dict[link.upstream]

                            if link.downstream in obin_dict:
                                olink.arrival = obin_dict[link.downstream]

                            olink.paths = (link.participants + 1).tolist()
                            ocomp.links.append(olink)


                    g = Graph()
                    vg = Namespace('http://biohackathon.org/resource/vg#')
                    faldo = Namespace('http://biohackathon.org/resource/faldo#')
                    g.bind('vg', vg)
                    g.bind('faldo', faldo)

                    zoom_level.add_to_graph(g, vg, faldo)  # here the magic happens

                    p = ontology_folder.joinpath(schematic.ttl_filename(i))
                    g.serialize(destination=str(p), format='turtle', encoding='utf-8')


        return bin2file_mapping

    def find_cut_points_in_file_split(self, bins_per_file, column_counts):
        """Use binary search bisect to find the component boundaries closest to
        target breakpoints for splitting files."""
        cut_points, prev_point = [0], 0
        for start_bin in range(0, column_counts[-1], bins_per_file):
            cut = bisect(column_counts, start_bin + bins_per_file)
            cut_points.append(max(prev_point + 1, cut))
            prev_point = cut
        cut_points.append(len(self.components))  # don't chop of dangling end
        self.total_nr_files = len(cut_points) - 1
        return cut_points

    def rolling_sum_column_count(self):
        """accounting for SVs in column count (cells per file) makes file size much more stable"""
        self.update_first_last_bin()
        column_counts, rolling_sum = [], 0
        for c in self.components:
            rolling_sum += c.column_count()
            column_counts.append(rolling_sum)
        return column_counts

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
