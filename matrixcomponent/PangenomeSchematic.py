import json
from collections import OrderedDict
from statistics import mean

import math
from bisect import bisect

from math import ceil
from typing import List

from dataclasses import dataclass

from matrixcomponent import JSON_VERSION
from matrixcomponent.matrix import Component, Bin, LinkColumn

from DNASkittleUtils.Contigs import Contig, write_contigs_to_file

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
                # todo: get rid of this once the JS side can work with sparse containers
                if self.json_version <= 14:
                    bools = [False] * len(self.path_names)
                    for i in obj.participants:
                        bools[i] = True
                    return {'upstream':obj.upstream, 'downstream':obj.downstream, 'participants':bools}
                else:
                    return {'upstream':obj.upstream, 'downstream':obj.downstream, 'participants':obj.participants.tolist()}
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

    def split_and_write(self, cells_per_file, folder, fasta : Contig, no_adjacent_links):
        # todo: get rid of this once the JS side can work with sparse containers
        if self.json_version <= 14:
            empty = []
            for comp in self.components:
                bools = [False] * len(self.path_names)
                for i in comp.occupants:
                    bools[i] = True
                comp.occupants = bools

                matrix = [empty] * len(self.path_names)
                fb, lb = comp.first_bin, comp.last_bin
                for item in comp.matrix:
                    padded = [empty] * (lb - fb + 1)
                    sliced = item[1]
                    for id, val in zip(sliced[0], sliced[1]):
                        padded[id] = val
                    matrix[item[0]] = padded

                comp.matrix = matrix


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
                                 "x": these_comp[0].x,
                                 "compressedX": these_comp[0].compressedX}
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

        # todo: get rid of this once the JS side can work with sparse containers
        if self.json_version <= 14:
            avg_paths = mean([sum(x.occupants) for x in samples])
        else:
            avg_paths = mean([len(x.occupants) for x in samples])

        return avg_paths

    def pad_file_nr(self, file_nr):
        return str(file_nr).zfill(int(math.log10(self.total_nr_files) + 1))

    def filename(self, nth_file):
        return f'chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.schematic.json'

    def fasta_filename(self, nth_file):
        return f'seq_chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.fa'

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

