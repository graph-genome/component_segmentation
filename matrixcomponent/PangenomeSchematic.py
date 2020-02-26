import json
import math
from bisect import bisect

from math import ceil
from typing import List

from dataclasses import dataclass

from matrixcomponent import JSON_VERSION
from matrixcomponent.matrix import Component, Bin


@dataclass
class PangenomeSchematic:
    json_version: int
    bin_width: int
    first_bin: int
    last_bin: int
    components: List[Component]
    path_names: List[str]
    total_nr_files: int

    def json_dump(self):
        def dumper(obj):
            if isinstance(obj, Bin):  # should be in Bin class def
                return [obj.coverage, obj.inversion, obj.first_nucleotide, obj.last_nucleotide]
            if isinstance(obj, set):
                return list(obj)
            try:
                return obj.__dict__
            except:
                return obj

        return json.dumps(self, default=dumper, indent=2)

    def update_first_last_bin(self):
        self.first_bin = 1  # these have not been properly initialized
        self.last_bin = self.components[-1].last_bin

    def split(self, cells_per_file):
        """Splits one Schematic into multiple files with their own
        unique first and last_bin based on the volume of data desired per
        file specified by cells_per_file.  """
        bins_per_file = ceil(cells_per_file / len(self.path_names))
        # bins_per_schematic = bins_per_row * self.bin_width
        self.update_first_last_bin()
        self.total_nr_files = ceil(self.last_bin / bins_per_file)
        partitions = []
        column_counts = []
        rolling_sum = 0
        # accounting for SVs in column count (cells per file) makes file size much more stable
        for c in self.components:
            rolling_sum += c.column_count()
            column_counts.append(rolling_sum)

        cut_points = [0]
        # variables cut and end_cut are componentIDs
        # binIDs are in components.{first,last}_bin

        prev_point = 0
        for start_bin in range(0, column_counts[-1] + bins_per_file, bins_per_file):
            cut = bisect(column_counts, start_bin + bins_per_file)
            cut_points.append(max(prev_point + 1, cut))
            prev_point = cut
        cut_points.append(len(self.components))  # don't chop of dangling end

        bin2file_mapping = []
        for i, cut in enumerate(cut_points[:-1]):
            end_cut = cut_points[i + 1]
            these_comp = self.components[cut:end_cut]
            if these_comp:  # when we're out of data because the last component is 1 wide
                partitions.append(
                    PangenomeSchematic(JSON_VERSION,
                                       self.bin_width,
                                       these_comp[0].first_bin,
                                       these_comp[-1].last_bin,
                                       these_comp, self.path_names, self.total_nr_files))
                bin2file_mapping.append({"file": self.filename(i),
                                         "first_bin": these_comp[0].first_bin,
                                         "last_bin": these_comp[-1].last_bin})
        return partitions, bin2file_mapping

    def pad_file_nr(self, file_nr):
        return str(file_nr).zfill(int(math.log10(self.total_nr_files)))

    def filename(self, nth_file):
        return f'chunk{self.pad_file_nr(nth_file)}_bin{self.bin_width}.schematic.json'

    def write_index_file(self, folder, bin2file_mapping):
        """Also write the file2bin mapping into a master json file
        eventually, this could have one list per bin size,
        with all bin size integrated into the same folder"""
        index_file = folder.joinpath(f'bin2file.json')
        file_contents = {'bin_width': self.bin_width,
                         'json_version': JSON_VERSION,
                         'last_bin': self.last_bin,
                         'files': bin2file_mapping}
        with index_file.open('w') as out:
            out.write(json.dumps(file_contents, indent=4))
            print("Saved file2bin mapping to", index_file)
