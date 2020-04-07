#!/usr/bin/env python3
"""
Tasks
Get example output for tests - Ted
Parser for ODGI Bin file format - Ted
Component Segmentation Detection - Josiah and Joerg
  Python memory object model - Josiah
Output format
"""
from typing import List, Tuple, Set, Dict
from pathlib import Path as osPath, PurePath
from nested_dict import nested_dict
from datetime import datetime
from sortedcontainers import SortedDict
from DNASkittleUtils.Contigs import Contig, read_contigs, write_contigs_to_file

from matrixcomponent.matrix import Path, Component, LinkColumn, Bin
from matrixcomponent.PangenomeSchematic import PangenomeSchematic
import matrixcomponent.utils as utils

import os
import logging
import argparse
import matrixcomponent

import matrixcomponent.JSONparser as JSONparser

import numpy as np
import pandas as pd

MAX_COMPONENT_SIZE = 100  # automatic calculation from cells_per_file did not go well
LOGGER = logging.getLogger(__name__)
"""logging.Logger: The logger for this module"""


def populate_component_occupancy(schematic: PangenomeSchematic):
    for component in schematic.components:
        # are matrix paths in the same order as schematic.path_names?
        # side effect instead of return
        component.occupants = [any([bin.coverage > 0.1 for bin in bins if bin])
                               for bins in component.matrix]
    print("Populated Occupancy per component per path.")

def populate_component_matrix(paths: List[Path], schematic: PangenomeSchematic):
    # the loops are 1) paths, and then 2) schematic.components
    # paths are in the same order as schematic.path_names
    for i, path in enumerate(paths):
        sorted_bins = SortedDict((bin.bin_id, bin) for bin in path.bins)
        values = list(sorted_bins.values())
        for component in schematic.components:
            from_id = sorted_bins.bisect_left (component.first_bin)
            to_id   = sorted_bins.bisect_right(component.last_bin)
            relevant = values[from_id:to_id]
            padded = []
            if relevant:
                padded = [[]] * (component.last_bin - component.first_bin + 1)
                for bin in relevant:
                    padded[bin.bin_id - component.first_bin] =  \
                        Bin(bin.coverage, bin.inversion_rate, bin.first_nucleotide, bin.last_nucleotide)
            component.matrix.append(padded)  # ensure there's always 1 entry for each path
    print("Populated Matrix per component per path.")
    populate_component_occupancy(schematic)


def segment_matrix(matrix: List[Path], bin_width, cells_per_file, pangenome_length) -> PangenomeSchematic:
    from matrixcomponent import JSON_VERSION
    print(f"Starting Segmentation process on {len(matrix)} Paths.")
    schematic = PangenomeSchematic(JSON_VERSION,
                                   bin_width,
                                   1,
                                   1,
                                   [], [p.name for p in matrix], 1, pangenome_length)
    connections, dividers = dividers_with_max_size(matrix, cells_per_file)
    start_pos = 0
    for valid_start in dividers:
        if valid_start != 0:
            current = Component(start_pos, valid_start - 1)
            # current.active_members = 1
            schematic.components.append(current)
        start_pos = valid_start
    print(f"Created {len(schematic.components)} components")

    # populate Component occupancy per Path
    populate_component_matrix(matrix, schematic)

    incoming = nested_dict(2, set)
    outgoing = nested_dict(2, set)
    for f, t, p in connections.itertuples(index=False):
        incoming[t][f].add(p)
        outgoing[f][t].add(p)

    # populate all link columns onto schematic
    nLinkColumns = 0
    for component in schematic.components:
        # TODO: order columns based on traversal patterns,
        # TODO: insert additional columns for higher copy number
        for origin_pos, participants in incoming[component.first_bin].items():
            phase_dots = [indiv in participants for indiv in schematic.path_names]
            entering = LinkColumn(origin_pos,
                                  component.first_bin,
                                  participants=phase_dots)
            component.arrivals.append(entering)
            nLinkColumns += 1
        for arriving_pos, participants in outgoing[component.last_bin].items():
            # phase_dots depends on row ordering of path names, optimized for display
            phase_dots = [indiv in participants for indiv in schematic.path_names]
            leaving = LinkColumn(component.last_bin,
                                 arriving_pos,
                                 participants=phase_dots)
            component.departures.append(leaving)
            nLinkColumns += 1

    for i in range(len(schematic.components)-1):
        component, next_component = schematic.components[i],schematic.components[i+1]
        add_adjacent_connector_column(component, next_component, schematic)

    print(f"Created {nLinkColumns} LinkColumns")

    return schematic


def dividers_with_max_size(matrix: List[Path], cells_per_file: int):
    """Adds in additional dividers to ensure very large components are split into
    multiple components with no Links."""
    connections, dividers = find_dividers(matrix)
    # estimate number of paths, x10 because most paths are empty
    dividers_extended = []
    prev = 0
    for div in sorted(list(dividers)):
        gap_size = div - prev
        if gap_size > MAX_COMPONENT_SIZE:
            for i in range(prev + MAX_COMPONENT_SIZE, div, MAX_COMPONENT_SIZE):
                dividers_extended.append(i)  # add a series of dividers spaced ^ apart
        prev = div
        dividers_extended.append(div)

    return connections, dividers_extended


def add_adjacent_connector_column(component, next_component, schematic):
    """The last Departure LinkColumn is to the adjacent component
    Use logic to decide on which rows need adjacent connectors
    Start with the easy subtractive case of occupancy - departures and move to more complex,
    multiple copy cases."""
    adjacents = []
    for row in range(len(schematic.path_names)):
        connection_exists = False
        if component.occupants[row] and next_component.occupants[row]:  # occupant present
            # n_arrivals = sum([column.participants[row] for column in component.arrivals])
            departed = sum([column.participants[row] for column in component.departures])
            # connection_exists = n_arrivals + 1 > departed
            connection_exists = not departed  # didn't depart
        adjacents.append(connection_exists)
    component.departures.append(LinkColumn(  # LinkColumn for adjacents
        component.last_bin,
        component.last_bin + 1,
        participants=adjacents))


def find_dividers(matrix: List[Path]) -> Tuple[pd.DataFrame, Set[int]]:
    max_bin = 1

    self_loops = []  # track self loops just in case component gets cut in half
    connection_dfs = []  # pandas dataframe with columns (from, to, path [name])

    for i, path in enumerate(matrix):
        links = path.links
        bin_ids = np.array([b.bin_id for b in path.bins])
        bin_ids.sort()

        if bin_ids.size > 0:
            max_bin = max(max_bin, int(bin_ids[-1]))

        # we don't want these to become dividers
        boundary_mask = utils.path_boundaries(links)
        self_loops_mask = utils.self_loops(links)

        if np.any(self_loops_mask):
            self_loops.append(links[self_loops_mask])

        links = links[~(boundary_mask | self_loops_mask)]

        path_dividers_mask = utils.path_dividers(links, bin_ids)
        path_dividers = links[path_dividers_mask]
        if path_dividers.size == 0:
            continue

        df = pd.DataFrame.from_dict({
            'from': path_dividers[:, 0],  # aka upstream
            'to': path_dividers[:, 1],    # aka downstream
            'path': path.name
        })

        connection_dfs.append(df)

        # <old comments applicable to each divider>
        #
        # if (upstream + 1) in leaving.keys() :
        #     print(f"Found inherited rearrangement {upstream+1}")
        #
        # TODO: insert prevarications about exact position
        # Divider should be somewhere in here
        # Tolerable range?
        # Stack up others using the same LinkColumn

    df = pd.concat(connection_dfs)
    df.drop_duplicates(inplace=True)

    # all start positions of components
    dividers = set(df["from"] + 1) | set(df["to"])
    dividers.add(1)

    print(f"Largest bin_id was {max_bin}\n"
          f"Found {len(dividers)} dividers.")
    dividers.add(max_bin + 1)  # end of pangenome

    if self_loops:
        n_self_loops = np.unique(np.concatenate(self_loops), axis=0).shape[0]
        print(f"Eliminated {n_self_loops} self-loops")

    n_uniq_links = df[df["to"] >= df["from"]]\
        .drop(columns=["path"]).drop_duplicates().shape[0]

    n_links = sum([len(p.links) for p in matrix])
    print(f"Input has {n_links} listed Links.  "
          f"Segmentation eliminated {(1-n_uniq_links/n_links)*100}% of them.")
    print(f"Found {n_uniq_links} unique links")

    return df, dividers


def setup_logging():
    """Setup the logging, add a log file"""
    log_name = osPath(args.json_file).with_suffix('.log')
    if args.output_folder:
        log_name = osPath(args.output_folder).joinpath('log')
        os.makedirs(args.output_folder, exist_ok=True)
    t = datetime.now()
    timestr = f"{t.year}{t.month:02}{t.day:02}-{t.hour:02}-{t.minute:02}-{t.second:02}"
    log_name = str(log_name) + '.' + timestr + '.log'

    handler = logging.FileHandler(os.path.join(log_name))
    handler.setLevel(args.log_level)
    handler.setFormatter(logging.Formatter(matrixcomponent.LOGGING_FORMAT_STR,
                                           datefmt=matrixcomponent.LOGGING_DATE_FORMAT))
    logging.getLogger().addHandler(handler)


# Helper class to allow multi-line help messages for argparse user parameters:
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)


def write_json_files(json_file, schematic: PangenomeSchematic):
    partitions, bin2file_mapping = schematic.split(args.cells_per_file)

    folder = osPath(json_file).parent
    if args.output_folder:
        folder = osPath(args.output_folder)
    os.makedirs(folder, exist_ok=True)  # make directory for all files

    for part in partitions:
        p = folder.joinpath(part.filename)
        with p.open('w') as fpgh9:
            fpgh9.write(part.json_dump())
        print("Saved results to", p)

    schematic.write_index_file(folder, bin2file_mapping)


def write_fasta_files(odgi_fasta, output_folder: Path, schematic: PangenomeSchematic):
    partitions, bin2file_mapping = schematic.split(args.cells_per_file)
    fasta = read_contigs(odgi_fasta)[0]
    for part in partitions:
        x = part.bin_width
        fa_first, fa_last = (part.first_bin * x), ((part.last_bin + 1) * x)
        header = f"first_bin: {part.first_bin} " + f"last_bin: {part.last_bin}"
        chunk = [Contig(header, fasta.seq[fa_first:fa_last])]
        c = PurePath(output_folder).joinpath(part.fasta_filename)
        write_contigs_to_file(c, chunk)


def get_arguments():
    """Create the command line interface and return the command line arguments

    Returns
    -------
    Namespace
        The command line arguments

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Example Command:\n"
                    "--json-file=data/run1.B1phi1.i1.seqwish.w100.json --cells-per-file=5000 --fasta=data/run1.B1phi1.i1.seqwish.fasta")

    parser.add_argument('-j', '--json-file',
                            dest='json_file',
                            required=True,
                            help='input JSON file')

    parser.add_argument('-f', '--fasta',
                        dest='fasta',
                        help='Optional: Fasta file containing the pangenome sequence generated by '
                             'odgi for this Graph.')

    parser.add_argument('-o', '--out-folder',
                        dest='output_folder',
                        help='output folder')

    parser.add_argument('-c', '--cells-per-file',
                        dest='cells_per_file',
                        default=5000,
                        type=int,
                        help='Tip: Adjust this number to get chunk files output close to 2MB. '
                             'Number of cells per file (#bins per file = #cells / #paths)')

    parser.add_argument('-l', '--log-level',
                        default='DEBUG',
                        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR'),
                        help='level of logging verbosity. DEBUG is most verbose')

    args = parser.parse_args()
    if not args.output_folder:
        # directory with the same name as the json
        args.output_folder = osPath(args.json_file).parent.joinpath(osPath(args.json_file).stem)
    else:
        args.output_folder = osPath(args.output_folder)
    os.makedirs(args.output_folder, exist_ok=True)

    return args


def main():
    global args
    args = get_arguments()
    setup_logging()
    LOGGER.info(f'reading {osPath(args.json_file)}...\n')
    paths, pangenome_length, bin_width = JSONparser.parse(args.json_file)
    schematic = segment_matrix(paths, bin_width, args.cells_per_file, pangenome_length)
    del paths
    write_json_files(args.output_folder, schematic)
    if args.fasta:  # optional
        write_fasta_files(args.fasta, args.output_folder, schematic)


if __name__ == '__main__':
    main()
#--json-file=data/run1.B1phi1.i1.seqwish.w100.json --cells-per-file=5000
# --fasta=data/run1.B1phi1.i1.seqwish.fasta
