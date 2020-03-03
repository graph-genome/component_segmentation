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
from pathlib import Path as osPath
from nested_dict import nested_dict

from matrixcomponent.matrix import Path, Component, LinkColumn, Bin
from matrixcomponent.PangenomeSchematic import PangenomeSchematic
import os
import logging
import argparse
import matrixcomponent

import matrixcomponent.JSONparser as JSONparser

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
    for component in schematic.components:
        # paths paths are in the same order as schematic.path_names
        for i, path in enumerate(paths):
            relevant = [bin for bin in path.bins if
                        component.first_bin <= bin.bin_id <= component.last_bin]  # very costly loop
            padded = []
            if relevant:
                padded = [[]] * (component.last_bin - component.first_bin + 1)
                for bin in relevant:
                    padded[bin.bin_id - component.first_bin] =  \
                        Bin(bin.coverage, bin.inversion_rate, bin.first_nucleotide, bin.last_nucleotide)
            component.matrix.append(padded)  # ensure there's always 1 entry for each path
    print("Populated Matrix per component per path.")
    populate_component_occupancy(schematic)


def segment_matrix(matrix: List[Path], bin_width, cells_per_file) -> PangenomeSchematic:
    from matrixcomponent import JSON_VERSION
    print(f"Starting Segmentation process on {len(matrix)} Paths.")
    schematic = PangenomeSchematic(JSON_VERSION,
                                   bin_width,
                                   1,
                                   1,
                                   [], [p.name for p in matrix], 1)
    incoming, outgoing, dividers = dividers_with_max_size(matrix, cells_per_file)
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
    incoming, outgoing, dividers = find_dividers(matrix)
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

    return incoming, outgoing, dividers_extended


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


def find_dividers(matrix: List[Path]) -> Tuple[Dict[int, Dict[int, set]],
                                               Dict[int, Dict[int, set]], Set[int]]:
    max_bin = 1
    leaving = nested_dict(2, set)  # containing the set of participating Paths on that link column
    entering = nested_dict(2, set)  # list of indices of new components
    dividers = {1}  # all start positions of components, start with st
    copy_arrivals = set()  # track self loops just in case component gets cut in half
    uniq_links = set()
    for i, path in enumerate(matrix):
        print(f"Segmenting {len(path.bins)}", '{0:.1%}'.format(i / len(matrix)))
        max_bin = max(max_bin, max(path._bin_set))
        for link in path.links:  # Links are generated by odgi based
            upstream, downstream = link.upstream, link.downstream
            # Links to 0 Bin indicate the beginning or end of a path.  0 Bin has no sequence
            if 0 in [upstream, downstream]:
                continue  # ignore links to the telomere.
            if upstream == downstream:
                copy_arrivals.add(upstream)
                continue  # we don't want these to become dividers
            # Is the gap range anywhere else in this individual?
            # What if downstream < upstream?
            divider_verified = downstream < upstream
            if not divider_verified:
                missing_range = list(range(upstream + 1, downstream))
                for i in missing_range:
                    if i in path:
                        divider_verified = True
                        uniq_links.add((upstream, downstream))
                        break  # stop as soon as we have confirmation
            if divider_verified:
                # if (upstream + 1) in leaving.keys() :
                #     print(f"Found inherited rearrangement {upstream+1}")
                leaving[upstream][downstream].add(path.name)  # the first position of the new component
                entering[downstream][upstream].add(path.name)
                dividers.add(upstream + 1)
                dividers.add(downstream)
                # TODO: insert prevarications about exact position
                # Divider should be somewhere in here
                # Tolerable range?
                # Stack up others using the same LinkColumn
    print(f"Largest bin_id was {max_bin}\n"
          f"Found {len(dividers)} dividers.")
    dividers.add(max_bin + 1)  # end of pangenome
    print(f"Eliminated {len(copy_arrivals)} self-loops")
    n_links = sum([len(p.links) for p in matrix])
    print(f"Input has {n_links} listed Links.  "
          f"Segmentation eliminated {(1-len(uniq_links)/n_links)*100}% of them.")
    print(f"Found {len(uniq_links)} unique links")

    return entering, leaving, dividers


def discard_useless_links(matrix: List[Path]):
    """https://github.com/vgteam/odgi/issues/48
    Links that simply span a gap in the matrix can be discarded"""
    for path in matrix:
        keep = []
        for link in path.links:  # Links are generated by odgi based
            missing_range = list(range(link.upstream + 1, link.downstream))
            # Is the gap range anywhere else in this individual?
            if any([i in path for i in missing_range if i > 0]):
                keep.append(link)
        path.links = keep  # all other Paths get deleted


def setup_logging(output_dir):
    """Setup the logging, add a log file"""
    log_name = os.path.join(output_dir, 'log')
    handler = logging.FileHandler(log_name)
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
    folder = osPath(json_file).with_suffix('')
    os.makedirs(folder, exist_ok=True)  # make directory for all files
    for part in partitions:
        p = folder.joinpath(part.filename)
        with p.open('w') as fpgh9:
            fpgh9.write(part.json_dump())
        print("Saved results to", p)

    schematic.write_index_file(folder, bin2file_mapping)


def get_arguments():
    """Create the command line interface and return the command line arguments

    Returns
    -------
    Namespace
        The command line arguments

    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="")

    parser.add_argument('-j', '--json-file',
                            dest='json_file',
                            required=True,
                            help='input JSON file')

    parser.add_argument('-o', '--out-folder',
                        dest='output_folder',
                        required=True,
                        help='output folder')

    parser.add_argument('-b', '--bin-width',
                        dest='bin_width',
                        required=True,
                        help='bin size: number of nucleotides per bin')

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

    return args


def main():
    global args
    args = get_arguments()
    setup_logging(args.output_folder)
    LOGGER.info(f'reading {osPath(args.json_file)}...\n')
    paths = JSONparser.parse(args.json_file)
    schematic = segment_matrix(paths, args.bin_width, args.cells_per_file)
    del paths
    write_json_files(args.json_file, schematic)


if __name__ == '__main__':
    main()

# --json-file=data/run1.B1phi1.i1.seqwish.w100.json --out-folder=data/ --bin-width=100 --cells-per-file=5000
# --json-file=data/Athaliana_12_individuals_w100000.json --out-folder=data/ --bin-width=100000 --cells-per-file=10000