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

from matrixcomponent.matrix import Path, PangenomeSchematic, Component, LinkColumn
import os
import logging
import argparse
import matrixcomponent

import matrixcomponent.JSONparser as JSONparser

LOGGER = logging.getLogger(__name__)
"""logging.Logger: The logger for this module"""


def populate_component_occupancy(matrix: List[Path], schematic: PangenomeSchematic):
    for component in schematic.components:
        occupants = [False] * len(schematic.path_names)
        # are matrix paths in the same order as schematic.path_names?
        for i, path in enumerate(matrix):
            relevant = [bin for bin in path.bins if bin.bin_id >= component.first_bin and bin.bin_id <= component.last_bin]  # very costly loop
            # TODO: When we move to breakpoints, this will be
            # TODO: occupants[component.region.pathnames.index(path.name)]
            occupants[i] = any([bin.coverage > 0.1 for bin in relevant])
        component.occupants = occupants  # side effect instead of return
    print("Populated Occupancy per component per path.")


def segment_matrix(matrix: List[Path]) -> PangenomeSchematic:
    print(f"Starting Segmentation process on {len(matrix)} Paths.")
    schematic = PangenomeSchematic(100000, 1, 140*1000000, # FIXME: parse bin_size from filename
                                   [], [p.name for p in matrix], [])
    incoming, outgoing, dividers = find_dividers(matrix)
    start_pos = 0
    for valid_start in sorted(list(dividers)):
        if valid_start != 0:
            current = Component(start_pos, valid_start - 1)
            # current.active_members = 1
            schematic.components.append(current)
        start_pos = valid_start
    print(f"Created {len(schematic.components)} components")

    #populate Component occupancy per Path
    populate_component_occupancy(matrix, schematic)

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
    print(uniq_links)

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
    LOGGER.info("starting...\n")
    Paths = JSONparser.parse(args.json_file)
    schematic = segment_matrix(Paths)
    p = osPath(args.json_file).with_suffix('.schematic.json')
    with p.open('w') as fpgh9:
        fpgh9.write(schematic.json_dump())
    print("Saved results to", p)


if __name__ == '__main__':

    main()

