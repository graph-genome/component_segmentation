#!/usr/bin/env python3
"""
Tasks
Get example output for tests - Ted
Parser for ODGI Bin file format - Ted
Component Segmentation Detection - Josiah and Joerg
  Python memory object model - Josiah
Output format
"""

from typing import List, Tuple
from pathlib import Path as osPath
import sys
from datetime import datetime
from DNASkittleUtils.Contigs import read_contigs
from joblib import Parallel

from matrixcomponent.matrix import Path, Component, LinkColumn
from matrixcomponent.PangenomeSchematic import PangenomeSchematic
import matrixcomponent.utils as utils

import os
from glob import glob
import logging
import argparse
import matrixcomponent

import matrixcomponent.JSONparser as JSONparser

import numpy as np

MAX_COMPONENT_SIZE = 100  # automatic calculation from cells_per_file did not go well
LOGGER = logging.getLogger(__name__)
"""logging.Logger: The logger for this module"""

# is set to 0 in the very end
# if exited prematurely will indicate an unsuccessful execution
exitcode = 1

def populate_component_matrix(paths: List[Path], schematic: PangenomeSchematic):
    # the loops are 1) paths, and then 2) schematic.components
    # paths are in the same order as schematic.path_names
    first_bins = np.asarray([component.first_bin for component in schematic.components])
    last_bins  = np.asarray([component.last_bin for component in schematic.components])
    comp_array = np.asarray([component for component in schematic.components])

    empty = []

    for p,path in enumerate(paths):
        sorted_bins = path.bins
        keys = np.asarray(sorted_bins.keys())
        values = list(sorted_bins.values())

        # Numpy vectorized binary search
        from_ids = np.searchsorted(keys, first_bins, side='left')
        to_ids   = np.searchsorted(keys, last_bins,  side='right')

        # the case "from+1 == to". Here we can save lots of unnecessary slicing and looping
        mask = from_ids+1 == to_ids
        if np.any(mask):
            comp_filtered = comp_array[mask]
            from_filtered = from_ids[mask]

            # synchron loop over all arrays
            # this case enforces first_bin == last_bin --- comp.matrix[p] has a single element
            for comp, fr in zip(comp_filtered, from_filtered):
                bin = values[fr]
                bin.path_id = p  # save for later
                comp.matrix.append([p, [[0], [bin]]])
                if bin.coverage > 0.1:
                    comp.occupants.add(p)


        # and a general one "from+1 < to"
        mask = from_ids+1 < to_ids
        if np.any(mask):
            comp_filtered = comp_array[mask]
            from_filtered = from_ids[mask]
            to_filtered   = to_ids[mask]

            # synchron loop over all arrays
            for comp,fr,to in zip(comp_filtered,from_filtered,to_filtered):
                fb = comp.first_bin
                sliced = values[fr:to]
                ids = [bin.bin_id - fb for bin in sliced]
                comp.matrix.append([p,[ids, sliced]])
                for bin in values[fr:to]:
                    bin.path_id = p  # save for later
                if any([bin.coverage > 0.1 for bin in sliced]):
                    comp.occupants.add(p)

    LOGGER.info("Populated Matrix and Occupancy per component per path.")


def segment_matrix(matrix: List[Path], bin_width, cells_per_file, pangenome_length, no_adjacent_links, parallel) -> PangenomeSchematic:
    from matrixcomponent import JSON_VERSION
    LOGGER.info(f"Starting Segmentation process on {len(matrix)} Paths.")
    schematic = PangenomeSchematic(JSON_VERSION,
                                   bin_width,
                                   1,
                                   1,
                                   not no_adjacent_links,
                                   [], [p.name for p in matrix], 1, pangenome_length)
    connections, dividers = dividers_with_max_size(matrix, cells_per_file)
    LOGGER.info(f"Created dividers")

    component_by_first_bin = {}
    component_by_last_bin = {}
    start_pos = 0
    for valid_start in dividers:
        if valid_start != 0:
            current = Component(start_pos, valid_start - 1)
            # current.active_members = 1
            schematic.components.append(current)
            component_by_first_bin[start_pos] = current
            component_by_last_bin[valid_start - 1] = current
        start_pos = valid_start
    LOGGER.info(f"Created {len(schematic.components)} components")

    # populate Component occupancy per Path
    populate_component_matrix(matrix, schematic)
    LOGGER.info(f"populated matrix")

    path_indices = connections['path_index']
    connections_from = connections['from']
    connections_to   = connections['to']
    groups = utils.find_groups(connections_from, connections_to)

    for i in range(len(groups) - 1):
        start, end = groups[i], groups[i+1]
        src, dst = int(connections_from[start]), int(connections_to[start]) # important to cast to int()

        link_column = LinkColumn(src, dst, participants=path_indices[start:end])

        src_component = component_by_last_bin.get(src)
        dst_component = component_by_first_bin.get(dst)

        if src_component:
            src_component.departures.append(link_column)

        if dst_component:
            dst_component.arrivals.append(link_column)

    if not no_adjacent_links:
        for i in range(len(schematic.components)-1):
            component, next_component = schematic.components[i],schematic.components[i+1]
            add_adjacent_connector_column(component, next_component, schematic)
        # add special case connectors for the last component in the file
        add_adjacent_connector_column(schematic.components[-1], None, schematic)

    num_link_columns = sum([(len(comp.departures) + len(comp.arrivals)) for comp in schematic.components])
    LOGGER.info(f"Created {num_link_columns} LinkColumns")

    schematic.prerender()

    return schematic


def dividers_with_max_size(matrix: List[Path], cells_per_file: int):
    """Adds in additional dividers to ensure very large components are split into
    multiple components with no Links."""
    connections, dividers = find_dividers(matrix)
    # estimate number of paths, x10 because most paths are empty
    dividers_extended = []
    prev = 0
    for div in dividers:
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

    ids = np.arange(len(schematic.path_names))
    common = component.occupants & next_component.occupants if (component and next_component) else []
    filtered_rows = np.asarray([ids[j] for j in common])
    adjacents = filtered_rows # we take all the filtered IDs if there are no departures

    if len(filtered_rows) > 0 and len(component.departures) > 0: # potentially there's work to do
        ids = np.concatenate([column.participants for column in component.departures])
        isin = np.isin(filtered_rows, ids, invert=True)
        adjacents = filtered_rows[isin]

    # if adjacents.size > 0:  # add linkcolumn as placeholder even when an empty list of participants
    component.departures.append(LinkColumn(  # LinkColumn for adjacents
        component.last_bin,
        component.last_bin + 1,
        participants=np.asarray(adjacents).astype(dtype='int32')))


def find_dividers(matrix: List[Path]) -> Tuple[dict, List[int]]:
    max_bin = max([p.max_bin_id for p in matrix])
    n_links = sum([p.num_links for p in matrix])
    n_remaining_links = sum([len(p.path_dividers) for p in matrix])

    n_self_loops = 0
    self_loops = [p.self_loops for p in matrix if p.self_loops.size > 0]
    if self_loops:
        n_self_loops = np.unique(np.concatenate(self_loops), axis=0).shape[0]

    # length of the binary representation
    path_shift = len(bin(len(matrix))) - 2
    shift = len(bin(max_bin)) - 2

    connection_dfs = []
    for i, path in enumerate(matrix):
        path.n_self_loops = np.array([], dtype='int32') # no need to keep this array - erase it

        path_dividers = path.path_dividers
        if path_dividers.size == 0:
            continue

        arr = np.stack( (path_dividers[:, 0], path_dividers[:, 1], i*np.ones(len(path_dividers)).astype(dtype='int32')) )
        connection_dfs.append(utils.compress_array(arr, shift, path_shift))
        path.path_dividers = np.array([], dtype='int32') # no need to keep this array - erase it

        # <old comments applicable to each divider>
        #
        # if (upstream + 1) in leaving.keys() :
        #     print(f"Found inherited rearrangement {upstream+1}")
        #
        # TODO: insert prevarications about exact position
        # Divider should be somewhere in here
        # Tolerable range?
        # Stack up others using the same LinkColumn

    df = utils.sort_and_drop_duplicates(connection_dfs, shift, path_shift)
    n_uniq_links = len(df['path_index'])

    # all start positions of components
    # (max_bin + 1) is end of pangenome
    if n_uniq_links:
        dividers = np.concatenate([np.asarray([1, max_bin + 1], dtype='int32'), df["from"] + 1, df["to"]])
        dividers = np.unique(dividers).tolist()
    else:
        dividers = [1, max_bin + 1]

    LOGGER.info(f"Largest bin_id was {max_bin}; Found {len(dividers)} dividers.")

    if n_self_loops:
        LOGGER.info(f"Eliminated {n_self_loops} self-loops")

    LOGGER.info(f"Input has {n_links} listed Links.  "
          f"Segmentation eliminated {(1-n_remaining_links/n_links)*100}% of them.")
    LOGGER.info(f"Found {n_uniq_links} unique links")

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


def write_files(folder, ontology_folder, odgi_fasta: Path, schematic: PangenomeSchematic, no_adjacent_links):
    os.makedirs(folder, exist_ok=True)  # make directory for all files
    if ontology_folder:
        os.makedirs(ontology_folder, exist_ok=True)

    fasta = None
    if odgi_fasta:
        fasta = read_contigs(odgi_fasta)[0]

    bin2file_mapping = schematic.split_and_write(args.cells_per_file, folder, fasta, no_adjacent_links, ontology_folder)

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

    parser.add_argument('-p', '--parallel-cores',
                        dest='parallel_cores',
                        default=os.cpu_count(),
                        type=int,
                        help='Tip: do not set this one to more than available CPU cores)')

    parser.add_argument('-nal', '--no-adjacent-links',
                        dest='no_adjacent_links',
                        default=False,
                        action='store_true',
                        help='Switches off the add_adjacent_connector_column() routine)')

    parser.add_argument('-t', '--do-ttl',
                        dest='do_ttl',
                        default=False,
                        action='store_true',
                        help='do the ontology turtle output or not)')

    args = parser.parse_args()

    # file path logic for single or list of files with wildcard *
    if not args.output_folder:
        if args.json_file.endswith("*"):  # directory is user provided prefix
            args.output_folder = args.json_file[:-1]
        elif args.json_file.endswith(".json"):  # single json file
            args.output_folder = osPath(args.json_file).parent.joinpath(osPath(args.json_file).stem)
        else:
            print("Please provide an --out-folder or end --json-file= prefix with a *", file=sys.stderr)
            exit(1)
    else:
        args.output_folder = osPath(args.output_folder)
    os.makedirs(args.output_folder, exist_ok=True)

    if args.parallel_cores <= 0:
        args.parallel_cores = os.cpu_count()

    return args


def main():
    global args
    args = get_arguments()
    setup_logging()

    if args.parallel_cores > 0:
        chunk_size = args.parallel_cores
    else:
        chunk_size = os.cpu_count()

    parallel = Parallel(n_jobs=chunk_size, prefer="processes")

    if args.json_file.endswith("*"):
        files = glob(args.json_file + '.json')
        print("===Input Files Found===\n", '\n'.join(files))
    else:
        files = [args.json_file]

    for json_file in files:
        LOGGER.info(f'reading {osPath(json_file)}...\n')
        paths, pangenome_length, bin_width = JSONparser.parse(json_file, chunk_size*2, parallel)  # give 2x jobs to do
        schematic = segment_matrix(paths, bin_width, args.cells_per_file,
                                   pangenome_length, args.no_adjacent_links, parallel)

        # this one spits out json and optionally other output files (fasta, ttl)
        path_name = str(bin_width)
        folder_path = osPath(args.output_folder).joinpath(path_name)  # full path

        ontology_folder_path = None
        if args.do_ttl:
            ontology_folder_path = osPath(args.output_folder).joinpath(path_name + '-turtle')

        write_files(folder_path, ontology_folder_path, args.fasta, schematic, args.no_adjacent_links)

        LOGGER.info("Finished processing the file " + json_file)


def graceful_exit():
    import psutil

    # force kill the child processes
    parent = psutil.Process(os.getpid())
    for child in parent.children(recursive=True):
        child.kill()

    os._exit(exitcode)


if __name__ == '__main__':
    import atexit
    import gc

    # https://instagram-engineering.com/dismissing-python-garbage-collection-at-instagram-4dca40b29172
    # gc.disable() doesn't work, because some random 3rd-party library will
    # enable it back implicitly.
    gc.set_threshold(0)
    # Suicide immediately after other atexit functions finishes.
    # CPython will do a bunch of cleanups in Py_Finalize which
    # will again cause Copy-on-Write, including a final GC
    atexit.register(graceful_exit)

    main()

    exitcode = 0

"""
--json-file=data/run1.B1phi1.i1.seqwish.w100.json --cells-per-file=5000
--fasta=data/run1.B1phi1.i1.seqwish.fasta
For multiple files, add '*' prefix e.g. -j data/run1.B1phi1.i1.seqwish*
python segmentation.py -j data/run1.B1phi1.i1.seqwish.w1.json -f data/run1.B1phi1.i1.seqwish.fasta --cells-per-file 25000
"""
