import json
import logging

from joblib import delayed

import matrixcomponent.matrix as matrix
from matrixcomponent import ODGI_VERSION, utils

import numpy as np

LOGGER = logging.getLogger(__name__)


def process_path(line=None):
    pangenome_length = -1
    bin_width = -1
    p = None

    path = json.loads(line)

    if "odgi_version" in path:
        # this is the header
        assert path["odgi_version"] == ODGI_VERSION, f"Expecting version {ODGI_VERSION}." \
                                                     f" This version added the nucleotide ranges for each bin for each path."
        print(f"Found file with {path['pangenome_length']} nucleotides in the pangenome and"
             f" a {path['bin_width']}bp bin width.", flush=True)
        pangenome_length = path['pangenome_length']
        bin_width = path['bin_width']

    if "path_name" in path:
        LOGGER.info("reading " + path['path_name'])

        p = matrix.Path(path['path_name'])

        for b in path['bins']:
            ranges = b[4]
            if type(ranges) is not list and len(b) >= 6:
                ranges = [[b[4], b[5]]]

            compressed_ranges = []
            for r in ranges:
                compressed_ranges.extend(r)

            bin = matrix.Bin(b[0], b[1], b[2], compressed_ranges, 0, b[3])
            p.bins.setdefault(bin.bin_id, bin)

        # do the major part of the segmentation.find_dividers() method
        links = np.asarray(path['links'], dtype='int32')
        p.num_links = len(links)

        p.path_dividers = np.array([], dtype='int32')
        p.self_loops = np.array([], dtype='int32')

        p.max_bin_id = 1
        bin_ids = np.asarray(p.bins.keys()) # already sorted
        if bin_ids.size > 0:
            p.max_bin_id = int(bin_ids[-1])

        if links.size > 0:
            # we don't want these to become dividers
            boundary_mask = utils.path_boundaries(links)
            self_loops_mask = utils.self_loops(links)

            if np.any(self_loops_mask):
                p.self_loops = links[self_loops_mask]

            links = links[~(boundary_mask | self_loops_mask)]

            path_dividers_mask = utils.path_dividers(links, bin_ids)
            p.path_dividers = links[path_dividers_mask]

    return [pangenome_length, bin_width, p]


def do_processing(parallel, lines, pangenome_length, bin_width, paths):
    if parallel is None:
        results = [process_path(line) for line in lines]  # serial version
    else:
        results = parallel(delayed(process_path)(line) for line in lines)

    for res in results:
        plen, bwidth, p = res[0], res[1], res[2]
        if plen > -1:
            pangenome_length[0] = plen
        if bwidth > -1:
            bin_width[0] = bwidth
        if p is not None:
            paths.append(p)


def parse(file, chunk_size, parallel):
    paths = []
    pangenome_length = [0]
    bin_width = [0]

    lines = []
    with open(file) as f:
        for line in f:
            # early filtering saves time
            if line.startswith("odgi", 2, 10) or line.startswith("path", 2, 10):
                lines.append(line)

            # wait until there's enough things to do
            if len(lines) < chunk_size:
                continue

            do_processing(parallel, lines, pangenome_length, bin_width, paths)
            lines.clear() # start collecting the next block

# and process the leftovers
    if len(lines) > 0:
        do_processing(parallel, lines, pangenome_length, bin_width, paths)

    return (paths, pangenome_length[0], bin_width[0])
