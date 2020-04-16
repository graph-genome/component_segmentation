import json
import logging

from joblib import delayed

import matrixcomponent.matrix as matrix
from matrixcomponent import ODGI_VERSION

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
                                                     f"This version added the header with pangenome nucleotide count."
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

            bin = matrix.Bin(b[0], b[1], b[2], ranges)
            p.bins.setdefault(bin.bin_id, bin)

        p.links = np.asarray(path['links'], dtype='int32')

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

    return (paths, pangenome_length[0], bin_width[0], parallel)
