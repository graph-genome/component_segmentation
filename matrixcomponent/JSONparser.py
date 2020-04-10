import json
import logging

import os
from joblib import Parallel, delayed

import matrixcomponent.matrix as matrix
from matrixcomponent import ODGI_VERSION
from itertools import islice

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

            p.bins.append(p.Bin(b[0], b[1], b[2], ranges))
        p.finalize_bins()

        p.links = np.array(path['links'])

    return [pangenome_length, bin_width, p]


def parse(file, parallel_cores):
    paths = []
    pangenome_length = 0
    bin_width = 0

    if parallel_cores > 0:
        chunk_size = parallel_cores
    else:
        chunk_size = os.cpu_count()

    parallel = Parallel(n_jobs=chunk_size, prefer="processes")
    with open(file) as f:
        while True:
            lines = list(islice(f, chunk_size))
            if len(lines) == 0:
                break

            # use this one for serial version
            # results = [process_path(line) for line in lines]
            results = parallel(delayed(process_path)(line) for line in lines)
            for res in results:
                plen, bwidth, p = res[0], res[1], res[2]
                if plen > -1:
                    pangenome_length = plen
                if bwidth > -1:
                    bin_width = bwidth
                if p is not None:
                    paths.append(p)

    return (paths, pangenome_length, bin_width)
