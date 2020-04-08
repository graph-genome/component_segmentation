import json
import logging

import ray
from ray.thirdparty_files import psutil

import matrixcomponent.matrix as matrix
from matrixcomponent import ODGI_VERSION
from itertools import islice

import numpy as np

LOGGER = logging.getLogger(__name__)


@ray.remote(num_return_vals=1)
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
            p.bins.append(p.Bin(b[0], b[1], b[2], b[4], b[5]))
        p.finalize_bins()

        p.links = np.array(path['links'])

    return [pangenome_length, bin_width, p]


def parse(file):
    paths = []
    pangenome_length = 0
    bin_width = 0

    chunk_size = psutil.cpu_count()

    ray.init()

    with open(file) as f:
        while True:
            lines = list(islice(f, chunk_size))
            if len(lines) == 0:
                break

            results = []
            for line in lines:
                res = process_path.remote(line)
                results.append(res)

            results = ray.get(results)
            for res in results:
                plen, bwidth, p = res[0], res[1], res[2]
                if plen > -1:
                    pangenome_length = plen
                if bwidth > -1:
                    bin_width = bwidth
                if p is not None:
                    paths.append(p)

    return (paths, pangenome_length, bin_width)
