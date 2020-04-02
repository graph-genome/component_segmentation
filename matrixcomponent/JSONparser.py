import json
import logging
import matrixcomponent.matrix as matrix
from matrixcomponent import ODGI_VERSION

LOGGER = logging.getLogger(__name__)


def parse(file):
        data = []
        with open(file) as f:
            for line in f:
                data.append(json.loads(line))

        paths = []
        pangenome_length = 0
        for path in data:
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
                
                for l in path['links']:
                    p.links.append(p.LinkEntry(l[0], l[1]))

                paths.append(p)

        return(paths, pangenome_length, bin_width)

