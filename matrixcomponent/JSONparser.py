import json
import logging
import matrixcomponent.matrix as matrix

LOGGER = logging.getLogger(__name__)

def parse(file):
        data = []
        with open(file) as f:
            for line in f:
                data.append(json.loads(line))

        paths = []
        for path in data:
            if "path_name" in path:
                LOGGER.info("reading " + path['path_name'])

                p = matrix.Path(path['path_name'])

                for b in path['bins']:
                    p.bins.append(p.Bin(b[0], b[1], b[2], b[3]))
                p.finalize_bins()
                
                for l in path['links']:
                    p.links.append(p.LinkEntry(l[0], l[1]))

                paths.append(p)

        return(paths)

