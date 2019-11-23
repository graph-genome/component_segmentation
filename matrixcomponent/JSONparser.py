import json
import logging
import matrixcomponent.matrix as matrix

LOGGER = logging.getLogger(__name__)

# class Parser(object):
#
#     infile = ''
#
#     def __init__(self, file):
#         self.infile = file

def parse(file):
        data = []
        with open(file) as f:
            for line in f:
                data.append(json.loads(line))

        paths = []
        for path in data:
            LOGGER.info("reading " + path['path_name'])

            p = matrix.Path()

            for b in path['bins']:
                p.bins.append(p.Bin(b[0], b[1], b[2], b[3]))
            p.finalize_bins()
            for l in path['links']:
                p.links.append(p.LinkEntry(l[0], l[1]))

            paths.append(p)

        return(paths)

