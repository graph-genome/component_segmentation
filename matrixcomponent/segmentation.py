from typing import List
from matrixcomponent.matrix import Path


#Tasks
# Get example output for tests - Ted
# Parser for ODGI Bin file format - Ted
# Component Segmentation Detection - Josiah and Joerg
#   Python memory object model
# Output format


def segment_matrix(matrix: List[Path]):
    dividers = set()  # list of indices of new components

    for path in matrix:  # TODO: brilliant except for the fact that bin_id is already sorted...
        for x, bin in enumerate(path.bins):  # For each node
            previous = path.bins[x - 1].bin_id
            if x and bin.bin_id > previous + 1:  # Is there a gap?
                missing_range = list(range(previous+1, bin.bin_id))
                # Is the gap range anywhere else in this individual?
                if any([i in path for i in missing_range if i > 0]):
                    dividers.add(bin.bin_id + 1)  # the first position of the new component
                    # TODO: insert prevarications about exact position
                    # Divider should be somewhere in here
                    # Tolerable range?
                    # Stack up others using the same Link
