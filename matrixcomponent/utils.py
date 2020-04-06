""" Utility functions """

import numpy as np


def path_boundaries(links: np.array) -> np.array:
    ''' boolean mask of links which have 0 as either upstream or downstream '''
    # Links to 0 Bin indicate the beginning or end of a path.  0 Bin has no sequence
    return np.any(links == 0, axis=1)


def self_loops(links: np.array) -> np.array:
    ''' boolean mask of links where upstream equals downstream '''
    return links[:, 0] == links[:, 1]


def path_dividers(links: np.array, bin_ids: np.array) -> np.array:
    '''
    Returns:
      subarray of links that are considered dividers
    Args:
      links(np.array) - array of path links
      bin_ids(np.array) - sorted array of path bin ids
    '''
    mask = links[:, 0] > links[:, 1]  # simple case: upstream > downstream

    # if downstream >= upstream, we consider a link to be a divider
    # iff there exists a bin id within [upstream + 1, downstream) range
    ranges = np.copy(links[~mask])
    ranges[:, 0] += 1   # each row is [upstream + 1, downstream]
    indices = np.searchsorted(bin_ids, ranges)

    # maintain the original order just in case
    mask[~mask] = indices[:, 1] > indices[:, 0]

    return links[mask]
