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
      boolean mask of links that are considered dividers
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

    return mask


# warning: this function is intended to be numba-compatible
def _split_numpy_arr(arr):
    groups = []
    src, dst = arr[0]
    group_start = 0

    for i in range(arr.shape[0]):
        item = arr[i]
        if item[0] != src or item[1] != dst:
            groups.append((group_start, i))
            group_start = i
            src, dst = item

    # add last group
    groups.append((group_start, arr.shape[0]))
    return groups

try:
    # provides ~7x speedup on large tables
    from numba import jit
    _split_numpy_arr = jit(_split_numpy_arr)
except ImportError:
    pass


def find_groups(data: np.array) -> 'List[(int, int)]':
    '''
    Returns:
      list of [start, end) indices such that for each item data[start:end] has constant value
    Args:
      data(np.array): [N x 2] array of data; in context of segmentation each row is (upstream, downstream)
    '''
    if data.size == 0:
        return []

    return _split_numpy_arr(data)