""" Utility functions """
from typing import List

import numpy as np
import pandas as pd


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
    from numba import jit, njit

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


def sort_and_drop_duplicates(connections: 'pd.DataFrame', shift=21, path_shift=10) -> 'pd.DataFrame':
    '''
    returns connections sorted by ["from", "to", "path_index"] without duplicate entries;
    see find_dividers in segmentation.py
    '''
    mask = (1 << shift) - 1
    lower_mask = (1 << path_shift) - 1

    if np.any(connections.max() > mask):
        # nigh impossible with the default limit: (1 << 21) = 2M bins / paths;
        # as such, this line of code is mostly for illustration purposes
        return connections.drop_duplicates().sort_values(by=["from", "to", "path_index"])

    # the columns are assumed to be (from, to, path_index)
    array = connections.to_numpy()

    # compress all columns into a single 64-bit integer
    compressed = (array[:, 0] << (shift + path_shift)) + (array[:, 1] << path_shift) + array[:, 2]
    compressed_no_dups = np.unique(compressed)

    return pd.DataFrame.from_dict({
        'from': (compressed_no_dups >> (shift + path_shift)).astype('int32', copy=False),
        'to': ((compressed_no_dups >> path_shift) & mask).astype('int32', copy=False),
        'path_index': (compressed_no_dups & lower_mask).astype('int32', copy=False)
    })