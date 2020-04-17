""" Utility functions """
from typing import List

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
def _split_numpy_arr(connections_from, connections_to):
    groups = []
    src, dst = connections_from[0], connections_to[0]
    group_start = 0

    for i in range(len(connections_from)):
        isrc, idst = connections_from[i], connections_to[i]
        if isrc != src or idst != dst:
            groups.append((group_start, i))
            group_start = i
            src, dst = isrc, idst

    # add last group
    groups.append((group_start, len(connections_from)))
    return groups

try:
    # provides ~7x speedup on large tables
    from numba import jit, njit

    _split_numpy_arr = jit(_split_numpy_arr)
except ImportError:
    pass


def find_groups(connections_from, connections_to: np.array) -> 'List[(int, int)]':
    '''
    Returns:
      list of [start, end) indices such that for each item data[start:end] has constant value
    Args:
      connections_from(np.array): [N x 1] array of upstream data
      connections_to(np.array): [N x 1] array of downstream data
    '''
    if connections_from.size == 0:
        return []

    return _split_numpy_arr(connections_from, connections_to)


def sort_and_drop_duplicates(connections: List, shift=21, path_shift=10) -> dict:
    '''
    returns connections sorted by ["from", "to", "path_index"] without duplicate entries;
    see find_dividers in segmentation.py
    '''
    mask = (1 << shift) - 1
    lower_mask = (1 << path_shift) - 1

    # the columns are assumed to be (from, to, path_index)
    array = np.concatenate(connections, axis=-1)

    # compress all columns into a single 64-bit integer
    compressed = (np.array(array[0, :], dtype='int64', copy=False) << (shift + path_shift)) + \
                 (np.array(array[1, :], dtype='int64', copy=False) << path_shift) + \
                 array[2, :]

    compressed_no_dups = np.unique(compressed)

    return dict({
        'from': (compressed_no_dups >> (shift + path_shift)).astype('int32', copy=False),
        'to': ((compressed_no_dups >> path_shift) & mask).astype('int32', copy=False),
        'path_index': (compressed_no_dups & lower_mask).astype('int32', copy=False)
    })