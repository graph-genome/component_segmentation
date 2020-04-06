""" Utility functions """

import numpy as np

def numpy_links(links: 'List[Path.LinkEntry]') -> np.array:
    np_links = np.zeros((len(links), 2), dtype=int)
    for i, link in enumerate(links):
        np_links[i, 0] = link.upstream
        np_links[i, 1] = link.downstream
    return np_links


def path_boundaries(links: np.array) -> np.array:
    ''' boolean mask of links which have 0 as either upstream or downstream '''
    # Links to 0 Bin indicate the beginning or end of a path.  0 Bin has no sequence
    return np.any(links == 0, axis=1)


def self_loops(links: np.array) -> np.array:
    ''' boolean mask of links where upstream equals downstream '''
    return links[:, 0] == links[:, 1]
