import numpy as np

from matrixcomponent.utils import (
    path_boundaries,
    self_loops,
    path_dividers
)


def test_path_boundaries():
    links = np.array([[0, 1], [1, 2], [2, 0]])
    assert np.array_equal(path_boundaries(links), [True, False, True])


def test_self_loops():
    links = np.array([[0, 1], [1, 1], [1, 2], [2, 1],
                      [1, 1], [1, 3], [3, 3], [3, 0]])
    loops = np.unique(links[self_loops(links)], axis=0)
    assert np.array_equal(loops, [[1, 1], [3, 3]])


def test_path_dividers():
    bins = [1, 3, 4, 8, 15]

    # test the simple case separately
    links = np.array([[3, 1], [4, 1]])
    assert np.array_equal(links[path_dividers(links, bins)], links)

    # test the general case
    links = np.array([
        [3, 1],  # u > d => yes
        [1, 3],  # no bins within [1 + 1, 3) range => no
        [3, 5],  # bin 4 is within [3 + 1, 5) range => yes
        [5, 9],  # bin 8 is within [5 + 1, 9) => yes
        [5, 8],  # no bins within [5 + 1, 9) => yes
        [3, 4],  # empty bin range [3 + 1, 4) => no
        [4, 1],  # u > d => yes
    ])
    assert np.array_equal(links[path_dividers(links, bins)],
                          [[3, 1], [3, 5], [5, 9], [4, 1]])
