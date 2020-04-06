import pytest
import numpy as np

from matrixcomponent.matrix import Path
from matrixcomponent.utils import (
    numpy_links,
    path_boundaries,
    self_loops
)

def test_numpy_links():
    links = numpy_links([Path.LinkEntry(1, 2)])
    assert links.shape == (1, 2)
    assert links[0, 0] == 1
    assert links[0, 1] == 2


def test_path_boundaries():
    links = np.array([[0, 1], [1, 2], [2, 0]])
    assert np.array_equal(path_boundaries(links), [True, False, True])


def test_self_loops():
    links = np.array([[0, 1], [1, 1], [1, 2], [2, 1],
                      [1, 1], [1, 3], [3, 3], [3, 0]])
    loops = np.unique(links[self_loops(links)], axis=0)
    assert np.array_equal(loops, [[1, 1], [3, 3]])