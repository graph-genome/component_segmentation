from matrixcomponent.matrix import Path


def test_construction():
    p = Path()
    p.bins = [Path.Bin(1, .9, 0, 0.11), Path.Bin(2, .8, 0, 0.15),]
    p.links = [Path.LinkEntry(1, 2)]

