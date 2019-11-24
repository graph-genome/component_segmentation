"""This module deals with classes that describe a range from begin to end.
Span can have sections in the middle removed, creating two or less new Spans.
This is used by UniqueOnlyChainParser to track which parts of the file are untouched.
AlignedSpans use a pair of Span objects to track the coordinate frames of the
original and gapped sequence as gaps are added."""

gap_char = '-'


class Span(object):
    """ Span can have sections in the middle removed, creating two or less new Spans.
    This is used by UniqueOnlyChainParser to track which parts of the file are untouched."""
    def __init__(self, begin, end, contig_name=None, strand='+', zero_ok=True):
        assert zero_ok or begin != end, "%s %s are the same means zero length" % (begin, end)
        if not (begin >= 0 and end >= 0):
            raise ValueError("No negative indices! %i to %i" % (begin, end))
        self.begin = begin
        self.end = end
        self.contig_name = contig_name
        self.strand = strand
        assert self.strand in '+-'

    def __lt__(self, other_int):
        return self.begin < other_int

    def __contains__(self, index):
        if isinstance(index, Span):
            return self.overlaps(index)
        return self.begin <= index < self.end

    def __eq__(self, other):
        return self.begin == other.begin and self.end == other.end

    def __repr__(self):
        return ">%s:%s-%s" % (self.contig_name, '{:,}'.format(self.begin), '{:,}'.format(self.end))

    def __len__(self):
        return self.size()

    def size(self):
        return self.end - self.begin

    def overlaps(self, other):
        boundaries_check = other.begin in self or other.end - 1 in self
        is_superset = self.begin in other
        # shared_start = other.begin == self.begin
        # right_before = other.end == self.begin and other.begin != other.end
        # begin_or_end_on_wrong_side = other.end < self.begin or other.begin >= self.end
        # a = shared_start or not (begin_or_end_on_wrong_side or right_before)
        return boundaries_check or is_superset

    def split(self, split_index):
        """Splits the Span so that split_index is the first index of the second Span.
        The second span starts at split_index.  The first valid split point is begin + 1"""
        assert isinstance(self, Span), "First argument should be a Span"
        if split_index in self and split_index != self.begin + 1:
                return Span(self.begin, split_index, self.contig_name, self.strand), \
                       Span(split_index, self.end, self.contig_name, self.strand)
        raise ValueError("split_index %i is not in Span %s" % (split_index, str(self)))

    def remove_from_range(self, remove_this):
        """self is a range defined by (start, end).  Remove a middle range 'remove_this'
        with a (start, end) and you end up with a pair of two ranges on either side of the removal.
        Special casing for the removal overlapping the beginning or end."""
        assert isinstance(self, Span) and isinstance(remove_this, Span)

        # doesn't even overlap
        if not self.overlaps(remove_this):
            if not remove_this.size():
                return None, self
            if self.size():
                raise IndexError("Remove_this doesn't overlap self at all %s %s" % (str(remove_this), str(self)))
            else:  # self has no size, just throw it away
                return None, None

        first = Span(self.begin, remove_this.begin, self.contig_name, self.strand)
        second = Span(remove_this.end, self.end, self.contig_name, self.strand)

        if remove_this.begin <= self.begin and remove_this.end >= self.end:  # delete the whole thing
            return None, None
        if remove_this.begin < self.begin < remove_this.end:  # overlaps start
            return None, second
        if remove_this.end >= self.end > remove_this.begin:  # overlaps ending
            return first, None

        return first, second  # happy path

    def sample(self, sequence, error_gap_okay=True):
        try:
            return sequence[self.begin: self.end]
        except IndexError as e:
            if error_gap_okay:
                return gap_char * len(self)
            else:
                raise e

    def set_of_points(self):
        return set(range(self.begin, self.end))


