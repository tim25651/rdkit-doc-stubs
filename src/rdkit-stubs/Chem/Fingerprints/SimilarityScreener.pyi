"""
rdkit.Chem.Fingerprints.SimilarityScreener module¶
class definitions for similarity screening
See _SimilarityScreener_ for overview of required API
"""
from _typeshed import Incomplete
from rdkit import DataStructs as DataStructs
from rdkit.DataStructs import TopNContainer as TopNContainer

class SimilarityScreener(object):
    """
    base class

    important attributes:
    probe: the probe fingerprint against which we screen.

    metric: a function that takes two arguments and returns a similaritymeasure between them

    dataSource: the source pool from which to draw, needs to supporta next() method

    fingerprinter: a function that takes a molecule and returns afingerprint of the appropriate format

    Notessubclasses must support either an iterator interface
    or __len__ and __getitem__"""

    def GetSingleFingerprint(self, probe):
        """
        returns a fingerprint for a single probe object
        This is potentially useful in initializing our internal
        probe object."""
        ...
    metric: Incomplete
    dataSource: Incomplete
    fingerprinter: Incomplete
    probe: Incomplete

    def __init__(
        self,
        probe: Incomplete | None = ...,
        metric: Incomplete | None = ...,
        dataSource: Incomplete | None = ...,
        fingerprinter: Incomplete | None = ...,
    ) -> None: ...
    def Reset(self):
        """
        used to reset screeners that behave as iterators"""
        ...
    def SetProbe(self, probeFingerprint):
        """
        sets our probe fingerprint"""
        ...
    def GetSingleFingerprint(self, probe):
        """
        returns a fingerprint for a single probe object
        This is potentially useful in initializing our internal
        probe object."""
        ...

class ThresholdScreener(SimilarityScreener):
    """
    Used to return all compounds that have a similarityto the probe beyond a threshold value

    Notes:

    This is as lazy as possible, so the data source isn’t
    queried until the client asks for a hit.
    In addition to being lazy, this class is as thin as possible.
    (Who’d have thought it was possible!)
    Hits are not stored locally, so if a client resets
    the iteration and starts over, the same amount of work must
    be done to retrieve the hits.
    The thinness and laziness forces us to support only forward
    iteration (not random access)"""

    threshold: Incomplete
    dataIter: Incomplete

    def __init__(self, threshold, **kwargs) -> None: ...
    def Reset(self):
        """
        used to reset our internal state so that iteration
        starts again from the beginning"""
        ...
    def __iter__(self): ...
    def next(self):
        """
        required part of iterator interface"""
        ...
    __next__ = next

class TopNScreener(SimilarityScreener):
    """
    A screener that only returns the top N hits found
    Notes

    supports forward iteration and getitem"""

    numToGet: Incomplete
    topN: Incomplete

    def __init__(self, num, **kwargs) -> None: ...
    def Reset(self):
        """
        used to reset screeners that behave as iterators"""
        ...
    def __iter__(self): ...
    def next(self): ...
    __next__ = next

    def __len__(self) -> int: ...
    def __getitem__(self, idx): ...
