"""
rdkit.DataStructs.TopNContainer module¶
"""
from _typeshed import Incomplete

class TopNContainer(object):
    """
    maintains a sorted list of a particular number of data elements.
    if size is negative, all entries will be kept in sorted order"""

    def GetExtras(self):
        """
        returns our set of extras"""
        ...
    def GetPts(self):
        """
        returns our set of points"""
        ...
    best: Incomplete
    extras: Incomplete

    def __init__(self, size, mostNeg=...) -> None: ...
    def Insert(self, val, extra=None):
        """
        only does the insertion if val fits"""
        ...
    def GetPts(self):
        """
        returns our set of points"""
        ...
    def GetExtras(self):
        """
        returns our set of extras"""
        ...
    def __len__(self) -> int: ...
    def __getitem__(self, which): ...
    def reverse(self): ...
