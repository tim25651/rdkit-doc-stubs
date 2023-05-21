"""
rdkit.DataStructs.BitEnsemble module¶
#DOC
"""
from _typeshed import Incomplete

class BitEnsemble(object):
    """
    used to store a collection of bits and score
    BitVects (or signatures) against them."""

    def __init__(self, bits: Incomplete | None = ...) -> None: ...
    def SetBits(self, bits): ...
    def AddBit(self, bit): ...
    def GetBits(self): ...
    def GetNumBits(self): ...
    def ScoreWithIndex(self, other):
        """
        other must support __getitem__()"""
        ...
    def ScoreWithOnBits(self, other):
        """
        other must support GetOnBits()"""
        ...
    def SetBits(self, bits): ...
    def ScoreWithIndex(self, other):
        """
        other must support __getitem__()"""
        ...