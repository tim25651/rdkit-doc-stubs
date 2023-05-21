"""
rdkit.Chem.Pharm3D.ExcludedVolume module¶
"""
from _typeshed import Incomplete

class ExcludedVolume(object):
    """
    featInfo should be a sequence of ([indices],min,max) tuples"""

    ...
    index: Incomplete
    featInfo: Incomplete
    exclusionDist: Incomplete
    pos: Incomplete

    def __init__(
        self, featInfo, index: int = ..., exclusionDist: float = ...
    ) -> None: ...
