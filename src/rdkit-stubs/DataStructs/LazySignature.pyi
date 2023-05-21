"""
rdkit.DataStructs.LazySignature module¶
"""
from _typeshed import Incomplete

class LazySig(object):
    """
    computeFunc should take a single argument, the integer bit id
    to compute"""

    ...
    computeFunc: Incomplete
    size: Incomplete

    def __init__(self, computeFunc, sigSize) -> None: ...
    def __len__(self) -> int: ...
    def __getitem__(self, which): ...
