"""
Module contents¶
Module containing an assortment of functionality for basic data structures.

At the moment the data structures defined are:
Bit Vector classes (for storing signatures, fingerprints and the like:

ExplicitBitVect: class for relatively small (10s of thousands of bits) ordense bit vectors.

SparseBitVect:   class for large, sparse bit vectors

DiscreteValueVect:   class for storing vectors of integers
SparseIntVect:       class for storing sparse vectors of integers
"""
from _typeshed import Incomplete
from rdkit import rdBase as rdBase
from rdkit.DataStructs import cDataStructs as cDataStructs
from rdkit.DataStructs.cDataStructs import *

__doc__: Incomplete
similarityFunctions: Incomplete

def FingerprintSimilarity(fp1, fp2, metric=...): ...
def FoldToTargetDensity(fp, density: float = ..., minLength: int = ...): ...
