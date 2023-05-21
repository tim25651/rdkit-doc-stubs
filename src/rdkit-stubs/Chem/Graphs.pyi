"""
rdkit.Chem.Graphs module¶
Python functions for manipulating molecular graphs
In theory much of the functionality in here should be migrating into the
C/C++ codebase.
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs

def CharacteristicPolynomial(self, mol, mat=None):
    """
    calculates the characteristic polynomial for a molecular graph
    if mat is not passed in, the molecule’s Weighted Adjacency Matrix will
    be used.
    The approach used is the Le Verrier-Faddeev-Frame method described
    in _Chemical Graph Theory, 2nd Edition_ by Nenad Trinajstic (CRC Press,
    1992), pg 76."""
    ...
