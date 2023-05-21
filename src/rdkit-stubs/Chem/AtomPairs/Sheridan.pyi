"""
rdkit.Chem.AtomPairs.Sheridan module¶
Contains an implementation of Physicochemical property fingerprints, as
described in:
Kearsley, S. K. et al.
“Chemical Similarity Using Physiochemical Property Descriptors.”
J. Chem.Inf. Model. 36, 118-127 (1996)
The fingerprints can be accessed through the following functions:
- GetBPFingerprint
- GetBTFingerprint
"""
from typing import Callable

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import GetAtomPairFingerprint as GetAtomPairFingerprint
from rdkit.Chem.rdMolDescriptors import (
    GetTopologicalTorsionFingerprint as GetTopologicalTorsionFingerprint,
)

numPathBits: Incomplete
numFpBits: Incomplete
fpLen: Incomplete

def AssignPattyTypes(self, mol, defns=None):
    """
    >>> from rdkit import Chem
    >>> AssignPattyTypes(Chem.MolFromSmiles('OCC(=O)O'))
    ['POL', 'HYD', 'OTH', 'ANI', 'ANI']"""
    ...

typMap: Incomplete

def GetBPFingerprint(self, mol, fpfn: Callable = ...):
    """
    >>> from rdkit import Chem
    >>> fp = GetBPFingerprint(Chem.MolFromSmiles('OCC(=O)O'))
    >>> fp.GetTotalVal()
    10
    >>> nze = fp.GetNonzeroElements()
    >>> sorted([(k, v) for k, v in nze.items()])
    [(32834, 1), (49219, 2), (98370, 2), (98401, 1), (114753, 2), (114786, 1), (114881, 1)]
    """
    ...

def GetBTFingerprint(self, mol, fpfn: Callable = ...):
    """
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('OCC(N)O')
    >>> AssignPattyTypes(mol)
    ['POL', 'HYD', 'HYD', 'CAT', 'POL']
    >>> fp = GetBTFingerprint(mol)
    >>> fp.GetTotalVal()
    2
    >>> nze = fp.GetNonzeroElements()
    >>> sorted([(k, v) for k, v in nze.items()])
    [(538446850..., 1), (538446852..., 1)]"""
    ...
