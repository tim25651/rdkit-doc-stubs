"""
rdkit.Chem.Lipinski module¶
Calculation of Lipinski parameters for molecules
"""
from typing import Callable, overload

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem.rdchem import Mol

@overload
def FractionCSP3(self, x, y: Callable = ...):
    """
    returns the fraction of C atoms that are SP3 hybridized

    C++ signature :double CalcFractionCSP3(RDKit::ROMol)"""
    ...

@overload
def FractionCSP3(self, mol: Mol) -> float: ...
def HeavyAtomCount(self, mol):
    """
    Number of heavy atoms a molecule."""
    ...

def NHOHCount(self, x):
    """
    Number of NHs or OHs"""
    ...

def NOCount(self, x):
    """
    Number of Nitrogens and Oxygens"""
    ...

@overload
def NumAliphaticCarbocycles(self, x, y: Callable = ...):
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) carbocycles for a molecule

    C++ signature :unsigned int CalcNumAliphaticCarbocycles(RDKit::ROMol)"""
    ...

@overload
def NumAliphaticCarbocycles(self, mol: Mol) -> int: ...
@overload
def NumAliphaticHeterocycles(self, x, y: Callable = ...):
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) heterocycles for a molecule

    C++ signature :unsigned int CalcNumAliphaticHeterocycles(RDKit::ROMol)"""
    ...

@overload
def NumAliphaticHeterocycles(self, mol: Mol) -> int: ...
@overload
def NumAliphaticRings(self, x, y: Callable = ...):
    """
    returns the number of aliphatic (containing at least one non-aromatic bond) rings for a molecule

    C++ signature :unsigned int CalcNumAliphaticRings(RDKit::ROMol)"""
    ...

@overload
def NumAliphaticRings(self, mol: Mol) -> int: ...
@overload
def NumAromaticCarbocycles(self, x, y: Callable = ...):
    """
    returns the number of aromatic carbocycles for a molecule

    C++ signature :unsigned int CalcNumAromaticCarbocycles(RDKit::ROMol)"""
    ...

@overload
def NumAromaticCarbocycles(self, mol: Mol) -> int: ...
@overload
def NumAromaticHeterocycles(self, x, y: Callable = ...):
    """
    returns the number of aromatic heterocycles for a molecule

    C++ signature :unsigned int CalcNumAromaticHeterocycles(RDKit::ROMol)"""
    ...

@overload
def NumAromaticHeterocycles(self, mol: Mol) -> int: ...
@overload
def NumAromaticRings(self, x, y: Callable = ...):
    """
    returns the number of aromatic rings for a molecule

    C++ signature :unsigned int CalcNumAromaticRings(RDKit::ROMol)"""
    ...

@overload
def NumAromaticRings(self, mol: Mol) -> int: ...
def NumHAcceptors(self, x):
    """
    Number of Hydrogen Bond Acceptors"""
    ...

HDonorSmarts: Incomplete
HAcceptorSmarts: Incomplete
HeteroatomSmarts: Incomplete
RotatableBondSmarts: Incomplete
NHOHSmarts: Incomplete
NOCountSmarts: Incomplete

def NumHDonors(self, x):
    """
    Number of Hydrogen Bond Donors"""
    ...

def NumHAcceptors(self, x):
    """
    Number of Hydrogen Bond Acceptors"""
    ...

def NumHeteroatoms(self, x):
    """
    Number of Heteroatoms"""
    ...

def NumRotatableBonds(self, x):
    """
    Number of Rotatable Bonds"""
    ...

@overload
def NumSaturatedCarbocycles(self, x, y: Callable = ...):
    """
    returns the number of saturated carbocycles for a molecule

    C++ signature :unsigned int CalcNumSaturatedCarbocycles(RDKit::ROMol)"""
    ...

@overload
def NumSaturatedCarbocycles(self, mol: Mol) -> int: ...
@overload
def NumSaturatedHeterocycles(self, x, y: Callable = ...):
    """
    returns the number of saturated heterocycles for a molecule

    C++ signature :unsigned int CalcNumSaturatedHeterocycles(RDKit::ROMol)"""
    ...

@overload
def NumSaturatedHeterocycles(self, mol: Mol) -> int: ...
@overload
def NumSaturatedRings(self, x, y: Callable = ...):
    """
    returns the number of saturated rings for a molecule

    C++ signature :unsigned int CalcNumSaturatedRings(RDKit::ROMol)"""
    ...

@overload
def NumSaturatedRings(self, mol: Mol) -> int: ...
def NOCount(self, x):
    """
    Number of Nitrogens and Oxygens"""
    ...

def NHOHCount(self, x):
    """
    Number of NHs or OHs"""
    ...

def RingCount(self, x): ...
def HeavyAtomCount(self, mol):
    """
    Number of heavy atoms a molecule."""
    ...

nm: Incomplete
