"""
rdkit.Chem.EState.AtomTypes module¶
contains SMARTS definitions and calculators for EState atom types
defined in: Hall and Kier JCICS _35_ 1039-1045 (1995)  Table 1
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem

esPatterns: Incomplete

def BuildPatts(self, rawV=None):
    """
    Internal Use Only"""
    ...

def TypeAtoms(self, mol):
    """
    assigns each atom in a molecule to an EState type
    Returns:

    list of tuples (atoms can possibly match multiple patterns) with atom types"""
    ...
