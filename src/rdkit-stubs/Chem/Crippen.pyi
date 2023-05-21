"""
rdkit.Chem.Crippen module¶
Atom-based calculation of LogP and MR using Crippen’s approach

Reference:

Wildman and G. M. Crippen JCICS _39_ 868-873 (1999)
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors

defaultPatternFileName: Incomplete

def MolLogP(self, *x, **y):
    """
    Wildman-Crippen LogP value

    Uses an atom-based scheme based on the values in the paper:

    Wildman and G. M. Crippen JCICS 39 868-873 (1999)

    Arguments

    inMol: a molecule
    addHs: (optional) toggles adding of Hs to the molecule for the calculation.
    If true, hydrogens will be added to the molecule and used in the calculation."""
    ...

def MolMR(self, *x, **y):
    """
    Wildman-Crippen MR value

    Uses an atom-based scheme based on the values in the paper:

    Wildman and G. M. Crippen JCICS 39 868-873 (1999)

    Arguments

    inMol: a molecule
    addHs: (optional) toggles adding of Hs to the molecule for the calculation.
    If true, hydrogens will be added to the molecule and used in the calculation."""
    ...
