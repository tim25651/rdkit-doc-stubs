"""
rdkit.Chem.EState.Fingerprinter module¶
EState fingerprinting
"""
from rdkit.Chem.EState import AtomTypes as AtomTypes
from rdkit.Chem.EState import EStateIndices as EStateIndices

def FingerprintMol(self, mol):
    """
    generates the EState fingerprints for the molecule
    Concept from the paper: Hall and Kier JCICS _35_ 1039-1045 (1995)

    two numeric arrays are returned:The first (of ints) contains the number of times each possible atom type is hit
    The second (of floats) contains the sum of the EState indices for atoms of

    each type."""
    ...
