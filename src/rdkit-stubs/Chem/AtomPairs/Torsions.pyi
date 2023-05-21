"""
rdkit.Chem.AtomPairs.Torsions module¶
Contains an implementation of Topological-torsion fingerprints, as
described in:
R. Nilakantan, N. Bauman, J. S. Dixon, R. Venkataraghavan;
“Topological Torsion: A New Molecular Descriptor for SAR Applications.
Comparison with Other Descriptors” JCICS 27, 82-85 (1987).
The fingerprints can be accessed through the following functions:
- GetTopologicalTorsionFingerprint
- GetHashedTopologicalTorsionFingerprint
- GetTopologicalTorsionFingerprintAsIntVect (identical to GetTopologicalTorsionFingerprint)
- GetTopologicalTorsionFingerprintAsIds
"""
from _typeshed import Incomplete
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem.AtomPairs import Utils as Utils
from rdkit.Chem.rdMolDescriptors import (
    GetHashedTopologicalTorsionFingerprint as GetHashedTopologicalTorsionFingerprint,
)
from rdkit.Chem.rdMolDescriptors import (
    GetTopologicalTorsionFingerprint as GetTopologicalTorsionFingerprint,
)

GetTopologicalTorsionFingerprintAsIntVect: Incomplete

def pyScorePath(self, mol, path, size, atomCodes=None):
    """
    Returns a score for an individual path.
    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCCC')
    >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0), 1)
    >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1), 2)
    >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2), 2)
    >>> c4 = Utils.GetAtomCode(m.GetAtomWithIdx(3), 1)
    >>> t = c1 | (c2 << rdMolDescriptors.AtomPairsParameters.codeSize) | (c3 << (rdMolDescriptors.AtomPairsParameters.codeSize * 2)) | (c4 << (rdMolDescriptors.AtomPairsParameters.codeSize * 3))
    >>> pyScorePath(m, (0, 1, 2, 3), 4) == t
    1

    The scores are path direction independent:
    >>> pyScorePath(m, (3, 2, 1, 0), 4) == t
    1

    >>> m = Chem.MolFromSmiles('C=CC(=O)O')
    >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0), 1)
    >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1), 2)
    >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2), 2)
    >>> c4 = Utils.GetAtomCode(m.GetAtomWithIdx(4), 1)
    >>> t = c1 | (c2 << rdMolDescriptors.AtomPairsParameters.codeSize) | (c3 << (rdMolDescriptors.AtomPairsParameters.codeSize * 2)) | (c4 << (rdMolDescriptors.AtomPairsParameters.codeSize * 3))
    >>> pyScorePath(m, (0, 1, 2, 4), 4) == t
    1"""
    ...

def ExplainPathScore(self, score, size=4):
    """
    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('C=CC')
    >>> score=pyScorePath(m, (0, 1, 2), 3)
    >>> ExplainPathScore(score, 3)
    (('C', 1, 0), ('C', 2, 1), ('C', 1, 1))

    Again, it’s order independent:
    >>> score = pyScorePath(m, (2, 1, 0), 3)
    >>> ExplainPathScore(score, 3)
    (('C', 1, 0), ('C', 2, 1), ('C', 1, 1))

    >>> m = Chem.MolFromSmiles('C=CO')
    >>> score=pyScorePath(m, (0, 1, 2), 3)
    >>> ExplainPathScore(score, 3)
    (('C', 1, 1), ('C', 2, 1), ('O', 1, 0))

    >>> m = Chem.MolFromSmiles('OC=CO')
    >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
    >>> ExplainPathScore(score, 4)
    (('O', 1, 0), ('C', 2, 1), ('C', 2, 1), ('O', 1, 0))

    >>> m = Chem.MolFromSmiles('CC=CO')
    >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
    >>> ExplainPathScore(score, 4)
    (('C', 1, 0), ('C', 2, 1), ('C', 2, 1), ('O', 1, 0))

    >>> m = Chem.MolFromSmiles('C=CC(=O)O')
    >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
    >>> ExplainPathScore(score, 4)
    (('C', 1, 1), ('C', 2, 1), ('C', 3, 1), ('O', 1, 1))
    >>> score=pyScorePath(m, (0, 1, 2, 4), 4)
    >>> ExplainPathScore(score, 4)
    (('C', 1, 1), ('C', 2, 1), ('C', 3, 1), ('O', 1, 0))

    >>> m = Chem.MolFromSmiles('OOOO')
    >>> score=pyScorePath(m, (0, 1, 2), 3)
    >>> ExplainPathScore(score, 3)
    (('O', 1, 0), ('O', 2, 0), ('O', 2, 0))
    >>> score=pyScorePath(m, (0, 1, 2, 3), 4)
    >>> ExplainPathScore(score, 4)
    (('O', 1, 0), ('O', 2, 0), ('O', 2, 0), ('O', 1, 0))"""
    ...

def GetTopologicalTorsionFingerprintAsIds(self, mol, targetSize=4): ...
def pyScorePath(self, mol, path, size, atomCodes=None):
    """
    Returns a score for an individual path.
    >>> from rdkit import Chem
    >>> m = Chem.MolFromSmiles('CCCCC')
    >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0), 1)
    >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1), 2)
    >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2), 2)
    >>> c4 = Utils.GetAtomCode(m.GetAtomWithIdx(3), 1)
    >>> t = c1 | (c2 << rdMolDescriptors.AtomPairsParameters.codeSize) | (c3 << (rdMolDescriptors.AtomPairsParameters.codeSize * 2)) | (c4 << (rdMolDescriptors.AtomPairsParameters.codeSize * 3))
    >>> pyScorePath(m, (0, 1, 2, 3), 4) == t
    1

    The scores are path direction independent:
    >>> pyScorePath(m, (3, 2, 1, 0), 4) == t
    1

    >>> m = Chem.MolFromSmiles('C=CC(=O)O')
    >>> c1 = Utils.GetAtomCode(m.GetAtomWithIdx(0), 1)
    >>> c2 = Utils.GetAtomCode(m.GetAtomWithIdx(1), 2)
    >>> c3 = Utils.GetAtomCode(m.GetAtomWithIdx(2), 2)
    >>> c4 = Utils.GetAtomCode(m.GetAtomWithIdx(4), 1)
    >>> t = c1 | (c2 << rdMolDescriptors.AtomPairsParameters.codeSize) | (c3 << (rdMolDescriptors.AtomPairsParameters.codeSize * 2)) | (c4 << (rdMolDescriptors.AtomPairsParameters.codeSize * 3))
    >>> pyScorePath(m, (0, 1, 2, 4), 4) == t
    1"""
    ...
