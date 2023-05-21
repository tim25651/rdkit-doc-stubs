"""
rdkit.Chem.EnumerateHeterocycles module¶
"""
from collections.abc import Generator

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem

def EnumerateHeterocycles(self, inputmol, depth=None):
    """
    Enumerate possible relevant heterocycles on the given input
    molecule.
    >>> from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles
    >>> from rdkit import Chem
    >>> for smi in sorted(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(Chem.MolFromSmiles('c1ccccc1'))):
    ...     print(smi)
    c1ccccc1
    c1ccncc1
    c1ccnnc1
    c1cnccn1
    c1cncnc1
    c1cnncn1
    c1cnnnc1
    c1ncncn1

    The algorithm works by mutating only one atom at a time. The depth
    parameter can be used to control the level of this recursion. For
    example, only enumerating aromatic rings that are one atom different:
    >>> for smi in sorted(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(Chem.MolFromSmiles('n1ccccc1'), depth=1)):
    ...     print(smi)
    c1ccccc1
    c1ccnnc1
    c1cnccn1
    c1cncnc1"""
    ...

def GetHeterocycleReactionSmarts(self):
    """
    Return a list of the individual patterns for mutating individual
    atoms in aromatic rings to generate a new aromatic ring. The
    elements are collections.namedtuple objects with the following
    fields:
    SMARTS - the left side of the reaction SMARTS pattern, matches the atom to mutate
    CONVERT_FROM - the element type being converted: c, n, o, s
    CONVERT_TO - the right side of the reaction SMARTS pattern, there can be multiple destination types separated with a comma, these should map to multiple actual reaction objects
    EXAMPLE - an example aromatic ring system that SMARTS should match against, used in test cases
    NEGATIVE_EXAMPLE - an example aromatic ring system that SMART should NOT match against, used in test cases
    DESCRIPTION - a human readable description of the SMARTS pattern matching"""
    ...

REACTION_CACHE: Incomplete

def GetHeterocycleReactions(self):
    """
    Return RDKit ChemicalReaction objects of the reaction SMARTS
    returned from GetHeterocyleReactionSmarts."""
    ...

def EnumerateHeterocycles(self, inputmol, depth=None):
    """
    Enumerate possible relevant heterocycles on the given input
    molecule.
    >>> from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles
    >>> from rdkit import Chem
    >>> for smi in sorted(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(Chem.MolFromSmiles('c1ccccc1'))):
    ...     print(smi)
    c1ccccc1
    c1ccncc1
    c1ccnnc1
    c1cnccn1
    c1cncnc1
    c1cnncn1
    c1cnnnc1
    c1ncncn1

    The algorithm works by mutating only one atom at a time. The depth
    parameter can be used to control the level of this recursion. For
    example, only enumerating aromatic rings that are one atom different:
    >>> for smi in sorted(Chem.MolToSmiles(m) for m in EnumerateHeterocycles(Chem.MolFromSmiles('n1ccccc1'), depth=1)):
    ...     print(smi)
    c1ccccc1
    c1ccnnc1
    c1cnccn1
    c1cncnc1"""
    ...
