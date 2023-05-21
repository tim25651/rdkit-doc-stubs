"""
rdkit.Chem.Scaffolds.MurckoScaffold module¶
Generation of Murcko scaffolds from a molecule
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem

def GetScaffoldForMol(self, mol):
    """
    Return molecule object containing scaffold of mol
    >>> m = Chem.MolFromSmiles('Cc1ccccc1')
    >>> GetScaffoldForMol(m)
    <rdkit.Chem.rdchem.Mol object at 0x...>
    >>> Chem.MolToSmiles(GetScaffoldForMol(m))
    'c1ccccc1'

    >>> m = Chem.MolFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
    >>> Chem.MolToSmiles(GetScaffoldForMol(m))
    'c1ccc(Oc2ccccn2)cc1'"""
    ...

murckoTransforms: Incomplete

def MakeScaffoldGeneric(self, mol):
    """
    Makes a Murcko scaffold generic (i.e. all atom types->C and all bonds ->single
    >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('c1ccccc1')))
    'C1CCCCC1'
    >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles('c1ncccc1')))
    'C1CCCCC1'

    The following were associated with sf.net issue 246
    >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles(‘c1[nH]ccc1’)))
    ‘C1CCCC1’
    >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles(‘C1[NH2+]C1’)))
    ‘C1CC1’
    >>> Chem.MolToSmiles(MakeScaffoldGeneric(Chem.MolFromSmiles(‘C1[C@](Cl)(F)O1’)))
    ‘CC1(C)CC1’"""
    ...

murckoPatts: Incomplete
murckoQ: Incomplete
aromaticNTransform: Incomplete

def GetScaffoldForMol(self, mol):
    """
    Return molecule object containing scaffold of mol
    >>> m = Chem.MolFromSmiles('Cc1ccccc1')
    >>> GetScaffoldForMol(m)
    <rdkit.Chem.rdchem.Mol object at 0x...>
    >>> Chem.MolToSmiles(GetScaffoldForMol(m))
    'c1ccccc1'

    >>> m = Chem.MolFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
    >>> Chem.MolToSmiles(GetScaffoldForMol(m))
    'c1ccc(Oc2ccccn2)cc1'"""
    ...

def MurckoScaffoldSmiles(self, smiles=None, mol=None, includeChirality=False):
    """
    Returns MurckScaffold Smiles from smiles
    >>> MurckoScaffoldSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
    'c1ccc(Oc2ccccn2)cc1'

    >>> MurckoScaffoldSmiles(mol=Chem.MolFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1'))
    'c1ccc(Oc2ccccn2)cc1'"""
    ...

def MurckoScaffoldSmilesFromSmiles(self, smiles, includeChirality=False):
    """
    Returns MurckScaffold Smiles from smiles
    >>> MurckoScaffoldSmilesFromSmiles('Cc1cc(Oc2nccc(CCC)c2)ccc1')
    'c1ccc(Oc2ccccn2)cc1'"""
    ...
