"""
rdkit.Chem.TorsionFingerprints module¶
Torsion Fingerprints (Deviation) (TFD)
According to a paper from Schulz-Gasch et al., JCIM, 52, 1499-1512 (2012).
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import Geometry as Geometry
from rdkit import RDConfig as RDConfig
from rdkit import rdBase as rdBase
from rdkit.Chem import rdchem as rdchem
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors

def CalculateTFD(self, torsions1, torsions2, weights=None):
    """
    Calculate the torsion deviation fingerprint (TFD) given two lists of
    torsion angles.
    Arguments:
    - torsions1:  torsion angles of conformation 1
    - torsions2:  torsion angles of conformation 2
    - weights:    list of torsion weights (default: None)
    Return: TFD value (float)"""
    ...

def CalculateTorsionAngles(self, mol, tors_list, tors_list_rings, confId=-1):
    """
    Calculate the torsion angles for a list of non-ring and
    a list of ring torsions.
    Arguments:
    - mol:       the molecule of interest
    - tors_list: list of non-ring torsions
    - tors_list_rings: list of ring torsions
    - confId:    index of the conformation (default: first conformer)
    Return: list of torsion angles"""
    ...

def CalculateTorsionLists(
    self, mol, maxDev="equal", symmRadius=2, ignoreColinearBonds=True
):
    """
    Calculate a list of torsions for a given molecule. For each torsion
    the four atom indices are determined and stored in a set.
    Arguments:
    - mol:      the molecule of interest
    - maxDev:   maximal deviation used for normalization

    ‘equal’: all torsions are normalized using 180.0 (default)
    ‘spec’:  each torsion is normalized using its specific

    maximal deviation as given in the paper

    symmRadius: radius used for calculating the atom invariants(default: 2)

    ignoreColinearBonds: if True (default), single bonds adjacent totriple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: two lists of torsions: non-ring and ring torsions"""
    ...

def CalculateTorsionAngles(self, mol, tors_list, tors_list_rings, confId=-1):
    """
    Calculate the torsion angles for a list of non-ring and
    a list of ring torsions.
    Arguments:
    - mol:       the molecule of interest
    - tors_list: list of non-ring torsions
    - tors_list_rings: list of ring torsions
    - confId:    index of the conformation (default: first conformer)
    Return: list of torsion angles"""
    ...

def CalculateTorsionWeights(self, mol, aid1=-1, aid2=-1, ignoreColinearBonds=True):
    """
    Calculate the weights for the torsions in a molecule.
    By default, the highest weight is given to the bond
    connecting the two most central atoms.
    If desired, two alternate atoms can be specified (must
    be connected by a bond).
    Arguments:
    - mol:   the molecule of interest
    - aid1:  index of the first atom (default: most central)
    - aid2:  index of the second atom (default: second most central)
    - ignoreColinearBonds: if True (default), single bonds adjacent to

    triple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: list of torsion weights (both non-ring and ring)"""
    ...

def GetBestTFDBetweenMolecules(
    self,
    mol1,
    mol2,
    confId1=-1,
    useWeights=True,
    maxDev="equal",
    symmRadius=2,
    ignoreColinearBonds=True,
):
    """
    Wrapper to calculate the best TFD between a single conformer of mol1 and all the conformers of mol2
    Important: The two molecules must be isomorphic
    Arguments:
    - mol1:     first instance of the molecule of interest
    - mol2:     second instance the molecule of interest
    - confId1:  conformer index for mol1 (default: first conformer)
    - useWeights: flag for using torsion weights in the TFD calculation
    - maxDev:   maximal deviation used for normalization

    ‘equal’: all torsions are normalized using 180.0 (default)
    ‘spec’:  each torsion is normalized using its specific

    maximal deviation as given in the paper

    symmRadius: radius used for calculating the atom invariants(default: 2)

    ignoreColinearBonds: if True (default), single bonds adjacent totriple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: TFD value"""
    ...

def CalculateTFD(self, torsions1, torsions2, weights=None):
    """
    Calculate the torsion deviation fingerprint (TFD) given two lists of
    torsion angles.
    Arguments:
    - torsions1:  torsion angles of conformation 1
    - torsions2:  torsion angles of conformation 2
    - weights:    list of torsion weights (default: None)
    Return: TFD value (float)"""
    ...

def GetTFDBetweenConformers(
    self,
    mol,
    confIds1,
    confIds2,
    useWeights=True,
    maxDev="equal",
    symmRadius=2,
    ignoreColinearBonds=True,
):
    """
    Wrapper to calculate the TFD between two list of conformers
    of a molecule
    Arguments:
    - mol:      the molecule of interest
    - confIds1:  first list of conformer indices
    - confIds2:  second list of conformer indices
    - useWeights: flag for using torsion weights in the TFD calculation
    - maxDev:   maximal deviation used for normalization

    ‘equal’: all torsions are normalized using 180.0 (default)
    ‘spec’:  each torsion is normalized using its specific

    maximal deviation as given in the paper

    symmRadius: radius used for calculating the atom invariants(default: 2)

    ignoreColinearBonds: if True (default), single bonds adjacent totriple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: list of TFD values"""
    ...

def GetTFDBetweenMolecules(
    self,
    mol1,
    mol2,
    confId1=-1,
    confId2=-1,
    useWeights=True,
    maxDev="equal",
    symmRadius=2,
    ignoreColinearBonds=True,
):
    """
    Wrapper to calculate the TFD between two molecules.
    Important: The two molecules must be isomorphic
    Arguments:
    - mol1:     first instance of the molecule of interest
    - mol2:     second instance the molecule of interest
    - confId1:  conformer index for mol1 (default: first conformer)
    - confId2:  conformer index for mol2 (default: first conformer)
    - useWeights: flag for using torsion weights in the TFD calculation
    - maxDev:   maximal deviation used for normalization

    ‘equal’: all torsions are normalized using 180.0 (default)
    ‘spec’:  each torsion is normalized using its specific

    maximal deviation as given in the paper

    symmRadius: radius used for calculating the atom invariants(default: 2)

    ignoreColinearBonds: if True (default), single bonds adjacent totriple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: TFD value"""
    ...

def GetBestTFDBetweenMolecules(
    self,
    mol1,
    mol2,
    confId1=-1,
    useWeights=True,
    maxDev="equal",
    symmRadius=2,
    ignoreColinearBonds=True,
):
    """
    Wrapper to calculate the best TFD between a single conformer of mol1 and all the conformers of mol2
    Important: The two molecules must be isomorphic
    Arguments:
    - mol1:     first instance of the molecule of interest
    - mol2:     second instance the molecule of interest
    - confId1:  conformer index for mol1 (default: first conformer)
    - useWeights: flag for using torsion weights in the TFD calculation
    - maxDev:   maximal deviation used for normalization

    ‘equal’: all torsions are normalized using 180.0 (default)
    ‘spec’:  each torsion is normalized using its specific

    maximal deviation as given in the paper

    symmRadius: radius used for calculating the atom invariants(default: 2)

    ignoreColinearBonds: if True (default), single bonds adjacent totriple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: TFD value"""
    ...

def GetTFDMatrix(
    self, mol, useWeights=True, maxDev="equal", symmRadius=2, ignoreColinearBonds=True
):
    """
    Wrapper to calculate the matrix of TFD values for the
    conformers of a molecule.
    Arguments:
    - mol:      the molecule of interest
    - useWeights: flag for using torsion weights in the TFD calculation
    - maxDev:   maximal deviation used for normalization

    ‘equal’: all torsions are normalized using 180.0 (default)
    ‘spec’:  each torsion is normalized using its specific

    maximal deviation as given in the paper

    symmRadius: radius used for calculating the atom invariants(default: 2)

    ignoreColinearBonds: if True (default), single bonds adjacent totriple bonds are ignored
    if False, alternative not-covalently bound
    atoms are used to define the torsion

    Return: matrix of TFD values
    Note that the returned matrix is symmetrical, i.e. it is the
    lower half of the matrix, e.g. for 5 conformers:
    matrix = [ a,

    b, c,
    d, e, f,
    g, h, i, j]"""
    ...