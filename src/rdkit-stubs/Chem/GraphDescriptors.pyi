"""
rdkit.Chem.GraphDescriptors module¶
Calculation of topological/topochemical descriptors.
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Graphs as Graphs
from rdkit.Chem import rdchem as rdchem
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.ML.InfoTheory import entropy as entropy

ptable: Incomplete
hallKierAlphas: Incomplete

def Ipc(self, mol, avg=False, dMat=None, forceDMat=False):
    """
    This returns the information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.
    ‘avg = True’ returns the information content divided by the total population.
    From Eq 6 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)"""
    ...

def AvgIpc(self, mol, dMat=None, forceDMat=False):
    """
    This returns the average information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.
    From Eq 7 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)"""
    ...

def BalabanJ(self, mol, dMat=None, forceDMat=0):
    """
    Calculate Balaban’s J value for a molecule
    Arguments

    mol: a molecule
    dMat: (optional) a distance/adjacency matrix for the molecule, if this
    is not provide, one will be calculated
    forceDMat: (optional) if this is set, the distance/adjacency matrix
    will be recalculated regardless of whether or not _dMat_ is provided
    or the molecule already has one

    Returns

    a float containing the J value

    We follow the notation of Balaban’s paper:Chem. Phys. Lett. vol 89, 399-404, (1982)
    """
    ...

def BertzCT(self, mol, cutoff=100, dMat=None, forceDMat=1):
    """
    A topological index meant to quantify “complexity” of molecules.
    Consists of a sum of two terms, one representing the complexity
    of the bonding, the other representing the complexity of the
    distribution of heteroatoms.
    From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)
    “cutoff” is an integer value used to limit the computational
    expense.  A cutoff value tells the program to consider vertices
    topologically identical if their distance vectors (sets of
    distances to all other vertices) are equal out to the “cutoff”th
    nearest-neighbor.

    NOTE  The original implementation had the following comment:
    > this implementation treats aromatic rings as the
    > corresponding Kekule structure with alternating bonds,
    > for purposes of counting “connections”.

    Upon further thought, this is the WRONG thing to do.  Itresults in the possibility of a molecule giving two different
    CT values depending on the kekulization.  For example, in the
    old implementation, these two SMILES:

    CC2=CN=C1C3=C(C(C)=C(C=N3)C)C=CC1=C2C
    CC3=CN=C2C1=NC=C(C)C(C)=C1C=CC2=C3C

    which correspond to differentk kekule forms, yield different
    values.

    The new implementation uses consistent (aromatic) bond ordersfor aromatic bonds.

    THIS MEANS THAT THIS IMPLEMENTATION IS NOT BACKWARDS COMPATIBLE.
    Any molecule containing aromatic rings will yield different
    values with this implementation.  The new behavior is the correct
    one, so we’re going to live with the breakage.

    NOTE this barfs if the molecule contains a second (ornth) fragment that is one atom.
    """
    ...

def HallKierAlpha(self, x): ...
def Kappa1(self, x): ...
def Kappa2(self, x): ...
def Kappa3(self, x): ...
def Chi0(self, mol):
    """
    From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)"""
    ...

def Chi0n(self, x): ...
def Chi0v(self, x): ...
def Chi1(self, mol):
    """
    From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)"""
    ...

def Chi1n(self, x): ...
def Chi0v(self, x): ...
def Chi1v(self, x): ...
def Chi2n(self, x): ...
def Chi2v(self, x): ...
def Chi3n(self, x): ...
def Chi3v(self, x): ...
def Chi4n(self, x): ...
def Chi4v(self, x): ...
def ChiNn_(self, x, y): ...
def ChiNv_(self, x, y): ...
def HallKierAlpha(self, x): ...
def Ipc(self, mol, avg=False, dMat=None, forceDMat=False):
    """
    This returns the information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.
    ‘avg = True’ returns the information content divided by the total population.
    From Eq 6 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)"""
    ...

def Kappa1(self, x): ...
def Kappa2(self, x): ...
def Kappa3(self, x): ...
def Chi0n(self, x): ...
def Chi1n(self, x): ...
def Chi2n(self, x): ...
def Chi3n(self, x): ...
def Chi4n(self, x): ...
def ChiNn_(self, x, y): ...
def BalabanJ(self, mol, dMat=None, forceDMat=0):
    """
    Calculate Balaban’s J value for a molecule
    Arguments

    mol: a molecule
    dMat: (optional) a distance/adjacency matrix for the molecule, if this
    is not provide, one will be calculated
    forceDMat: (optional) if this is set, the distance/adjacency matrix
    will be recalculated regardless of whether or not _dMat_ is provided
    or the molecule already has one

    Returns

    a float containing the J value

    We follow the notation of Balaban’s paper:Chem. Phys. Lett. vol 89, 399-404, (1982)
    """
    ...

def BertzCT(self, mol, cutoff=100, dMat=None, forceDMat=1):
    """
    A topological index meant to quantify “complexity” of molecules.
    Consists of a sum of two terms, one representing the complexity
    of the bonding, the other representing the complexity of the
    distribution of heteroatoms.
    From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)
    “cutoff” is an integer value used to limit the computational
    expense.  A cutoff value tells the program to consider vertices
    topologically identical if their distance vectors (sets of
    distances to all other vertices) are equal out to the “cutoff”th
    nearest-neighbor.

    NOTE  The original implementation had the following comment:
    > this implementation treats aromatic rings as the
    > corresponding Kekule structure with alternating bonds,
    > for purposes of counting “connections”.

    Upon further thought, this is the WRONG thing to do.  Itresults in the possibility of a molecule giving two different
    CT values depending on the kekulization.  For example, in the
    old implementation, these two SMILES:

    CC2=CN=C1C3=C(C(C)=C(C=N3)C)C=CC1=C2C
    CC3=CN=C2C1=NC=C(C)C(C)=C1C=CC2=C3C

    which correspond to differentk kekule forms, yield different
    values.

    The new implementation uses consistent (aromatic) bond ordersfor aromatic bonds.

    THIS MEANS THAT THIS IMPLEMENTATION IS NOT BACKWARDS COMPATIBLE.
    Any molecule containing aromatic rings will yield different
    values with this implementation.  The new behavior is the correct
    one, so we’re going to live with the breakage.

    NOTE this barfs if the molecule contains a second (ornth) fragment that is one atom.
    """
    ...
