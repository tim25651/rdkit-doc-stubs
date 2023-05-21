"""
rdkit.Chem.rdDepictor module¶
Module containing the functionality to compute 2D coordinates for a molecule
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class UsingCoordGen(Boost.Python.instance):
    """
    Context manager to temporarily set CoordGen library preference in RDKit depiction.
    Constructor

    C++ signature :void __init__(_object*,bool)"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __enter__(cls, RDDepict) -> Any: ...
    @classmethod
    def __exit__(cls, type, value, traceback) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def AddRingSystemTemplates(self, templatePath: str) -> None:
    """
    Adds the ring system templates from the specified file to be used in 2D coordinate generation. If there are duplicates, the most recently added template will be used. Each template must be a single line in the file represented using CXSMILES, and the structure should be a single ring system. Throws a DepictException if any templates are invalid.

    C++ signature :void AddRingSystemTemplates(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def Compute2DCoords(
    self,
    mol: Mol,
    canonOrient: bool = True,
    clearConfs: bool = True,
    coordMap: dict = {},
    nFlipsPerSample: int = 0,
    nSample: int = 0,
    sampleSeed: int = 0,
    permuteDeg4Nodes: bool = False,
    bondLength: float = -1.0,
    forceRDKit: bool = False,
    useRingTemplates: bool = False,
) -> int:
    """
    Compute 2D coordinates for a molecule. The resulting coordinates are stored on each atom of the molecule
    ARGUMENTS:

    mol - the molecule of interest
    canonOrient - orient the molecule in a canonical way
    clearConfs - if true, all existing conformations on the molecule

    will be cleared

    coordMap - a dictionary mapping atom Ids -> Point2D objects with starting coordinates for atoms that should
    have their positions locked.

    nFlipsPerSample - number of rotatable bonds that areflipped at random at a time.

    nSample - Number of random samplings of rotatable bonds.
    sampleSeed - seed for the random sampling process.
    permuteDeg4Nodes - allow permutation of bonds at a degree 4

    node during the sampling process

    bondLength - change the default bond length for depiction
    forceRDKit - use RDKit to generate coordinates even if

    preferCoordGen is set to true

    useRingTemplates - use templates to generate coordinates of complexring systems

    RETURNS:

    ID of the conformation added to the molecule

    C++ signature :unsigned int Compute2DCoords(RDKit::ROMol {lvalue} [,bool=True [,bool=True [,boost::python::dict {lvalue}={} [,unsigned int=0 [,unsigned int=0 [,int=0 [,bool=False [,double=-1.0 [,bool=False [,bool=False]]]]]]]]]])
    """
    ...

def Compute2DCoordsMimicDistmat(
    self,
    mol: Mol,
    distMat: AtomPairsParameters,
    canonOrient: bool = False,
    clearConfs: bool = True,
    weightDistMat: float = 0.5,
    nFlipsPerSample: int = 3,
    nSample: int = 100,
    sampleSeed: int = 100,
    permuteDeg4Nodes: bool = True,
    bondLength: float = -1.0,
    forceRDKit: bool = False,
) -> int:
    """
    Compute 2D coordinates for a molecule such that the inter-atom distances mimic those in a user-provided
    distance matrix.
    The resulting coordinates are stored on each atom of the molecule
    ARGUMENTS:

    mol - the molecule of interest
    distMat - distance matrix that we want the 2D structure to mimic
    canonOrient - orient the molecule in a canonical way
    clearConfs - if true, all existing conformations on the molecule

    will be cleared

    weightDistMat - weight assigned in the cost function to mimickingthe distance matrix.
    This must be between (0.0,1.0). (1.0-weightDistMat)
    is then the weight assigned to improving
    the density of the 2D structure i.e. try to
    make it spread out

    nFlipsPerSample - number of rotatable bonds that areflipped at random at a time.

    nSample - Number of random samplings of rotatable bonds.
    sampleSeed - seed for the random sampling process.
    permuteDeg4Nodes - allow permutation of bonds at a degree 4

    node during the sampling process

    bondLength - change the default bond length for depiction
    forceRDKit - use RDKit to generate coordinates even if

    preferCoordGen is set to true

    RETURNS:

    ID of the conformation added to the molecule

    C++ signature :unsigned int Compute2DCoordsMimicDistmat(RDKit::ROMol {lvalue},boost::python::api::object [,bool=False [,bool=True [,double=0.5 [,unsigned int=3 [,unsigned int=100 [,int=100 [,bool=True [,double=-1.0 [,bool=False]]]]]]]]])
    """
    ...

@overload
def GenerateDepictionMatching2DStructure(
    self,
    mol: Mol,
    reference: Mol,
    confId: int = -1,
    refPatt: AtomPairsParameters = None,
    acceptFailure: bool = False,
    forceRDKit: bool = False,
    allowRGroups: bool = False,
) -> tuple:
    """

    Generate a depiction for a molecule where a piece of the molecule is constrained to have the same coordinates as a reference.
    This is useful for, for example, generating depictions of SAR data
    sets so that the cores of the molecules are all oriented the same way.
    ARGUMENTS:

    mol -    the molecule to be aligned, this will come back with a single conformer.

    reference -    a molecule with the reference atoms to align to; this should have a depiction.

    atomMap -      a sequence of (queryAtomIdx, molAtomIdx) pairs that will be used to generate the atom mapping between the molecule
    and the reference.

    confId -       (optional) the id of the reference conformation to use
    forceRDKit -   (optional) use RDKit to generate coordinates even if

    preferCoordGen is set to true

    C++ signature :void GenerateDepictionMatching2DStructure(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue},boost::python::api::object [,int=-1 [,bool=False]])
    """
    ...

@overload
def GenerateDepictionMatching2DStructure(
    self,
    mol: Mol,
    reference: Mol,
    atomMap: AtomPairsParameters,
    confId: int = -1,
    forceRDKit: bool = False,
) -> None: ...
def GenerateDepictionMatching3DStructure(
    self,
    mol: Mol,
    reference: Mol,
    confId: int = -1,
    refPatt: AtomPairsParameters = None,
    acceptFailure: bool = False,
    forceRDKit: bool = False,
) -> None:
    """
    Generate a depiction for a molecule where a piece of the molecule is constrained to have coordinates similar to those of a 3D reference
    structure.
    ARGUMENTS:

    mol -    the molecule to be aligned, this will come back with a single conformer containing the 2D coordinates.

    reference -    a molecule with the reference atoms to align to. By default this should be the same as mol, but with
    3D coordinates

    confId -       (optional) the id of the reference conformation to use
    referencePattern -  (optional) a query molecule to map a subset of

    the reference onto the mol, so that only some of the
    atoms are aligned.

    acceptFailure - (optional) if True, standard depictions will be generated for molecules that don’t match the reference or the
    referencePattern; if False, throws a DepictException.

    forceRDKit -    (optional) use RDKit to generate coordinates even if preferCoordGen is set to true

    C++ signature :void GenerateDepictionMatching3DStructure(RDKit::ROMol {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,boost::python::api::object=None [,bool=False [,bool=False]]]])
    """
    ...

def GetPreferCoordGen(self) -> bool:
    """
    Return whether or not the CoordGen library is used for coordinate generation in the RDKit depiction library.

    C++ signature :bool GetPreferCoordGen()"""
    ...

def IsCoordGenSupportAvailable(self) -> bool:
    """
    Returns whether RDKit was built with CoordGen support.

    C++ signature :bool IsCoordGenSupportAvailable()"""
    ...

def LoadDefaultRingSystemTemplates(self) -> None:
    """
    Loads the default ring system templates and removes existing ones, if present.

    C++ signature :void LoadDefaultRingSystemTemplates()"""
    ...

def NormalizeDepiction(
    self, mol: Mol, confId: int = -1, canonicalize: int = 1, scaleFactor: float = -1.0
) -> float:
    """
    Normalizes the 2D depiction.
    If canonicalize is != 0, the depiction is subjected to a canonical
    transformation such that its main axis is aligned along the X axis
    (canonicalize >0, the default) or the Y axis (canonicalize <0).
    If canonicalize is 0, no canonicalization takes place.
    If scaleFactor is <0.0 (the default) the depiction is scaled such
    that bond lengths conform to RDKit standards. The applied scaling
    factor is returned.
    ARGUMENTS:
    mol          - the molecule to be normalized
    confId       - (optional) the id of the reference conformation to use
    canonicalize - (optional) if != 0, a canonical transformation is

    applied: if >0 (the default), the main molecule axis is
    aligned to the X axis, if <0 to the Y axis.
    If 0, no canonical transformation is applied.

    scaleFactor  - (optional) if >0.0, the scaling factor to apply. The default(-1.0) means that the depiction is automatically scaled
    such that bond lengths are the standard RDKit ones.

    RETURNS: the applied scaling factor.

    C++ signature :double NormalizeDepiction(RDKit::ROMol {lvalue} [,int=-1 [,int=1 [,double=-1.0]]])
    """
    ...

def SetPreferCoordGen(self, val: bool) -> None:
    """
    Sets whether or not the CoordGen library should be preferred to the RDKit depiction library.

    C++ signature :void SetPreferCoordGen(bool)"""
    ...

def SetRingSystemTemplates(self, templatePath: str) -> None:
    """
    Loads the ring system templates from the specified file to be used in 2D coordinate generation. Each template must be a single line in the file represented using CXSMILES, and the structure should be a single ring system. Throws a DepictException if any templates are invalid.

    C++ signature :void SetRingSystemTemplates(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def StraightenDepiction(
    self, mol: Mol, confId: int = -1, minimizeRotation: bool = False
) -> None:
    """
    Rotate the 2D depiction such that the majority of bonds have a30-degree angle with the X axis.
    ARGUMENTS:
    mol              - the molecule to be rotated.
    confId           - (optional) the id of the reference conformation to use.
    minimizeRotation - (optional) if False (the default), the molecule

    is rotated such that the majority of bonds have an angle
    with the X axis of 30 or 90 degrees. If True, the minimum
    rotation is applied such that the majority of bonds have
    an angle with the X axis of 0, 30, 60, or 90 degrees,
    with the goal of altering the initial orientation as
    little as possible .

    C++ signature :void StraightenDepiction(RDKit::ROMol {lvalue} [,int=-1 [,bool=False]])
    """
    ...
