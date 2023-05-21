"""
rdkit.Chem.rdDistGeom module¶
Module containing functions to compute atomic coordinates in 3D using distance geometry
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.rdBase import _vecti

CHECK_CHIRAL_CENTERS: EmbedFailureCauses
CHECK_TETRAHEDRAL_CENTERS: EmbedFailureCauses
ETK_MINIMIZATION: EmbedFailureCauses
FINAL_CENTER_IN_VOLUME: EmbedFailureCauses
FINAL_CHIRAL_BOUNDS: EmbedFailureCauses
FIRST_MINIMIZATION: EmbedFailureCauses
INITIAL_COORDS: EmbedFailureCauses
MINIMIZE_FOURTH_DIMENSION: EmbedFailureCauses

class EmbedFailureCauses(Boost.Python.enum):
    CHECK_CHIRAL_CENTERS: EmbedFailureCauses = ...
    CHECK_TETRAHEDRAL_CENTERS: EmbedFailureCauses = ...
    ETK_MINIMIZATION: EmbedFailureCauses = ...
    FINAL_CENTER_IN_VOLUME: EmbedFailureCauses = ...
    FINAL_CHIRAL_BOUNDS: EmbedFailureCauses = ...
    FIRST_MINIMIZATION: EmbedFailureCauses = ...
    INITIAL_COORDS: EmbedFailureCauses = ...
    MINIMIZE_FOURTH_DIMENSION: EmbedFailureCauses = ...
    names: dict[str, EmbedFailureCauses] = ...
    values: dict[int, EmbedFailureCauses] = ...
    __slots__: ClassVar[tuple] = ...

class EmbedParameters(Boost.Python.instance):
    """
    Parameters controlling embedding

    C++ signature :void __init__(_object*)

    property ETversion¶
    version of the experimental torsion-angle preferences"""

    __instance_size__: ClassVar[int] = ...
    ETversion: Any
    boundsMatForceScaling: Any
    boxSizeMult: Any
    clearConfs: Any
    embedFragmentsSeparately: Any
    enforceChirality: Any
    forceTransAmides: Any
    ignoreSmoothingFailures: Any
    maxIterations: Any
    numThreads: Any
    numZeroFail: Any
    onlyHeavyAtomsForRMS: Any
    optimizerForceTol: Any
    pruneRmsThresh: Any
    randNegEig: Any
    randomSeed: Any
    trackFailures: Any
    useBasicKnowledge: Any
    useExpTorsionAnglePrefs: Any
    useMacrocycleTorsions: Any
    useRandomCoords: Any
    useSmallRingTorsions: Any
    useSymmetryForPruning: Any
    verbose: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetFailureCounts(self, arg1: EmbedParameters) -> tuple:
        """
        returns the counts of eacu

        C++ signature :boost::python::tuple GetFailureCounts(RDKit::DGeomHelpers::EmbedParameters*)
        """
        ...
    def SetBoundsMat(self, arg1: EmbedParameters, arg2: AtomPairsParameters) -> None:
        """
        set the distance-bounds matrix to be used (no triangle smoothing will be done on this) from a Numpy array

        C++ signature :void SetBoundsMat(RDKit::DGeomHelpers::EmbedParameters*,boost::python::api::object)
        """
        ...
    def SetCPCI(self, arg1: EmbedParameters, arg2: dict) -> None:
        """
        set the customised pairwise Columb-like interaction to atom pairs.used during structural minimisation stage

        C++ signature :void SetCPCI(RDKit::DGeomHelpers::EmbedParameters*,boost::python::dict {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def ETDG(self) -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETDG method.

    C++ signature :RDKit::DGeomHelpers::EmbedParameters* ETDG()"""
    ...

def ETKDG(self) -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETKDG method - version 1.

    C++ signature :RDKit::DGeomHelpers::EmbedParameters* ETKDG()"""
    ...

def ETKDGv2(self) -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETKDG method - version 2.

    C++ signature :RDKit::DGeomHelpers::EmbedParameters* ETKDGv2()"""
    ...

def ETKDGv3(self) -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETKDG method - version 3 (macrocycles).

    C++ signature :RDKit::DGeomHelpers::EmbedParameters* ETKDGv3()"""
    ...

@overload
def EmbedMolecule(
    self,
    mol: Mol,
    maxAttempts: int = 0,
    randomSeed: int = -1,
    clearConfs: bool = True,
    useRandomCoords: bool = False,
    boxSizeMult: float = 2.0,
    randNegEig: bool = True,
    numZeroFail: int = 1,
    coordMap: dict = {},
    forceTol: float = 0.001,
    ignoreSmoothingFailures: bool = False,
    enforceChirality: bool = True,
    useExpTorsionAnglePrefs: bool = True,
    useBasicKnowledge: bool = True,
    printExpTorsionAngles: bool = False,
    useSmallRingTorsions: bool = False,
    useMacrocycleTorsions: bool = False,
    ETversion: int = 1,
) -> int:
    """

    Use distance geometry to obtain intial coordinates for a molecule
    ARGUMENTS:

    mol : the molecule of interest
    params : an EmbedParameters object

    RETURNS:

    ID of the new conformation added to the molecule

    C++ signature :int EmbedMolecule(RDKit::ROMol {lvalue},RDKit::DGeomHelpers::EmbedParameters {lvalue})
    """
    ...

@overload
def EmbedMolecule(self, mol: Mol, params: EmbedParameters) -> int: ...
@overload
def EmbedMultipleConfs(
    self,
    mol: Mol,
    numConfs: int = 10,
    maxAttempts: int = 0,
    randomSeed: int = -1,
    clearConfs: bool = True,
    useRandomCoords: bool = False,
    boxSizeMult: float = 2.0,
    randNegEig: bool = True,
    numZeroFail: int = 1,
    pruneRmsThresh: float = -1.0,
    coordMap: dict = {},
    forceTol: float = 0.001,
    ignoreSmoothingFailures: bool = False,
    enforceChirality: bool = True,
    numThreads: int = 1,
    useExpTorsionAnglePrefs: bool = True,
    useBasicKnowledge: bool = True,
    printExpTorsionAngles: bool = False,
    useSmallRingTorsions: bool = False,
    useMacrocycleTorsions: bool = False,
    ETversion: int = 1,
) -> _vecti:
    """

    Use distance geometry to obtain multiple sets of coordinates for a molecule
    ARGUMENTS:

    mol : the molecule of interest
    numConfs : the number of conformers to generate
    params : an EmbedParameters object

    RETURNS:

    List of new conformation IDs

    C++ signature :std::vector<int, std::allocator<int> > EmbedMultipleConfs(RDKit::ROMol {lvalue},unsigned int,RDKit::DGeomHelpers::EmbedParameters {lvalue})
    """
    ...

@overload
def EmbedMultipleConfs(
    self, mol: Mol, numConfs: int, params: EmbedParameters
) -> _vecti: ...
@overload
def GetExperimentalTorsions(
    self,
    mol: Mol,
    useExpTorsionAnglePrefs: bool = True,
    useSmallRingTorsions: bool = False,
    useMacrocycleTorsions: bool = True,
    useBasicKnowledge: bool = True,
    ETversion: int = 2,
    printExpTorsionAngles: bool = False,
) -> tuple:
    """
    returns information about the bonds corresponding to experimental torsions

    C++ signature :boost::python::tuple GetExperimentalTorsions(RDKit::ROMol,RDKit::DGeomHelpers::EmbedParameters)
    """
    ...

@overload
def GetExperimentalTorsions(self, mol: Mol, embedParams: EmbedParameters) -> tuple: ...
def GetMoleculeBoundsMatrix(
    self,
    mol: Mol,
    set15bounds: bool = True,
    scaleVDW: bool = False,
    doTriangleSmoothing: bool = True,
    useMacrocycle14config: bool = False,
) -> object:
    """
    Returns the distance bounds matrix for a molecule

    ARGUMENTS:

    mol : the molecule of interest

    set15boundsset bounds for 1-5 atom distances based on topology (otherwise stop at 1-4s)

    scaleVDWscale down the sum of VDW radii when setting the lower bounds for atoms less than 5 bonds apart

    doTriangleSmoothingrun triangle smoothing on the bounds matrix before returning it

    RETURNS:

    the bounds matrix as a Numeric array with lower bounds in
    the lower triangle and upper bounds in the upper triangle

    C++ signature :_object* GetMoleculeBoundsMatrix(RDKit::ROMol {lvalue} [,bool=True [,bool=False [,bool=True [,bool=False]]]])
    """
    ...

def KDG(self) -> EmbedParameters:
    """
    Returns an EmbedParameters object for the KDG method.

    C++ signature :RDKit::DGeomHelpers::EmbedParameters* KDG()"""
    ...

def srETKDGv3(self) -> EmbedParameters:
    """
    Returns an EmbedParameters object for the ETKDG method - version 3 (small rings).

    C++ signature :RDKit::DGeomHelpers::EmbedParameters* srETKDGv3()"""
    ...
