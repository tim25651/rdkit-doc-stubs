"""
rdkit.Chem.rdmolops module¶
Module containing RDKit functionality for manipulating molecules.
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Bond, Conformer, Mol, MolBundle, StereoGroupType
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Chem.rdmolops import _vectN5RDKit9Chirality10StereoInfoE
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, ULongSparseIntVect
from rdkit.Geometry import rdGeometry
from rdkit.Geometry.rdGeometry import Point3D
from rdkit.rdBase import _listSt6vectorIiSaIiEE, _vecti, _vectSt6vectorIiSaIiEE

ADJUST_IGNOREALL: AdjustQueryWhichFlags
ADJUST_IGNORECHAINS: AdjustQueryWhichFlags
ADJUST_IGNOREDUMMIES: AdjustQueryWhichFlags
ADJUST_IGNORENONDUMMIES: AdjustQueryWhichFlags
ADJUST_IGNORENONE: AdjustQueryWhichFlags
ADJUST_IGNORERINGS: AdjustQueryWhichFlags
AROMATICITY_CUSTOM: AromaticityModel
AROMATICITY_DEFAULT: AromaticityModel
AROMATICITY_MDL: AromaticityModel
AROMATICITY_RDKIT: AromaticityModel
AROMATICITY_SIMPLE: AromaticityModel
LayeredFingerprint_substructLayers: int
SANITIZE_ADJUSTHS: SanitizeFlags
SANITIZE_ALL: SanitizeFlags
SANITIZE_CLEANUP: SanitizeFlags
SANITIZE_CLEANUPCHIRALITY: SanitizeFlags
SANITIZE_FINDRADICALS: SanitizeFlags
SANITIZE_KEKULIZE: SanitizeFlags
SANITIZE_NONE: SanitizeFlags
SANITIZE_PROPERTIES: SanitizeFlags
SANITIZE_SETAROMATICITY: SanitizeFlags
SANITIZE_SETCONJUGATION: SanitizeFlags
SANITIZE_SETHYBRIDIZATION: SanitizeFlags
SANITIZE_SYMMRINGS: SanitizeFlags
_LayeredFingerprint_version: str
_PatternFingerprint_version: str
_RDKFingerprint_version: str

class AdjustQueryParameters(Boost.Python.instance):
    """
    Parameters controlling which components of the query atoms/bonds are adjusted.

    Note that some of the options here are either directly contradictory or makeno sense when combined with each other. We generally assume that client code
    is doing something sensible and don’t attempt to detect possible conflicts or
    problems.

    A note on the flags controlling which atoms/bonds are modified: These generally limit the set of atoms/bonds to be modified.
    For example:

    ADJUST_IGNORERINGS atoms/bonds in rings will not be modified.
    ADJUST_IGNORENONE causes all atoms/bonds to be modified
    ADJUST_IGNOREALL no atoms/bonds will be modified

    Some of the options obviously make no sense for bonds

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...
    adjustConjugatedFiveRings: Any
    adjustDegree: Any
    adjustDegreeFlags: Any
    adjustHeavyDegree: Any
    adjustHeavyDegreeFlags: Any
    adjustRingChain: Any
    adjustRingChainFlags: Any
    adjustRingCount: Any
    adjustRingCountFlags: Any
    adjustSingleBondsBetweenAromaticAtoms: Any
    adjustSingleBondsToDegreeOneNeighbors: Any
    aromatizeIfPossible: Any
    makeAtomsGeneric: Any
    makeAtomsGenericFlags: Any
    makeBondsGeneric: Any
    makeBondsGenericFlags: Any
    makeDummiesQueries: Any
    setMDLFiveRingAromaticity: Any
    useStereoCareForBonds: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def NoAdjustments(self) -> AdjustQueryParameters:
        """
        Returns an AdjustQueryParameters object with all parameters set to false

        C++ signature :RDKit::MolOps::AdjustQueryParameters NoAdjustments()"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AdjustQueryWhichFlags(Boost.Python.enum):
    ADJUST_IGNOREALL: AdjustQueryWhichFlags = ...
    ADJUST_IGNORECHAINS: AdjustQueryWhichFlags = ...
    ADJUST_IGNOREDUMMIES: AdjustQueryWhichFlags = ...
    ADJUST_IGNORENONDUMMIES: AdjustQueryWhichFlags = ...
    ADJUST_IGNORENONE: AdjustQueryWhichFlags = ...
    ADJUST_IGNORERINGS: AdjustQueryWhichFlags = ...
    names: dict[str, AdjustQueryWhichFlags] = ...
    values: dict[int, AdjustQueryWhichFlags] = ...
    __slots__: ClassVar[tuple] = ...

class AromaticityModel(Boost.Python.enum):
    AROMATICITY_CUSTOM: AromaticityModel = ...
    AROMATICITY_DEFAULT: AromaticityModel = ...
    AROMATICITY_MDL: AromaticityModel = ...
    AROMATICITY_RDKIT: AromaticityModel = ...
    AROMATICITY_SIMPLE: AromaticityModel = ...
    names: dict[str, AromaticityModel] = ...
    values: dict[int, AromaticityModel] = ...
    __slots__: ClassVar[tuple] = ...

class BondWedgingParameters(Boost.Python.instance):
    """
    Parameters controlling how bond wedging is done.

    C++ signature :void __init__(_object*)

    property wedgeTwoBondsIfPossible¶
    If this is enabled then two bonds will be wedged at chiral
    centers subject to the following constraints:

    ring bonds will not be wedged
    bonds to chiral centers will not be wedged

    bonds separated by more than 120 degrees will not bewedged"""

    ...
    __instance_size__: ClassVar[int] = ...
    wedgeTwoBondsIfPossible: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolzipLabel(Boost.Python.enum):
    AtomMapNumber: MolzipLabel = ...
    AtomType: MolzipLabel = ...
    FragmentOnBonds: MolzipLabel = ...
    Isotope: MolzipLabel = ...
    names: dict[str, MolzipLabel] = ...
    values: dict[int, MolzipLabel] = ...
    __slots__: ClassVar[tuple] = ...

class MolzipParams(Boost.Python.instance):
    """
    Parameters controllnig how to zip molecules together

    OPTIONS:label : set the MolzipLabel option [default MolzipLabel.AtomMapNumber]

    MolzipLabel.AtomMapNumber: atom maps are on dummy atoms, zip together the correspondingattaced atoms, i.e.  zip ‘C[:1]’ ‘N[:1]’ results in ‘CN’

    MolzipLabel.Isotope: isotope labels are on dummy atoms, zip together the correspondingattaced atoms, i.e.  zip ‘C[1*]’ ‘N[1*]’ results in ‘CN’

    MolzipLabel.FragmentOnBonds: zip together molecules generated by fragment on bonds.Note the atom indices cannot change or be reorderd from the output of fragmentOnBonds

    MolzipLabel.AtomTypes: choose the atom types to act as matching dummy atoms.i.e.  ‘C[V]’ and ‘N[Xe]’ with atoms pairs [(‘V’, ‘Xe’)] results in ‘CN’

    C++ signature :void __init__(_object*)

    property enforceValenceRules¶
    If true (default) enforce valences after zipping
    Setting this to false allows assembling chemically incorrect fragments.

    property generateCoordinates¶
    If true will add depiction coordinates to input molecules and
    zipped molecule (for molzipFragments only)

    property label¶
    Set the atom labelling system to zip together"""

    __instance_size__: ClassVar[int] = ...
    enforceValenceRules: Any
    generateCoordinates: Any
    label: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def setAtomSymbols(self, arg1: MolzipParams, arg2: AtomPairsParameters) -> None:
        """
        Set the atom symbols used to zip mols together when using AtomType labeling

        C++ signature :void setAtomSymbols(RDKit::MolzipParams {lvalue},boost::python::api::object)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RemoveHsParameters(Boost.Python.instance):
    """
    Parameters controlling which Hs are removed.

    C++ signature :void __init__(_object*)

    property removeAndTrackIsotopes¶
    hydrogens with non-default isotopes and store them in the _isotopicHs atom property such that AddHs() can add the same isotope at a later stage

    property removeDefiningBondStereo¶
    hydrogens defining bond stereochemistry

    property removeDegreeZero¶
    hydrogens that have no bonds

    property removeDummyNeighbors¶
    hydrogens with at least one dummy-atom neighbor

    property removeHigherDegrees¶
    hydrogens with two (or more) bonds

    property removeHydrides¶
    hydrogens with formal charge -1

    property removeInSGroups¶
    hydrogens involved in SubstanceGroups

    property removeIsotopes¶
    hydrogens with non-default isotopes

    property removeMapped¶
    mapped hydrogens

    property removeNonimplicit¶
    DEPRECATED

    property removeNontetrahedralNeighbors¶
    hydrogens with neighbors that have non-tetrahedral stereochemistry

    property removeOnlyHNeighbors¶
    hydrogens with bonds only to other hydrogens

    property removeWithQuery¶
    hydrogens with queries defined

    property removeWithWedgedBond¶
    hydrogens with wedged bonds to them

    property showWarnings¶
    display warning messages for some classes of removed Hs

    property updateExplicitCount¶
    DEPRECATED"""

    ...
    __instance_size__: ClassVar[int] = ...
    removeAndTrackIsotopes: Any
    removeDefiningBondStereo: Any
    removeDegreeZero: Any
    removeDummyNeighbors: Any
    removeHigherDegrees: Any
    removeHydrides: Any
    removeInSGroups: Any
    removeIsotopes: Any
    removeMapped: Any
    removeNonimplicit: Any
    removeNontetrahedralNeighbors: Any
    removeOnlyHNeighbors: Any
    removeWithQuery: Any
    removeWithWedgedBond: Any
    showWarnings: Any
    updateExplicitCount: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class SanitizeFlags(Boost.Python.enum):
    SANITIZE_ADJUSTHS: SanitizeFlags = ...
    SANITIZE_ALL: SanitizeFlags = ...
    SANITIZE_CLEANUP: SanitizeFlags = ...
    SANITIZE_CLEANUPCHIRALITY: SanitizeFlags = ...
    SANITIZE_FINDRADICALS: SanitizeFlags = ...
    SANITIZE_KEKULIZE: SanitizeFlags = ...
    SANITIZE_NONE: SanitizeFlags = ...
    SANITIZE_PROPERTIES: SanitizeFlags = ...
    SANITIZE_SETAROMATICITY: SanitizeFlags = ...
    SANITIZE_SETCONJUGATION: SanitizeFlags = ...
    SANITIZE_SETHYBRIDIZATION: SanitizeFlags = ...
    SANITIZE_SYMMRINGS: SanitizeFlags = ...
    names: dict[str, SanitizeFlags] = ...
    values: dict[int, SanitizeFlags] = ...
    __slots__: ClassVar[tuple] = ...

class StereoBondThresholds(Boost.Python.instance):
    """
    Constants used to set the thresholds for which single bonds can be made wavy.
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    CHIRAL_ATOM: property[int] = ...
    DBL_BOND_NO_STEREO: property[int] = ...
    DBL_BOND_SPECIFIED_STEREO: property[int] = ...
    DIRECTION_SET: property[int] = ...

class _vectN5RDKit9Chirality10StereoInfoE(Boost.Python.instance):
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def append(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def extend(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __contains__(cls, other) -> Any: ...
    @classmethod
    def __delitem__(cls, other) -> Any: ...
    @classmethod
    def __getitem__(cls, index) -> Any: ...
    @classmethod
    def __iter__(cls, boost, std) -> Any: ...
    @classmethod
    def __len__(cls) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, index, object) -> Any: ...

def AddHs(
    self,
    mol: Mol,
    explicitOnly: bool = False,
    addCoords: bool = False,
    onlyOnAtoms: AtomPairsParameters = None,
    addResidueInfo: bool = False,
) -> Mol:
    """
    Adds hydrogens to the graph of a molecule.

    ARGUMENTS:

    mol: the molecule to be modified
    explicitOnly: (optional) if this toggle is set, only explicit Hs will
    be added to the molecule.  Default value is 0 (add implicit and explicit Hs).
    addCoords: (optional) if this toggle is set, The Hs will have 3D coordinates
    set.  Default value is 0 (no 3D coords).
    onlyOnAtoms: (optional) if this sequence is provided, only these atoms will be
    considered to have Hs added to them
    addResidueInfo: (optional) if this is true, add residue info to
    hydrogen atoms (useful for PDB files).

    RETURNS: a new molecule with added Hs
    NOTES:

    The original molecule is not modified.
    Much of the code assumes that Hs are not included in the molecular
    topology, so be very careful with the molecule that comes back from
    this function.

    C++ signature :RDKit::ROMol* AddHs(RDKit::ROMol [,bool=False [,bool=False [,boost::python::api::object=None [,bool=False]]]])
    """
    ...

def AddRecursiveQuery(
    self, mol: Mol, query: Mol, atomIdx: int, preserveExistingQuery: bool = True
) -> None:
    """
    Adds a recursive query to an atom

    ARGUMENTS:

    mol: the molecule to be modified
    query: the molecule to be used as the recursive query (this will be copied)
    atomIdx: the atom to modify
    preserveExistingQuery: (optional) if this is set, existing query information on the atom will be preserved

    RETURNS: None

    C++ signature :void AddRecursiveQuery(RDKit::ROMol {lvalue},RDKit::ROMol,unsigned int [,bool=True])
    """
    ...

def AddWavyBondsForStereoAny(
    self, mol: Mol, clearDoubleBondFlags: bool = True, addWhenImpossible: int = 1000
) -> None:
    """
    set wavy bonds around double bonds with STEREOANY stereo
    ARGUMENTS :
    molecule : the molecule to updaten -
    conformer : the conformer to use to determine wedge direction

    C++ signature :void AddWavyBondsForStereoAny(RDKit::ROMol {lvalue} [,bool=True [,unsigned int=1000]])
    """
    ...

def AdjustQueryProperties(self, mol: Mol, params: AtomPairsParameters = None) -> Mol:
    """
    Returns a new molecule where the query properties of atoms have been modified.

    C++ signature :RDKit::ROMol* AdjustQueryProperties(RDKit::ROMol [,boost::python::api::object=None])
    """
    ...

def AssignAtomChiralTagsFromMolParity(
    self, mol: Mol, replaceExistingTags: bool = True
) -> None:
    """
    Sets the chiral tags on a molecule’s atoms based onthe molParity atom property.
    ARGUMENTS:

    mol: the molecule to use
    replaceExistingTags: if True, existing stereochemistry information will be cleared

    before running the calculation.

    C++ signature :void AssignAtomChiralTagsFromMolParity(RDKit::ROMol {lvalue} [,bool=True])
    """
    ...

def AssignAtomChiralTagsFromStructure(
    self, mol: Mol, confId: int = -1, replaceExistingTags: bool = True
) -> None:
    """
    Sets the chiral tags on a molecule’s atoms based ona 3D conformation.
    NOTE that this does not check to see if atoms are chiral centers (i.e. all
    substituents are different), it merely sets the chiral type flags based on the
    coordinates and atom ordering. Use AssignStereochemistryFrom3D() if you
    want chiral flags only on actual stereocenters.
    ARGUMENTS:

    mol: the molecule to use
    confId: the conformer id to use, -1 for the default
    replaceExistingTags: if True, existing stereochemistry information will be cleared

    before running the calculation.

    C++ signature :void AssignAtomChiralTagsFromStructure(RDKit::ROMol {lvalue} [,int=-1 [,bool=True]])
    """
    ...

def AssignChiralTypesFromBondDirs(
    self, mol: Mol, confId: int = -1, replaceExistingTags: bool = True
) -> None:
    """
    Uses bond directions to assign ChiralTypes to a molecule’s atoms.

    ARGUMENTS:

    mol: the molecule to use
    confId: (optional) the conformation to use
    replaceExistingTags: (optional) replace any existing information about stereochemistry

    C++ signature :void AssignChiralTypesFromBondDirs(RDKit::ROMol {lvalue} [,int=-1 [,bool=True]])
    """
    ...

def AssignRadicals(self, mol: Mol) -> None:
    """
    Assigns radical counts to atoms

    ARGUMENTS:

    mol: the molecule to use

    NOTES:

    The molecule is modified in place.

    C++ signature :void AssignRadicals(RDKit::ROMol {lvalue})"""
    ...

def AssignStereochemistry(
    self,
    mol: Mol,
    cleanIt: bool = False,
    force: bool = False,
    flagPossibleStereoCenters: bool = False,
) -> None:
    """
    Does the CIP stereochemistry assignment for the molecule’s atoms (R/S) and double bond (Z/E).
    Chiral atoms will have a property ‘_CIPCode’ indicating
    their chiral code.
    ARGUMENTS:

    mol: the molecule to use

    cleanIt: (optional) if provided, any existing values of the property _CIPCodewill be cleared, atoms with a chiral specifier that aren’t

    actually chiral (e.g. atoms with duplicate substituents or only 2 substituents,
    etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with
    STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP
    atom ranks will be marked STEREONONE.

    force: (optional) causes the calculation to be repeated, even if it has already
    been done
    flagPossibleStereoCenters (optional)   set the _ChiralityPossible property on
    atoms that are possible stereocenters

    C++ signature :void AssignStereochemistry(RDKit::ROMol {lvalue} [,bool=False [,bool=False [,bool=False]]])
    """
    ...

def AssignStereochemistryFrom3D(
    self, mol: Mol, confId: int = -1, replaceExistingTags: bool = True
) -> None:
    """
    Uses a conformer (should be 3D) to assign ChiralTypes to a molecule’s atoms
    and stereo flags to its bonds

    ARGUMENTS:

    mol: the molecule to use
    confId: (optional) the conformation to use
    replaceExistingTags: (optional) replace any existing information about stereochemistry

    C++ signature :void AssignStereochemistryFrom3D(RDKit::ROMol {lvalue} [,int=-1 [,bool=True]])
    """
    ...

def Cleanup(self, mol: Mol) -> None:
    """
    cleans up certain common bad functionalities in the molecule

    ARGUMENTS:

    mol: the molecule to use

    NOTES:

    The molecule is modified in place.

    C++ signature :void Cleanup(RDKit::ROMol {lvalue})"""
    ...

def CleanupOrganometallics(self, mol: Mol) -> None:
    """
    cleans up certain common bad functionalities in the organometallic molecule

    Note that this function is experimental and may either change in behavior
    or be replaced with something else in future releases.

    ARGUMENTS :

    mol : the molecule to use

    NOTES :

    The molecule is modified in place.

    C++ signature :void CleanupOrganometallics(RDKit::ROMol {lvalue})"""
    ...

def CombineMols(self, mol1: Mol, mol2: Mol, offset: rdGeometry.Point3D = ...) -> Mol:
    """
    Combine the atoms from two molecules to produce a third

    C++ signature :RDKit::ROMol* CombineMols(RDKit::ROMol,RDKit::ROMol [,RDGeom::Point3D=<rdkit.Geometry.rdGeometry.Point3D object at 0x7fd481d02740>])
    """
    ...

def ConvertGenericQueriesToSubstanceGroups(self, mol: Mol) -> None:
    """
    documentation

    C++ signature :void ConvertGenericQueriesToSubstanceGroups(RDKit::ROMol {lvalue})"""
    ...

def DativeBondsToHaptic(self, mol: Mol) -> Mol:
    """
    Does the reverse of hapticBondsToDative.  If there are multiple
    contiguous atoms attached by dative bonds to an atom (probably a metal
    atom), the dative bonds will be replaced by a dummy atom in their
    centre attached to the (metal) atom by a dative bond, which is
    labelled with ENDPTS of the atoms that had the original dative bonds.
    ARGUMENTS:

    mol: the molecule to use

    RETURNS:a modified copy of the molecule

    C++ signature :RDKit::ROMol* DativeBondsToHaptic(RDKit::ROMol)"""
    ...

def DeleteSubstructs(
    self, mol: Mol, query: Mol, onlyFrags: bool = False, useChirality: bool = False
) -> Mol:
    """
    Removes atoms matching a substructure query from a molecule

    ARGUMENTS:

    mol: the molecule to be modified
    query: the molecule to be used as a substructure query
    onlyFrags: (optional) if this toggle is set, atoms will only be removed if
    the entire fragment in which they are found is matched by the query.
    See below for examples.
    Default value is 0 (remove the atoms whether or not the entire fragment matches)
    useChirality: (optional) match the substructure query using chirality

    RETURNS: a new molecule with the substructure removed
    NOTES:

    The original molecule is not modified.

    EXAMPLES:

    The following examples substitute SMILES/SMARTS strings for molecules, you’d have
    to actually use molecules:

    DeleteSubstructs(‘CCOC’,’OC’) -> ‘CC’
    DeleteSubstructs(‘CCOC’,’OC’,1) -> ‘CCOC’
    DeleteSubstructs(‘CCOCCl.Cl’,’Cl’,1) -> ‘CCOCCl’
    DeleteSubstructs(‘CCOCCl.Cl’,’Cl’) -> ‘CCOC’

    C++ signature :RDKit::ROMol* DeleteSubstructs(RDKit::ROMol,RDKit::ROMol [,bool=False [,bool=False]])
    """
    ...

def DetectBondStereoChemistry(self, mol: Mol, conformer: Conformer) -> None:
    """
    Assign stereochemistry to bonds based on coordinates and a conformer.
    DEPRECATED

    ARGUMENTS:

    mol: the molecule to be modified
    conformer: Conformer providing the coordinates

    C++ signature :void DetectBondStereoChemistry(RDKit::ROMol {lvalue},RDKit::Conformer const*)
    """
    ...

def DetectBondStereochemistry(self, mol: Mol, confId: int = -1) -> None:
    """
    DEPRECATED, use SetDoubleBondNeighborDirections() insteadARGUMENTS:

    mol: the molecule to be modified
    confId: Conformer to use for the coordinates

    C++ signature :void DetectBondStereochemistry(RDKit::ROMol {lvalue} [,int=-1])"""
    ...

def DetectChemistryProblems(
    self, mol: Mol, sanitizeOps: int = SanitizeFlags.SANITIZE_ALL
) -> tuple:
    """
    checks for chemistry problems

    C++ signature :boost::python::tuple DetectChemistryProblems(RDKit::ROMol [,unsigned int=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL])
    """
    ...

def FastFindRings(self, arg1: Mol) -> None:
    """
    Does a non-SSSR ring finding for a molecule.

    ARGUMENTS:

    mol: the molecule to use.

    RETURNS: Nothing

    C++ signature :void FastFindRings(RDKit::ROMol)"""
    ...

def FindAllPathsOfLengthN(
    self,
    mol: Mol,
    length: int,
    useBonds: bool = True,
    useHs: bool = False,
    rootedAtAtom: int = -1,
    onlyShortestPaths: bool = False,
) -> _listSt6vectorIiSaIiEE:
    """
    Finds all paths of a particular length in a molecule

    ARGUMENTS:

    mol: the molecule to use
    length: an integer with the target length for the paths.
    useBonds: (optional) toggles the use of bond indices in the paths.
    Otherwise atom indices are used.  Note this behavior is different
    from that for subgraphs.
    Defaults to 1.
    rootedAtAtom: (optional) if nonzero, only paths from the specified
    atom will be returned.
    onlyShortestPaths: (optional) if set then only paths which are <= the shortest
    path between the begin and end atoms will be included in the results

    RETURNS: a tuple of tuples with IDs for the bonds.
    NOTES:

    Difference between _subgraphs_ and _paths_
    Subgraphs are potentially branched, whereas paths (in our
    terminology at least) cannot be.  So, the following graph:

         C--0--C--1--C--3--C
               |
               2
               |
               C

    has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
    but only 2 _paths_ of length 3: (0,1,3),(2,1,3)

    C++ signature :std::__cxx11::list<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > FindAllPathsOfLengthN(RDKit::ROMol,unsigned int [,bool=True [,bool=False [,int=-1 [,bool=False]]]])
    """
    ...

def FindAllSubgraphsOfLengthMToN(
    self, mol: Mol, min: int, max: int, useHs: bool = False, rootedAtAtom: int = -1
) -> object:
    """
    Finds all subgraphs of a particular length in a moleculeSee documentation for FindAllSubgraphsOfLengthN for definitions

    C++ signature :boost::python::api::object FindAllSubgraphsOfLengthMToN(RDKit::ROMol,unsigned int,unsigned int [,bool=False [,int=-1]])
    """
    ...

def FindAllSubgraphsOfLengthN(
    self, mol: Mol, length: int, useHs: bool = False, rootedAtAtom: int = -1
) -> _listSt6vectorIiSaIiEE:
    """
    Finds all subgraphs of a particular length in a molecule

    ARGUMENTS:

    mol: the molecule to use
    length: an integer with the target number of bonds for the subgraphs.
    useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
    should be included in the results.
    Defaults to 0.
    rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
    atom will be returned.

    RETURNS: a tuple of 2-tuples with bond IDs
    NOTES:

    Difference between _subgraphs_ and _paths_
    Subgraphs are potentially branched, whereas paths (in our
    terminology at least) cannot be.  So, the following graph:

         C--0--C--1--C--3--C
               |
               2
               |
               C

    has 3 _subgraphs_ of length 3: (0,1,2),(0,1,3),(2,1,3)
    but only 2 _paths_ of length 3: (0,1,3),(2,1,3)

    C++ signature :std::__cxx11::list<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > FindAllSubgraphsOfLengthN(RDKit::ROMol,unsigned int [,bool=False [,int=-1]])
    """
    ...

def FindAtomEnvironmentOfRadiusN(
    self,
    mol: Mol,
    radius: int,
    rootedAtAtom: int,
    useHs: bool = False,
    enforceSize: bool = True,
    atomMap: AtomPairsParameters = None,
) -> _vecti:
    """
    Find bonds of a particular radius around an atom.
    Return empty result if there is no bond at the requested radius.

    ARGUMENTS:

    mol: the molecule to use
    radius: an integer with the target radius for the environment.
    rootedAtAtom: the atom to consider
    useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
    should be included in the results.
    Defaults to 0.
    enforceSize (optional) If set to False, all bonds within the requested radius is
    collected. Defaults to 1.
    atomMap: (optional) If provided, it will measure the minimum distance of the atom
    from the rooted atom (start with 0 from the rooted atom). The result is a pair of
    the atom ID and the distance.

    RETURNS: a vector of bond IDs

    C++ signature :std::vector<int, std::allocator<int> > FindAtomEnvironmentOfRadiusN(RDKit::ROMol,unsigned int,unsigned int [,bool=False [,bool=True [,boost::python::api::object=None]]])
    """
    ...

def FindPotentialStereo(
    self, mol: Mol, cleanIt: bool = False, flagPossible: bool = True
) -> _vectN5RDKit9Chirality10StereoInfoE:
    """
    find potential stereo elements in a molecule and returns them as StereoInfo objects
    Note that this function is still somewhat experimental and the API
    and results may change in a future release.

    C++ signature :std::vector<RDKit::Chirality::StereoInfo, std::allocator<RDKit::Chirality::StereoInfo> > FindPotentialStereo(RDKit::ROMol {lvalue} [,bool=False [,bool=True]])
    """
    ...

def FindPotentialStereoBonds(self, mol: Mol, cleanIt: bool = False) -> None:
    """
    Find bonds than can be cis/trans in a molecule and mark them as ‘any’.
    This function finds any double bonds that can potentially be part
    of a cis/trans system. No attempt is made here to mark them cis or trans

    ARGUMENTS:

    mol: the molecule to use

    cleanIt: (optional) if this option is set to true, any previous marking of _CIPCodeon the bond is cleared - otherwise it is left untouched

    C++ signature :void FindPotentialStereoBonds(RDKit::ROMol {lvalue} [,bool=False])"""
    ...

def FindRingFamilies(self, arg1: Mol) -> None:
    """
    generate Unique Ring Families

    C++ signature :void FindRingFamilies(RDKit::ROMol)"""
    ...

def FindUniqueSubgraphsOfLengthN(
    self,
    mol: Mol,
    length: int,
    useHs: bool = False,
    useBO: bool = True,
    rootedAtAtom: int = -1,
) -> _listSt6vectorIiSaIiEE:
    """
    Finds unique subgraphs of a particular length in a molecule

    ARGUMENTS:

    mol: the molecule to use
    length: an integer with the target number of bonds for the subgraphs.
    useHs: (optional) toggles whether or not bonds to Hs that are part of the graph
    should be included in the results.
    Defaults to 0.
    useBO: (optional) Toggles use of bond orders in distinguishing one subgraph from
    another.
    Defaults to 1.
    rootedAtAtom: (optional) if nonzero, only subgraphs from the specified
    atom will be returned.

    RETURNS: a tuple of tuples with bond IDs

    C++ signature :std::__cxx11::list<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > FindUniqueSubgraphsOfLengthN(RDKit::ROMol,unsigned int [,bool=False [,bool=True [,int=-1]]])
    """
    ...

def FragmentOnBRICSBonds(self, mol: Mol) -> Mol:
    """
    Return a new molecule with all BRICS bonds broken

    C++ signature :RDKit::ROMol* FragmentOnBRICSBonds(RDKit::ROMol)"""
    ...

def FragmentOnBonds(
    self,
    mol: Mol,
    bondIndices: AtomPairsParameters,
    addDummies: bool = True,
    dummyLabels: AtomPairsParameters = None,
    bondTypes: AtomPairsParameters = None,
    cutsPerAtom: list = [],
) -> Mol:
    """
    Return a new molecule with all specified bonds broken

    ARGUMENTS:

    mol            - the molecule to be modified
    bondIndices    - indices of the bonds to be broken
    addDummies  - toggles addition of dummy atoms to indicate where
    bonds were broken
    dummyLabels - used to provide the labels to be used for the dummies.
    the first element in each pair is the label for the dummy
    that replaces the bond’s beginAtom, the second is for the
    dummy that replaces the bond’s endAtom. If not provided, the
    dummies are labeled with atom indices.
    bondTypes - used to provide the bond type to use between the
    fragments and the dummy atoms. If not provided, defaults to single.
    cutsPerAtom - used to return the number of cuts made at each atom.

    RETURNS:a new Mol with the modifications

    C++ signature :RDKit::ROMol* FragmentOnBonds(RDKit::ROMol,boost::python::api::object [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None [,boost::python::list=[]]]]])
    """
    ...

def FragmentOnSomeBonds(
    self,
    mol: Mol,
    bondIndices: AtomPairsParameters,
    numToBreak: int = 1,
    addDummies: bool = True,
    dummyLabels: AtomPairsParameters = None,
    bondTypes: AtomPairsParameters = None,
    returnCutsPerAtom: bool = False,
) -> tuple:
    """
    fragment on some bonds

    C++ signature :boost::python::tuple FragmentOnSomeBonds(RDKit::ROMol,boost::python::api::object [,unsigned int=1 [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None [,bool=False]]]]])
    """
    ...

def Get3DDistanceMatrix(
    self,
    mol: Mol,
    confId: int = -1,
    useAtomWts: bool = False,
    force: bool = False,
    prefix: str = "",
) -> object:
    """
    Returns the molecule’s 3D distance matrix.

    ARGUMENTS:

    mol: the molecule to use
    confId: (optional) chooses the conformer Id to use
    Default value is -1.
    useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
    matrix (to return a “Balaban” distance matrix).
    Default value is 0.
    force: (optional) forces the calculation to proceed, even if there is a cached value.
    Default value is 0.
    prefix: (optional, internal use) sets the prefix used in the property cache
    Default value is .

    RETURNS: a Numeric array of floats with the distance matrix

    C++ signature :_object* Get3DDistanceMatrix(RDKit::ROMol {lvalue} [,int=-1 [,bool=False [,bool=False [,char const*=’’]]]])
    """
    ...

def GetAdjacencyMatrix(
    self,
    mol: Mol,
    useBO: bool = False,
    emptyVal: int = 0,
    force: bool = False,
    prefix: str = "",
) -> object:
    """
    Returns the molecule’s adjacency matrix.

    ARGUMENTS:

    mol: the molecule to use
    useBO: (optional) toggles use of bond orders in calculating the matrix.
    Default value is 0.
    emptyVal: (optional) sets the elements of the matrix between non-adjacent atoms
    Default value is 0.
    force: (optional) forces the calculation to proceed, even if there is a cached value.
    Default value is 0.
    prefix: (optional, internal use) sets the prefix used in the property cache
    Default value is .

    RETURNS: a Numeric array of floats containing the adjacency matrix

    C++ signature :_object* GetAdjacencyMatrix(RDKit::ROMol {lvalue} [,bool=False [,int=0 [,bool=False [,char const*=’’]]]])
    """
    ...

def GetAllowNontetrahedralChirality(self) -> bool:
    """
    returns whether or not recognition of non-tetrahedral chirality from 3D structures is enabled

    C++ signature :bool GetAllowNontetrahedralChirality()"""
    ...

def GetDistanceMatrix(
    self,
    mol: Mol,
    useBO: bool = False,
    useAtomWts: bool = False,
    force: bool = False,
    prefix: str = "",
) -> object:
    """
    Returns the molecule’s topological distance matrix.

    ARGUMENTS:

    mol: the molecule to use
    useBO: (optional) toggles use of bond orders in calculating the distance matrix.
    Default value is 0.
    useAtomWts: (optional) toggles using atom weights for the diagonal elements of the
    matrix (to return a “Balaban” distance matrix).
    Default value is 0.
    force: (optional) forces the calculation to proceed, even if there is a cached value.
    Default value is 0.
    prefix: (optional, internal use) sets the prefix used in the property cache
    Default value is .

    RETURNS: a Numeric array of floats with the distance matrix

    C++ signature :_object* GetDistanceMatrix(RDKit::ROMol {lvalue} [,bool=False [,bool=False [,bool=False [,char const*=’’]]]])
    """
    ...

def GetFormalCharge(self, arg1: Mol) -> int:
    """
    Returns the formal charge for the molecule.

    ARGUMENTS:

    mol: the molecule to use

    C++ signature :int GetFormalCharge(RDKit::ROMol)"""
    ...

def GetMolFrags(
    self,
    mol: Mol,
    asMols: bool = False,
    sanitizeFrags: bool = True,
    frags: AtomPairsParameters = None,
    fragsMolAtomMapping: AtomPairsParameters = None,
) -> tuple:
    """
    Finds the disconnected fragments from a molecule.

    For example, for the molecule ‘CC(=O)[O-].[NH3+]C’ GetMolFrags() returns
    ((0, 1, 2, 3), (4, 5))
    ARGUMENTS:

    mol: the molecule to use
    asMols: (optional) if this is provided and true, the fragments
    will be returned as molecules instead of atom ids.
    sanitizeFrags: (optional) if this is provided and true, the fragments
    molecules will be sanitized before returning them.

    frags: (optional, defaults to None) if asMols is true and this is providedas an empty list, the result will be mol.GetNumAtoms() long on return and
    will contain the fragment assignment for each Atom

    fragsMolAtomMapping: (optional, defaults to None) if asMols is true and this
    is provided as an empty list, the result will be numFrags long on
    return, and each entry will contain the indices of the Atoms in that fragment:
    [(0, 1, 2, 3), (4, 5)]

    RETURNS: a tuple of tuples with IDs for the atoms in each fragmentor a tuple of molecules.

    C++ signature :boost::python::tuple GetMolFrags(RDKit::ROMol [,bool=False [,bool=True [,boost::python::api::object=None [,boost::python::api::object=None]]]])
    """
    ...

def GetMostSubstitutedCoreMatch(
    self, mol: Mol, core: Mol, matches: AtomPairsParameters
) -> object:
    """
    Postprocesses the results of a mol.GetSubstructMatches(core) call
    where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups).
    It returns the match with the largest number of non-hydrogen matches to
    the terminal dummy atoms.

    ARGUMENTS:

    mol: the molecule GetSubstructMatches was run on
    core: the molecule used as a substructure query
    matches: the result returned by GetSubstructMatches

    RETURNS: the tuple where terminal dummy atoms in the core match the largest number of non-hydrogen atoms in mol

    C++ signature :_object* GetMostSubstitutedCoreMatch(RDKit::ROMol,RDKit::ROMol,boost::python::api::object)
    """
    ...

def GetSSSR(self, mol: Mol, includeDativeBonds: bool = False) -> _vectSt6vectorIiSaIiEE:
    """
    Get the smallest set of simple rings for a molecule.

    ARGUMENTS:

    mol: the molecule to use.
    includeDativeBonds: whether or not dative bonds should be included in the ring finding.

    RETURNS: a sequence of sequences containing the rings found as atom idsThe length of this will be equal to NumBonds-NumAtoms+1 for single-fragment molecules.

    C++ signature :std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > GetSSSR(RDKit::ROMol {lvalue} [,bool=False])
    """
    ...

def GetShortestPath(self, arg1: Mol, arg2: int, arg3: int) -> tuple:
    """
    Find the shortest path between two atoms using the Bellman-Ford algorithm.

    ARGUMENTS:

    mol: the molecule to use
    idx1: index of the first atom
    idx2: index of the second atom

    C++ signature :boost::python::tuple GetShortestPath(RDKit::ROMol,int,int)"""
    ...

def GetSymmSSSR(
    self, mol: Mol, includeDativeBonds: bool = False
) -> _vectSt6vectorIiSaIiEE:
    """
    Get a symmetrized SSSR for a molecule.

    The symmetrized SSSR is at least as large as the SSSR for a molecule.
    In certain highly-symmetric cases (e.g. cubane), the symmetrized SSSR can be
    a bit larger (i.e. the number of symmetrized rings is >= NumBonds-NumAtoms+1).
    ARGUMENTS:

    mol: the molecule to use.
    includeDativeBonds: whether or not dative bonds should be included in the ring finding.

    RETURNS: a sequence of sequences containing the rings found as atom ids

    C++ signature :std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > GetSymmSSSR(RDKit::ROMol {lvalue} [,bool=False])
    """
    ...

def GetUseLegacyStereoPerception(self) -> bool:
    """
    returns whether or not the legacy stereo perception code is being used

    C++ signature :bool GetUseLegacyStereoPerception()"""
    ...

def HapticBondsToDative(self, mol: Mol) -> Mol:
    """
    One way of showing haptic bonds (such as cyclopentadiene to
    iron in ferrocene) is to use a dummy atom with a dative bond to the
    iron atom with the bond labelled with the atoms involved in the
    organic end of the bond.  Another way is to have explicit dative
    bonds from the atoms of the haptic group to the metal atom.  This
    function converts the former representation to the latter.
    ARGUMENTS:

    mol: the molecule to use

    RETURNS:a modified copy of the molecule

    C++ signature :RDKit::ROMol* HapticBondsToDative(RDKit::ROMol)"""
    ...

def Kekulize(self, mol: Mol, clearAromaticFlags: bool = False) -> None:
    """
    Kekulizes the molecule

    ARGUMENTS:

    mol: the molecule to use
    clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
    molecule will be marked non-aromatic following the kekulization.
    Default value is False.

    NOTES:

    The molecule is modified in place.
    this does not modify query bonds which have bond type queries (like those
    which come from SMARTS) or rings containing them.
    even if clearAromaticFlags is False the BondType for all modified
    aromatic bonds will be changed from AROMATIC to SINGLE or DOUBLE
    Kekulization.

    C++ signature :void Kekulize(RDKit::ROMol {lvalue} [,bool=False])"""
    ...

def KekulizeIfPossible(self, mol: Mol, clearAromaticFlags: bool = False) -> None:
    """
    Kekulizes the molecule if possible. Otherwise the molecule is not modified

    ARGUMENTS:

    mol: the molecule to use
    clearAromaticFlags: (optional) if this toggle is set, all atoms and bonds in the
    molecule will be marked non-aromatic if the kekulization succeds.
    Default value is False.

    NOTES:

    The molecule is modified in place.

    C++ signature :void KekulizeIfPossible(RDKit::ROMol {lvalue} [,bool=False])"""
    ...

def LayeredFingerprint(
    self,
    mol: Mol,
    layerFlags: int = 4294967295,
    minPath: int = 1,
    maxPath: int = 7,
    fpSize: int = 2048,
    atomCounts: list = [],
    setOnlyBits: ExplicitBitVect = None,
    branchedPaths: bool = True,
    fromAtoms: AtomPairsParameters = 0,
) -> ExplicitBitVect:
    """
    Returns a layered fingerprint for a molecule

    NOTE: This function is experimental. The API or results may change fromrelease to release.

    Explanation of the algorithm below.
    ARGUMENTS:

    mol: the molecule to use
    layerFlags: (optional) which layers to include in the fingerprint
    See below for definitions. Defaults to all.
    minPath: (optional) minimum number of bonds to include in the subgraphs
    Defaults to 1.
    maxPath: (optional) maximum number of bonds to include in the subgraphs
    Defaults to 7.
    fpSize: (optional) number of bits in the fingerprint
    Defaults to 2048.
    atomCounts: (optional)
    if provided, this should be a list at least as long as the number of atoms
    in the molecule. It will be used to provide the count of the number
    of paths that set bits each atom is involved in.
    NOTE: the list is not zeroed out here.
    setOnlyBits: (optional)
    if provided, only bits that are set in this bit vector will be set
    in the result. This is essentially the same as doing:

    res &= setOnlyBits

    but also has an impact on the atomCounts (if being used)

    branchedPaths: (optional) if set both branched and unbranched paths will be
    used in the fingerprint.
    Defaults to True.
    fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs
    starting from these atoms will be used.
    Defaults to empty.

    RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits

    Layer definitions:
    0x01: pure topology
    0x02: bond order
    0x04: atom types
    0x08: presence of rings
    0x10: ring sizes
    0x20: aromaticity

    C++ signature :ExplicitBitVect* LayeredFingerprint(RDKit::ROMol [,unsigned int=4294967295 [,unsigned int=1 [,unsigned int=7 [,unsigned int=2048 [,boost::python::list=[] [,ExplicitBitVect*=None [,bool=True [,boost::python::api::object=0]]]]]]]])
    """
    ...

def MergeQueryHs(
    self, mol: Mol, mergeUnmappedOnly: bool = False, mergeIsotopes: bool = False
) -> Mol:
    """
    merges hydrogens into their neighboring atoms as queries

    C++ signature :RDKit::ROMol* MergeQueryHs(RDKit::ROMol [,bool=False [,bool=False]])
    """
    ...

def MolAddRecursiveQueries(self, mol: Mol, queries: dict, propName: str) -> None:
    """
    Adds named recursive queries to atoms

    C++ signature :void MolAddRecursiveQueries(RDKit::ROMol {lvalue},boost::python::dict,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def MurckoDecompose(self, mol: Mol) -> Mol:
    """
    Do a Murcko decomposition and return the scaffold

    C++ signature :RDKit::ROMol* MurckoDecompose(RDKit::ROMol)"""
    ...

def ParseMolQueryDefFile(
    self,
    fileobj: AtomPairsParameters,
    standardize: bool = True,
    delimiter: str = "\t",
    comment: str = "//",
    nameColumn: int = 0,
    smartsColumn: int = 1,
) -> dict:
    """
    reads query definitions from a simply formatted file

    C++ signature :boost::python::dict ParseMolQueryDefFile(boost::python::api::object {lvalue} [,bool=True [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’t’ [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’//’ [,unsigned int=0 [,unsigned int=1]]]]])
    """
    ...

def PathToSubmol(
    self,
    mol: Mol,
    path: AtomPairsParameters,
    useQuery: bool = False,
    atomMap: AtomPairsParameters = None,
) -> Mol:
    """
    C++ signature :RDKit::ROMol* PathToSubmol(RDKit::ROMol,boost::python::api::object {lvalue} [,bool=False [,boost::python::api::object=None]])
    """
    ...

@overload
def PatternFingerprint(
    self,
    mol: Mol,
    fpSize: int = 2048,
    atomCounts: list = [],
    setOnlyBits: ExplicitBitVect = None,
    tautomerFingerprints: bool = False,
) -> ExplicitBitVect:
    """
    A fingerprint using SMARTS patterns

    NOTE: This function is experimental. The API or results may change fromrelease to release.

    C++ signature :ExplicitBitVect* PatternFingerprint(RDKit::MolBundle [,unsigned int=2048 [,ExplicitBitVect*=None [,bool=False]]])
    """
    ...

@overload
def PatternFingerprint(
    self,
    mol: MolBundle,
    fpSize: int = 2048,
    setOnlyBits: ExplicitBitVect = None,
    tautomerFingerprints: bool = False,
) -> ExplicitBitVect: ...
def RDKFingerprint(
    self,
    mol: Mol,
    minPath: int = 1,
    maxPath: int = 7,
    fpSize: int = 2048,
    nBitsPerHash: int = 2,
    useHs: bool = True,
    tgtDensity: float = 0.0,
    minSize: int = 128,
    branchedPaths: bool = True,
    useBondOrder: bool = True,
    atomInvariants: AtomPairsParameters = 0,
    fromAtoms: AtomPairsParameters = 0,
    atomBits: AtomPairsParameters = None,
    bitInfo: AtomPairsParameters = None,
) -> ExplicitBitVect:
    """
    Returns an RDKit topological fingerprint for a molecule

    Explanation of the algorithm below.
    ARGUMENTS:

    mol: the molecule to use
    minPath: (optional) minimum number of bonds to include in the subgraphs
    Defaults to 1.
    maxPath: (optional) maximum number of bonds to include in the subgraphs
    Defaults to 7.
    fpSize: (optional) number of bits in the fingerprint
    Defaults to 2048.
    nBitsPerHash: (optional) number of bits to set per path
    Defaults to 2.
    useHs: (optional) include paths involving Hs in the fingerprint if the molecule
    has explicit Hs.
    Defaults to True.
    tgtDensity: (optional) fold the fingerprint until this minimum density has
    been reached
    Defaults to 0.
    minSize: (optional) the minimum size the fingerprint will be folded to when
    trying to reach tgtDensity
    Defaults to 128.
    branchedPaths: (optional) if set both branched and unbranched paths will be
    used in the fingerprint.
    Defaults to True.
    useBondOrder: (optional) if set both bond orders will be used in the path hashes
    Defaults to True.
    atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
    Defaults to empty.
    fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs
    starting from these atoms will be used.
    Defaults to empty.
    atomBits: (optional) an empty list. If provided, the result will contain a list
    containing the bits each atom sets.
    Defaults to empty.
    bitInfo: (optional) an empty dict. If provided, the result will contain a dict
    with bits as keys and corresponding bond paths as values.
    Defaults to empty.

    RETURNS: a DataStructs.ExplicitBitVect with _fpSize_ bits
    ALGORITHM:

    This algorithm functions by find all subgraphs between minPath and maxPath in
    length.  For each subgraph:

    A hash is calculated.
    The hash is used to seed a random-number generator
    _nBitsPerHash_ random numbers are generated and used to set the corresponding
    bits in the fingerprint

    C++ signature :ExplicitBitVect* RDKFingerprint(RDKit::ROMol [,unsigned int=1 [,unsigned int=7 [,unsigned int=2048 [,unsigned int=2 [,bool=True [,double=0.0 [,unsigned int=128 [,bool=True [,bool=True [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=None [,boost::python::api::object=None]]]]]]]]]]]]])
    """
    ...

def ReapplyMolBlockWedging(self, arg1: Mol) -> None:
    """
    Set the wedging to that which was read from the originalMolBlock, over-riding anything that was originally there.

    ARGUMENTS:

    molecule: the molecule to update

    C++ signature :void ReapplyMolBlockWedging(RDKit::ROMol {lvalue})"""
    ...

def RemoveAllHs(self, mol: Mol, sanitize: bool = True) -> Mol:
    """
    Returns a copy of the molecule with all Hs removed.

    C++ signature :RDKit::ROMol* RemoveAllHs(RDKit::ROMol [,bool=True])"""
    ...

@overload
def RemoveHs(
    self,
    mol: Mol,
    implicitOnly: bool = False,
    updateExplicitCount: bool = False,
    sanitize: bool = True,
) -> Mol:
    """
    Returns a copy of the molecule with Hs removed. Which Hs are removed is controlled by the params argument

    C++ signature :RDKit::ROMol* RemoveHs(RDKit::ROMol,RDKit::MolOps::RemoveHsParameters [,bool=True])
    """
    ...

@overload
def RemoveHs(
    self, mol: Mol, params: RemoveHsParameters, sanitize: bool = True
) -> Mol: ...
def RemoveStereochemistry(self, mol: Mol) -> None:
    """
    Removes all stereochemistry info from the molecule.

    C++ signature :void RemoveStereochemistry(RDKit::ROMol {lvalue})"""
    ...

def RenumberAtoms(self, mol: Mol, newOrder: AtomPairsParameters) -> Mol:
    """
    Returns a copy of a molecule with renumbered atoms

    ARGUMENTS:

    mol: the molecule to be modified
    newOrder: the new ordering the atoms (should be numAtoms long)
    for example: if newOrder is [3,2,0,1], then atom 3 in the original
    molecule will be atom 0 in the new one

    C++ signature :RDKit::ROMol* RenumberAtoms(RDKit::ROMol,boost::python::api::object {lvalue})
    """
    ...

@overload
def ReplaceCore(
    self,
    mol: Mol,
    core: Mol,
    matches: AtomPairsParameters,
    replaceDummies: bool = True,
    labelByIndex: bool = False,
    requireDummyMatch: bool = False,
) -> Mol:
    """
    Removes the core of a molecule and labels the sidechains with dummy atoms.

    ARGUMENTS:

    mol: the molecule to be modified
    coreQuery: the molecule to be used as a substructure query for recognizing the core
    replaceDummies: toggles replacement of atoms that match dummies in the query
    labelByIndex: toggles labeling the attachment point dummy atoms with
    the index of the core atom they’re attached to.
    requireDummyMatch: if the molecule has side chains that attach at points not
    flagged with a dummy, it will be rejected (None is returned)
    useChirality: use chirality matching in the coreQuery

    RETURNS: a new molecule with the core removed
    NOTES:

    The original molecule is not modified.

    EXAMPLES:

    >>> from rdkit.Chem import MolToSmiles, MolFromSmiles, MolFromSmarts, ReplaceCore

    Basic usage: remove a core as specified by SMILES (or another molecule).
    To get the atom labels which are stored as an isotope of the matched atom,
    the output must be written as isomeric smiles.
    A small confusion is that atom isotopes of 0 aren’t shown in smiles strings.
    Here we remove a ring and leave the decoration (r-group) behind.
    >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCCC1CCC1'),MolFromSmiles('C1CCC1')),
    ...             isomericSmiles=True)
    '[1*]CCC'

    The isotope label by default is matched by the first connection found. In order to
    indicate which atom the decoration is attached in the core query, use labelByIndex=True.
    Here the attachment is from the third atom in the smiles string, which is indexed by 3
    in the core, like all good computer scientists expect, atoms indices start at 0.
    >>> MolToSmiles(ReplaceCore(MolFromSmiles('CCN1CCC1'),MolFromSmiles('C1CCN1'),
    ...                         labelByIndex=True),
    ...   isomericSmiles=True)
    '[3*]CC'

    Non-core matches just return None
    >>> ReplaceCore(MolFromSmiles('CCC1CC1'),MolFromSmiles('C1CCC1'))

    The bond between atoms are considered part of the core and are removed as well
    >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CC2C1CCC2'),MolFromSmiles('C1CCC1')),
    ...             isomericSmiles=True)
    '[1*]CCC[2*]'
    >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmiles('N')),
    ...             isomericSmiles=True)
    '[1*]CCCC[2*]'

    When using dummy atoms, cores should be read in as SMARTS.  When read as SMILES
    dummy atoms only match other dummy atoms.
    The replaceDummies flag indicates whether matches to the dummy atoms should be considered as part
    of the core or as part of the decoration (r-group)
    >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),
    ...                         replaceDummies=True),
    ...             isomericSmiles=True)
    '[1*]CC[2*]'
    >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CNCC1'),MolFromSmarts('[*]N[*]'),
    ...                         replaceDummies=False),
    ...             isomericSmiles=True)
    '[1*]CCCC[2*]'

    >>> MolToSmiles(ReplaceCore(MolFromSmiles('C1CCC1CN'),MolFromSmarts('C1CCC1[*]'),
    ...                         replaceDummies=False),
    ...             isomericSmiles=True)
    '[1*]CN'

    C++ signature :RDKit::ROMol* ReplaceCore(RDKit::ROMol,RDKit::ROMol [,bool=True [,bool=False [,bool=False [,bool=False]]]])
    """
    ...

@overload
def ReplaceCore(
    self,
    mol: Mol,
    coreQuery: Mol,
    replaceDummies: bool = True,
    labelByIndex: bool = False,
    requireDummyMatch: bool = False,
    useChirality: bool = False,
) -> Mol: ...
def ReplaceSidechains(
    self, mol: Mol, coreQuery: Mol, useChirality: bool = False
) -> Mol:
    """
    Replaces sidechains in a molecule with dummy atoms for their attachment points.

    ARGUMENTS:

    mol: the molecule to be modified
    coreQuery: the molecule to be used as a substructure query for recognizing the core
    useChirality: (optional) match the substructure query using chirality

    RETURNS: a new molecule with the sidechains removed
    NOTES:

    The original molecule is not modified.

    EXAMPLES:

    The following examples substitute SMILES/SMARTS strings for molecules, you’d have
    to actually use molecules:

    ReplaceSidechains(‘CCC1CCC1’,’C1CCC1’) -> ‘[Xa]C1CCC1’
    ReplaceSidechains(‘CCC1CC1’,’C1CCC1’) -> ‘’
    ReplaceSidechains(‘C1CC2C1CCC2’,’C1CCC1’) -> ‘[Xa]C1CCC1[Xb]’

    C++ signature :RDKit::ROMol* ReplaceSidechains(RDKit::ROMol,RDKit::ROMol [,bool=False])
    """
    ...

def ReplaceSubstructs(
    self,
    mol: Mol,
    query: Mol,
    replacement: Mol,
    replaceAll: bool = False,
    replacementConnectionPoint: int = 0,
    useChirality: bool = False,
) -> object:
    """
    Replaces atoms matching a substructure query in a molecule

    ARGUMENTS:

    mol: the molecule to be modified
    query: the molecule to be used as a substructure query
    replacement: the molecule to be used as the replacement
    replaceAll: (optional) if this toggle is set, all substructures matching
    the query will be replaced in a single result, otherwise each result will
    contain a separate replacement.
    Default value is False (return multiple replacements)
    replacementConnectionPoint: (optional) index of the atom in the replacement that
    the bond should be made to.
    useChirality: (optional) match the substructure query using chirality

    RETURNS: a tuple of new molecules with the substructures replaced removed
    NOTES:

    The original molecule is not modified.
    A bond is only formed to the remaining atoms, if any, that were bonded
    to the first atom in the substructure query. (For finer control over
    substructure replacement, consider using ChemicalReaction.)

    EXAMPLES:

    The following examples substitute SMILES/SMARTS strings for molecules, you’d have
    to actually use molecules:

    ReplaceSubstructs(‘CCOC’,’O[CH3]’,’NC’) -> (‘CCNC’,)
    ReplaceSubstructs(‘COCCOC’,’O[CH3]’,’NC’) -> (‘COCCNC’,’CNCCOC’)
    ReplaceSubstructs(‘COCCOC’,’O[CH3]’,’NC’,True) -> (‘CNCCNC’,)
    ReplaceSubstructs(‘COCCOC’,’O[CH3]’,’CN’,True,1) -> (‘CNCCNC’,)
    ReplaceSubstructs(‘CCOC’,’[CH3]O’,’NC’) -> (‘CC.CN’,)

    C++ signature :_object* ReplaceSubstructs(RDKit::ROMol,RDKit::ROMol,RDKit::ROMol [,bool=False [,unsigned int=0 [,bool=False]]])
    """
    ...

def SanitizeMol(
    self,
    mol: Mol,
    sanitizeOps: int = SanitizeFlags.SANITIZE_ALL,
    catchErrors: bool = False,
) -> SanitizeFlags:
    """
    Kekulize, check valencies, set aromaticity, conjugation and hybridization

    The molecule is modified in place.
    If sanitization fails, an exception will be thrown unless catchErrors is set

    ARGUMENTS:

    mol: the molecule to be modified
    sanitizeOps: (optional) sanitization operations to be carried out
    these should be constructed by or’ing together the
    operations in rdkit.Chem.SanitizeFlags
    catchErrors: (optional) if provided, instead of raising an exception
    when sanitization fails (the default behavior), the
    first operation that failed (as defined in rdkit.Chem.SanitizeFlags)
    is returned. Zero is returned on success.

    C++ signature :RDKit::MolOps::SanitizeFlags SanitizeMol(RDKit::ROMol {lvalue} [,unsigned long=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL [,bool=False]])
    """
    ...

def SetAllowNontetrahedralChirality(self, arg1: bool) -> None:
    """
    toggles recognition of non-tetrahedral chirality from 3D structures

    C++ signature :void SetAllowNontetrahedralChirality(bool)"""
    ...

def SetAromaticity(
    self, mol: Mol, model: AromaticityModel = AromaticityModel.AROMATICITY_DEFAULT
) -> None:
    """
    does aromaticity perception

    ARGUMENTS:

    mol: the molecule to use
    model: the model to use

    NOTES:

    The molecule is modified in place.

    C++ signature :void SetAromaticity(RDKit::ROMol {lvalue} [,RDKit::MolOps::AromaticityModel=rdkit.Chem.rdmolops.AromaticityModel.AROMATICITY_DEFAULT])
    """
    ...

def SetBondStereoFromDirections(self, mol: Mol) -> None:
    """
    Uses the directions of neighboring bonds to set cis/trans stereo on double bonds.

    ARGUMENTS:

    mol: the molecule to be modified

    C++ signature :void SetBondStereoFromDirections(RDKit::ROMol {lvalue})"""
    ...

def SetConjugation(self, mol: Mol) -> None:
    """
    finds conjugated bonds

    ARGUMENTS:

    mol: the molecule to use

    NOTES:

    The molecule is modified in place.

    C++ signature :void SetConjugation(RDKit::ROMol {lvalue})"""
    ...

def SetDoubleBondNeighborDirections(
    self, mol: Mol, conf: AtomPairsParameters = None
) -> None:
    """
    Uses the stereo info on double bonds to set the directions of neighboring single bonds

    ARGUMENTS:

    mol: the molecule to be modified

    C++ signature :void SetDoubleBondNeighborDirections(RDKit::ROMol {lvalue} [,boost::python::api::object=None])
    """
    ...

def SetGenericQueriesFromProperties(
    self, mol: Mol, useAtomLabels: bool = True, useSGroups: bool = True
) -> None:
    """
    documentation

    C++ signature :void SetGenericQueriesFromProperties(RDKit::ROMol {lvalue} [,bool=True [,bool=True]])
    """
    ...

def SetHybridization(self, mol: Mol) -> None:
    """
    Assigns hybridization states to atoms

    ARGUMENTS:

    mol: the molecule to use

    NOTES:

    The molecule is modified in place.

    C++ signature :void SetHybridization(RDKit::ROMol {lvalue})"""
    ...

def SetTerminalAtomCoords(self, arg1: Mol, arg2: int, arg3: int) -> None:
    """
    Sets Cartesian coordinates for a terminal atom.

    Useful for growing an atom off a molecule with sensible
    coordinates based on the geometry of the neighbor.
    NOTE: this sets the appropriate coordinates in all of the molecule’s conformers
    ARGUMENTS:

    mol: the molecule the atoms belong to.
    idx: index of the terminal atom whose coordinates are set.
    mol: index of the bonded neighbor atom.

    RETURNS: Nothing

    C++ signature :void SetTerminalAtomCoords(RDKit::ROMol {lvalue},unsigned int,unsigned int)
    """
    ...

def SetUseLegacyStereoPerception(self, arg1: bool) -> None:
    """
    toggles usage of the legacy stereo perception code

    C++ signature :void SetUseLegacyStereoPerception(bool)"""
    ...

def SortMatchesByDegreeOfCoreSubstitution(
    self, mol: Mol, core: Mol, matches: AtomPairsParameters
) -> object:
    """
    Postprocesses the results of a mol.GetSubstructMatches(core) call
    where mol has explicit Hs and core bears terminal dummy atoms (i.e., R groups).
    It returns a copy of matches sorted by decreasing number of non-hydrogen matches
    to the terminal dummy atoms.

    ARGUMENTS:

    mol: the molecule GetSubstructMatches was run on
    core: the molecule used as a substructure query
    matches: the result returned by GetSubstructMatches

    RETURNS: a copy of matches sorted by decreasing number of non-hydrogen matches to the terminal dummy atoms

    C++ signature :_object* SortMatchesByDegreeOfCoreSubstitution(RDKit::ROMol,RDKit::ROMol,boost::python::api::object)
    """
    ...

def SplitMolByPDBChainId(
    self, mol: Mol, whiteList: AtomPairsParameters = None, negateList: bool = False
) -> dict:
    """
    Splits a molecule into pieces based on PDB chain information.

    ARGUMENTS:

    mol: the molecule to use
    whiteList: only residues in this list will be returned
    negateList: if set, negates the white list inclusion logic

    RETURNS: a dictionary keyed by chain id with molecules as the values

    C++ signature :boost::python::dict SplitMolByPDBChainId(RDKit::ROMol [,boost::python::api::object=None [,bool=False]])
    """
    ...

def SplitMolByPDBResidues(
    self, mol: Mol, whiteList: AtomPairsParameters = None, negateList: bool = False
) -> dict:
    """
    Splits a molecule into pieces based on PDB residue information.

    ARGUMENTS:

    mol: the molecule to use
    whiteList: only residues in this list will be returned
    negateList: if set, negates the white list inclusion logic

    RETURNS: a dictionary keyed by residue name with molecules as the values

    C++ signature :boost::python::dict SplitMolByPDBResidues(RDKit::ROMol [,boost::python::api::object=None [,bool=False]])
    """
    ...

def TranslateChiralFlagToStereoGroups(
    self, mol: Mol, zeroFlagGroupType: StereoGroupType = StereoGroupType.STEREO_AND
) -> None:
    """
    Generate enhanced stereo groups based on the status of the chiral flag property.

    Arguments:
    mol: molecule to be modified

    zeroFlagGroupType: how to handle non-grouped stereo centers when thechiral flag is set to zero

    If the chiral flag is set to a value of 1 then all specified tetrahedral
    chiral centers which are not already in StereoGroups will be added to an
    ABS StereoGroup.
    If the chiral flag is set to a value of 0 then all specified tetrahedral
    chiral centers will be added to a StereoGroup of the type zeroFlagGroupType
    If there is no chiral flag set (i.e. the property is not present), the
    molecule will not be modified.

    C++ signature :void TranslateChiralFlagToStereoGroups(RDKit::ROMol {lvalue} [,RDKit::StereoGroupType=rdkit.Chem.rdchem.StereoGroupType.STEREO_AND])
    """
    ...

def UnfoldedRDKFingerprintCountBased(
    self,
    mol: Mol,
    minPath: int = 1,
    maxPath: int = 7,
    useHs: bool = True,
    branchedPaths: bool = True,
    useBondOrder: bool = True,
    atomInvariants: AtomPairsParameters = 0,
    fromAtoms: AtomPairsParameters = 0,
    atomBits: AtomPairsParameters = None,
    bitInfo: AtomPairsParameters = None,
) -> ULongSparseIntVect:
    """
    Returns an unfolded count-based version of the RDKit fingerprint for a molecule
    ARGUMENTS:

    mol: the molecule to use
    minPath: (optional) minimum number of bonds to include in the subgraphs
    Defaults to 1.
    maxPath: (optional) maximum number of bonds to include in the subgraphs
    Defaults to 7.
    useHs: (optional) include paths involving Hs in the fingerprint if the molecule
    has explicit Hs.
    Defaults to True.
    branchedPaths: (optional) if set both branched and unbranched paths will be
    used in the fingerprint.
    Defaults to True.
    useBondOrder: (optional) if set both bond orders will be used in the path hashes
    Defaults to True.
    atomInvariants: (optional) a sequence of atom invariants to use in the path hashes
    Defaults to empty.
    fromAtoms: (optional) a sequence of atom indices. If provided, only paths/subgraphs
    starting from these atoms will be used.
    Defaults to empty.
    atomBits: (optional) an empty list. If provided, the result will contain a list
    containing the bits each atom sets.
    Defaults to empty.
    bitInfo: (optional) an empty dict. If provided, the result will contain a dict
    with bits as keys and corresponding bond paths as values.
    Defaults to empty.

    C++ signature :RDKit::SparseIntVect<unsigned long>* UnfoldedRDKFingerprintCountBased(RDKit::ROMol [,unsigned int=1 [,unsigned int=7 [,bool=True [,bool=True [,bool=True [,boost::python::api::object=0 [,boost::python::api::object=0 [,boost::python::api::object=None [,boost::python::api::object=None]]]]]]]]])
    """
    ...

def WedgeBond(self, arg1: Bond, arg2: int, arg3: Conformer) -> None:
    """
    Set the wedging on an individual bond from a molecule.
    The wedging scheme used is that from Mol files.

    ARGUMENTS:
    bond: the bond to update
    atom ID: the atom from which to do the wedging
    conformer: the conformer to use to determine wedge direction

    C++ signature :void WedgeBond(RDKit::Bond*,unsigned int,RDKit::Conformer const*)"""
    ...

def WedgeMolBonds(
    self, mol: Mol, conformer: Conformer, params: BondWedgingParameters = None
) -> None:
    """
    Set the wedging on single bonds in a molecule.
    The wedging scheme used is that from Mol files.

    ARGUMENTS:

    molecule: the molecule to update
    conformer: the conformer to use to determine wedge direction

    C++ signature :void WedgeMolBonds(RDKit::ROMol {lvalue},RDKit::Conformer const* [,RDKit::Chirality::BondWedgingParameters const*=None])
    """
    ...

@overload
def molzip(self, a: Mol, b: Mol, params: MolzipParams = ...) -> Mol:
    """
    zip together two molecules using the given matching parameters

    C++ signature :RDKit::ROMol* molzip(RDKit::ROMol [,RDKit::MolzipParams=<rdkit.Chem.rdmolops.MolzipParams object at 0x7fd481d097d0>])
    """
    ...

@overload
def molzip(self, a: Mol, params: MolzipParams = ...) -> Mol: ...
def molzipFragments(self, mols: AtomPairsParameters, params: MolzipParams = ...) -> Mol:
    """
    zip together multiple molecules from an R group decomposition
    using the given matching parameters.  The first molecule in the list
    must be the core

    C++ signature :RDKit::ROMol* molzipFragments(boost::python::api::object {lvalue} [,RDKit::MolzipParams=<rdkit.Chem.rdmolops.MolzipParams object at 0x7fd481d09880>])
    """
    ...
