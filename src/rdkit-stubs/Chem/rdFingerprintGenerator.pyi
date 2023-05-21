"""
rdkit.Chem.rdFingerprintGenerator module¶
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.DataStructs.cDataStructs import (
    ExplicitBitVect,
    SparseBitVect,
    UIntSparseIntVect,
    ULongSparseIntVect,
)

AtomPairFP: FPType
MorganFP: FPType
RDKitFP: FPType
TopologicalTorsionFP: FPType

class AdditionalOutput(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AllocateAtomCounts(self, arg1: AdditionalOutput) -> None:
        """
        synonym for CollectAtomCounts()

        C++ signature :void AllocateAtomCounts(RDKit::AdditionalOutput {lvalue})"""
        ...
    def AllocateAtomToBits(self, arg1: AdditionalOutput) -> None:
        """
        synonym for CollectAtomToBits()

        C++ signature :void AllocateAtomToBits(RDKit::AdditionalOutput {lvalue})"""
        ...
    def AllocateBitInfoMap(self, arg1: AdditionalOutput) -> None:
        """
        synonym for CollectBitInfoMap()

        C++ signature :void AllocateBitInfoMap(RDKit::AdditionalOutput {lvalue})"""
        ...
    def AllocateBitPaths(self, arg1: AdditionalOutput) -> None:
        """
        synonym for CollectBitPaths()

        C++ signature :void AllocateBitPaths(RDKit::AdditionalOutput {lvalue})"""
        ...
    def CollectAtomCounts(self, arg1: AdditionalOutput) -> None:
        """
        toggles collection of information about the number of bits each atom is involved in

        C++ signature :void CollectAtomCounts(RDKit::AdditionalOutput {lvalue})"""
        ...
    def CollectAtomToBits(self, arg1: AdditionalOutput) -> None:
        """
        toggle collection of information mapping each atom to the bits it is involved in.

        C++ signature :void CollectAtomToBits(RDKit::AdditionalOutput {lvalue})"""
        ...
    def CollectBitInfoMap(self, arg1: AdditionalOutput) -> None:
        """
        toggles collection of information mapping each atom to more detail about the atom environment (not available from all fingerprints)

        C++ signature :void CollectBitInfoMap(RDKit::AdditionalOutput {lvalue})"""
        ...
    def CollectBitPaths(self, arg1: AdditionalOutput) -> None:
        """
        toggles collection of information matching each atom to information about the paths it is involved in (not available from all fingerprints).

        C++ signature :void CollectBitPaths(RDKit::AdditionalOutput {lvalue})"""
        ...
    def GetAtomCounts(self, arg1: AdditionalOutput) -> object:
        """
        C++ signature :boost::python::api::object GetAtomCounts(RDKit::AdditionalOutput)
        """
        ...
    def GetAtomToBits(self, arg1: AdditionalOutput) -> object:
        """
        C++ signature :boost::python::api::object GetAtomToBits(RDKit::AdditionalOutput)
        """
        ...
    def GetBitInfoMap(self, arg1: AdditionalOutput) -> object:
        """
        C++ signature :boost::python::api::object GetBitInfoMap(RDKit::AdditionalOutput)
        """
        ...
    def GetBitPaths(self, arg1: AdditionalOutput) -> object:
        """
        C++ signature :boost::python::api::object GetBitPaths(RDKit::AdditionalOutput)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AtomInvariantsGenerator(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class AtomPairFingerprintOptions(FingerprintOptions):
    """
    Raises an exception
    This class cannot be instantiated from Python

    property maxDistance¶
    maximum distance to be included

    property minDistance¶
    minimum distance to be included

    property use2D¶
    use 2D distances"""

    ...
    maxDistance: Any
    minDistance: Any
    use2D: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class BondInvariantsGenerator(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FPType(Boost.Python.enum):
    AtomPairFP: FPType = ...
    MorganFP: FPType = ...
    RDKitFP: FPType = ...
    TopologicalTorsionFP: FPType = ...
    names: dict[str, FPType] = ...
    values: dict[int, FPType] = ...
    __slots__: ClassVar[tuple] = ...

class FingeprintGenerator32(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetCountFingerprint(
        self,
        arg1: FingeprintGenerator32,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> UIntSparseIntVect:
        """
        Generates a count fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint

        C++ signature :RDKit::SparseIntVect<unsigned int>* GetCountFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetCountFingerprintAsNumPy(
        self,
        arg1: FingeprintGenerator32,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> object:
        """
        Generates a count fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint

        C++ signature :boost::python::api::object GetCountFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetFingerprint(
        self,
        arg1: FingeprintGenerator32,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> ExplicitBitVect:
        """
        Generates a fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a ExplicitBitVect containing fingerprint

        C++ signature :ExplicitBitVect* GetFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetFingerprintAsNumPy(
        self,
        arg1: FingeprintGenerator32,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> object:
        """
        Generates a fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint

        C++ signature :boost::python::api::object GetFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetInfoString(self, arg1: FingeprintGenerator32) -> str:
        """
        Returns a string containing information about the fingerprint generator

        RETURNS: an information string

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInfoString(RDKit::FingerprintGenerator<unsigned int> const*)
        """
        ...
    def GetOptions(self, arg1: FingeprintGenerator32) -> FingerprintOptions:
        """
        return the fingerprint options object

        C++ signature :RDKit::FingerprintArguments* GetOptions(RDKit::FingerprintGenerator<unsigned int>*)
        """
        ...
    def GetSparseCountFingerprint(
        self,
        arg1: FingeprintGenerator32,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> UIntSparseIntVect:
        """
        Generates a sparse count fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint

        C++ signature :RDKit::SparseIntVect<unsigned int>* GetSparseCountFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetSparseFingerprint(
        self,
        arg1: FingeprintGenerator32,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> SparseBitVect:
        """
        Generates a sparse fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseBitVect containing fingerprint

        C++ signature :SparseBitVect* GetSparseFingerprint(RDKit::FingerprintGenerator<unsigned int> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FingeprintGenerator64(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetCountFingerprint(
        self,
        arg1: FingeprintGenerator64,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> UIntSparseIntVect:
        """
        Generates a count fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint

        C++ signature :RDKit::SparseIntVect<unsigned int>* GetCountFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetCountFingerprintAsNumPy(
        self,
        arg1: FingeprintGenerator64,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> object:
        """
        Generates a count fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint

        C++ signature :boost::python::api::object GetCountFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetFingerprint(
        self,
        arg1: FingeprintGenerator64,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> ExplicitBitVect:
        """
        Generates a fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a ExplicitBitVect containing fingerprint

        C++ signature :ExplicitBitVect* GetFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetFingerprintAsNumPy(
        self,
        arg1: FingeprintGenerator64,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> object:
        """
        Generates a fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a numpy array containing the fingerprint

        C++ signature :boost::python::api::object GetFingerprintAsNumPy(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetInfoString(self, arg1: FingeprintGenerator64) -> str:
        """
        Returns a string containing information about the fingerprint generator

        RETURNS: an information string

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetInfoString(RDKit::FingerprintGenerator<unsigned long> const*)
        """
        ...
    def GetOptions(self, arg1: FingeprintGenerator64) -> FingerprintOptions:
        """
        return the fingerprint options object

        C++ signature :RDKit::FingerprintArguments* GetOptions(RDKit::FingerprintGenerator<unsigned long>*)
        """
        ...
    def GetSparseCountFingerprint(
        self,
        arg1: FingeprintGenerator64,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> ULongSparseIntVect:
        """
        Generates a sparse count fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseIntVect containing fingerprint

        C++ signature :RDKit::SparseIntVect<unsigned long>* GetSparseCountFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    def GetSparseFingerprint(
        self,
        arg1: FingeprintGenerator64,
        mol: Mol,
        fromAtoms: AtomPairsParameters = [],
        ignoreAtoms: AtomPairsParameters = [],
        confId: int = -1,
        customAtomInvariants: AtomPairsParameters = [],
        customBondInvariants: AtomPairsParameters = [],
        additionalOutput: AtomPairsParameters = None,
    ) -> SparseBitVect:
        """
        Generates a sparse fingerprint

        ARGUMENTS:
        mol: molecule to be fingerprinted
        fromAtoms: indices of atoms to use while generating the fingerprint
        ignoreAtoms: indices of atoms to exclude while generating the fingerprint
        confId: 3D confirmation to use, only used by AtomPair fingerprint
        customAtomInvariants: custom atom invariants to be used, overrides invariants from the invariant generator
        customBondInvariants: custom bond invariants to be used, overrides invariants from the invariant generator
        additionalOutput: AdditionalOutput instance used to return extra information about the bits

        RETURNS: a SparseBitVect containing fingerprint

        C++ signature :SparseBitVect* GetSparseFingerprint(RDKit::FingerprintGenerator<unsigned long> const*,RDKit::ROMol [,boost::python::api::object=[] [,boost::python::api::object=[] [,int=-1 [,boost::python::api::object=[] [,boost::python::api::object=[] [,boost::python::api::object=None]]]]]])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FingerprintOptions(Boost.Python.instance):
    """
    Raises an exception
    This class cannot be instantiated from Python"""

    countSimulation: Any
    fpSize: Any
    includeChirality: Any
    numBitsPerFeature: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def SetCountBounds(
        self, arg1: FingerprintOptions, arg2: AtomPairsParameters
    ) -> None:
        """
        set the bins for the count bounds

        C++ signature :void SetCountBounds(RDKit::FingerprintArguments {lvalue},boost::python::api::object)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MorganFingerprintOptions(FingerprintOptions):
    """
    Raises an exception
    This class cannot be instantiated from Python

    property includeRedundantEnvironments¶
    include redundant environments in the fingerprint

    property onlyNonzeroInvariants¶
    use include atoms which have nonzero invariants

    property radius¶
    the radius of the fingerprints to generate"""

    ...
    includeRedundantEnvironments: Any
    onlyNonzeroInvariants: Any
    radius: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RDKitFingerprintOptions(FingerprintOptions):
    """
    Raises an exception
    This class cannot be instantiated from Python

    property branchedPaths¶
    generate branched subgraphs, not just linear ones

    property maxPath¶
    maximum path length (in bonds) to be included

    property minPath¶
    minimum path length (in bonds) to be included

    property useBondOrder¶
    include bond orders in the path hashes

    property useHs¶
    use explicit Hs in the paths (if molecule has explicit Hs)"""

    ...
    branchedPaths: Any
    maxPath: Any
    minPath: Any
    useBondOrder: Any
    useHs: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class TopologicalTorsionFingerprintOptions(FingerprintOptions):
    """
    Raises an exception
    This class cannot be instantiated from Python

    property onlyShortestPaths¶
    whether or not to only include paths which are the shortest path between the start and end atoms

    property torsionAtomCount¶
    number of atoms to be included in the paths"""

    ...
    onlyShortestPaths: Any
    torsionAtomCount: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def GetAtomPairAtomInvGen(
    self, includeChirality: bool = False
) -> AtomInvariantsGenerator:
    """
    Get an atom pair atom-invariant generator

    ARGUMENTS:
    includeChirality: if set, chirality will be taken into account for invariants

    RETURNS: AtomInvariantsGenerator

    C++ signature :RDKit::AtomInvariantsGenerator* GetAtomPairAtomInvGen([ bool=False])
    """
    ...

def GetAtomPairGenerator(
    self,
    minDistance: int = 1,
    maxDistance: int = 30,
    includeChirality: bool = False,
    use2D: bool = True,
    countSimulation: bool = True,
    countBounds: AtomPairsParameters = None,
    fpSize: int = 2048,
    atomInvariantsGenerator: AtomPairsParameters = None,
) -> FingeprintGenerator64:
    """
    Get an atom pair fingerprint generator

    ARGUMENTS:
    minDistance: minimum distance between atoms to be considered in a pair, default is 1 bond
    maxDistance: maximum distance between atoms to be considered in a pair, default is maxPathLen-1 bonds
    includeChirality: if set, chirality will be used in the atom  invariants, this is ignored if atomInvariantsGenerator is provided
    use2D: if set, the 2D (topological) distance matrix  will be used
    countSimulation:  if set, use count simulation while  generating the fingerprint
    countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
    fpSize: size of the generated fingerprint, does not affect the sparse versions
    atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:

    atomToBits: which bits each atom is involved in
    atomCounts: how many bits each atom sets
    bitInfoMap: map from bitId to (atomId, radius) pairs

    RETURNS: FingerprintGenerator

    C++ signature :RDKit::FingerprintGenerator<unsigned long>* GetAtomPairGenerator([ unsigned int=1 [,unsigned int=30 [,bool=False [,bool=True [,bool=True [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,boost::python::api::object {lvalue}=None]]]]]]]])
    """
    ...

def GetCountFPs(self, molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    C++ signature :boost::python::list GetCountFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
    ...

def GetFPs(self, molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    C++ signature :boost::python::list GetFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
    ...

def GetMorganAtomInvGen(
    self, includeRingMembership: bool = False
) -> AtomInvariantsGenerator:
    """
    Get a morgan atom invariants generator

    ARGUMENTS:
    includeRingMembership: if set, whether or not the atom is in a ring will be used in the invariant list

    RETURNS: AtomInvariantsGenerator

    C++ signature :RDKit::AtomInvariantsGenerator* GetMorganAtomInvGen([ bool=False])"""
    ...

def GetMorganBondInvGen(
    self, useBondTypes: bool = True, useChirality: bool = False
) -> BondInvariantsGenerator:
    """
    Get a morgan bond invariants generator

    ARGUMENTS:
    useBondTypes: if set, bond types will be included as a part of the bond invariants
    useChirality: if set, chirality information will be included as a part of the bond invariants

    RETURNS: BondInvariantsGenerator

    C++ signature :RDKit::BondInvariantsGenerator* GetMorganBondInvGen([ bool=True [,bool=False]])
    """
    ...

def GetMorganFeatureAtomInvGen(
    self, patterns: AtomPairsParameters = None
) -> AtomInvariantsGenerator:
    """
    Get a morgan feature atom invariants generator

    ARGUMENTS:
    patterns: if provided should contain the queries used to assign atom-types. if not provided, feature definitions adapted from reference: Gobbi and Poppinger, Biotech. Bioeng. _61_ 47-54 (1998) will be used for Donor, Acceptor, Aromatic, Halogen, Basic, Acidic.

    RETURNS: AtomInvariantsGenerator

    C++ signature :RDKit::AtomInvariantsGenerator* GetMorganFeatureAtomInvGen([ boost::python::api::object {lvalue}=None])
    """
    ...

def GetMorganGenerator(
    self,
    radius: int = 3,
    countSimulation: bool = False,
    includeChirality: bool = False,
    useBondTypes: bool = True,
    onlyNonzeroInvariants: bool = False,
    includeRingMembership: bool = True,
    countBounds: AtomPairsParameters = None,
    fpSize: int = 2048,
    atomInvariantsGenerator: AtomPairsParameters = None,
    bondInvariantsGenerator: AtomPairsParameters = None,
    includeRedundantEnvironments: bool = False,
) -> FingeprintGenerator64:
    """
    Get a morgan fingerprint generator

    ARGUMENTS:
    radius:  the number of iterations to grow the fingerprint
    countSimulation: if set, use count simulation while generating the fingerprint
    includeChirality: if set, chirality information will be added to the generated fingerprint
    useBondTypes: if set, bond types will be included as a part of the default bond invariants
    countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
    fpSize: size of the generated fingerprint, does not affect the sparse versions
    atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:

    atomToBits: which bits each atom is the center of
    atomCounts: how many bits each atom sets
    bitInfoMap: map from bitId to (atomId1, radius) pairs

    RETURNS: FingerprintGenerator

    C++ signature :RDKit::FingerprintGenerator<unsigned long>* GetMorganGenerator([ unsigned int=3 [,bool=False [,bool=False [,bool=True [,bool=False [,bool=True [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,boost::python::api::object {lvalue}=None [,boost::python::api::object {lvalue}=None [,bool=False]]]]]]]]]]])
    """
    ...

def GetRDKitAtomInvGen(self) -> AtomInvariantsGenerator:
    """
    Get an RDKit atom invariants generator

    RETURNS: AtomInvariantsGenerator

    C++ signature :RDKit::AtomInvariantsGenerator* GetRDKitAtomInvGen()"""
    ...

def GetRDKitFPGenerator(
    self,
    minPath: int = 1,
    maxPath: int = 7,
    useHs: bool = True,
    branchedPaths: bool = True,
    useBondOrder: bool = True,
    countSimulation: bool = False,
    countBounds: AtomPairsParameters = None,
    fpSize: int = 2048,
    numBitsPerFeature: int = 2,
    atomInvariantsGenerator: AtomPairsParameters = None,
) -> FingeprintGenerator64:
    """
    Get an RDKit fingerprint generator

    ARGUMENTS:
    minPath: the minimum path length (in bonds) to be included
    maxPath: the maximum path length (in bonds) to be included
    useHs: toggles inclusion of Hs in paths (if the molecule has explicit Hs)
    branchedPaths: toggles generation of branched subgraphs, not just linear paths
    useBondOrder: toggles inclusion of bond orders in the path hashes
    countSimulation:  if set, use count simulation while  generating the fingerprint
    countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
    fpSize: size of the generated fingerprint, does not affect the sparse versions
    numBitsPerFeature: the number of bits set per path/subgraph found
    atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:

    atomToBits: which bits each atom is involved in
    atomCounts: how many bits each atom sets
    bitPaths: map from bitId to vectors of bond indices for the individual subgraphs

    RETURNS: FingerprintGenerator

    C++ signature :RDKit::FingerprintGenerator<unsigned long>* GetRDKitFPGenerator([ unsigned int=1 [,unsigned int=7 [,bool=True [,bool=True [,bool=True [,bool=False [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,unsigned int=2 [,boost::python::api::object {lvalue}=None]]]]]]]]]])
    """
    ...

def GetSparseCountFPs(
    self, molecules: list = [], fpType: FPType = FPType.MorganFP
) -> list:
    """
    C++ signature :boost::python::list GetSparseCountFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
    ...

def GetSparseFPs(self, molecules: list = [], fpType: FPType = FPType.MorganFP) -> list:
    """
    C++ signature :boost::python::list GetSparseFPs([ boost::python::list {lvalue}=[] [,RDKit::FPType=rdkit.Chem.rdFingerprintGenerator.FPType.MorganFP]])
    """
    ...

def GetTopologicalTorsionGenerator(
    self,
    includeChirality: bool = False,
    torsionAtomCount: int = 4,
    countSimulation: bool = True,
    countBounds: AtomPairsParameters = None,
    fpSize: int = 2048,
    atomInvariantsGenerator: AtomPairsParameters = None,
) -> FingeprintGenerator64:
    """
    Get an atom pair fingerprint generator

    ARGUMENTS:
    includeChirality: includeChirality argument for both the default atom invariants generator and the fingerprint arguments
    torsionAtomCount: the number of atoms to include in the “torsions”
    countSimulation:  if set, use count simulation while  generating the fingerprint
    countBounds: boundaries for count simulation, corresponding bit will be  set if the count is higher than the number provided for that spot
    fpSize: size of the generated fingerprint, does not affect the sparse versions
    atomInvariantsGenerator: atom invariants to be used during fingerprint generation

    This generator supports the following AdditionalOutput types:

    atomToBits: which bits each atom is involved in
    atomCounts: how many bits each atom sets
    bitPaths: map from bitId to vectors of atom indices

    RETURNS: FingerprintGenerator

    C++ signature :RDKit::FingerprintGenerator<unsigned long>* GetTopologicalTorsionGenerator([ bool=False [,unsigned int=4 [,bool=True [,boost::python::api::object {lvalue}=None [,unsigned int=2048 [,boost::python::api::object {lvalue}=None]]]]]])
    """
    ...
