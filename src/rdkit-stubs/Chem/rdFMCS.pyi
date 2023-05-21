"""
rdkit.Chem.rdFMCS module¶
Module containing a C++ implementation of the FMCS algorithm
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class AtomCompare(Boost.Python.enum):
    CompareAny: AtomCompare = ...
    CompareAnyHeavyAtom: AtomCompare = ...
    CompareElements: AtomCompare = ...
    CompareIsotopes: AtomCompare = ...
    names: dict[str, AtomCompare] = ...
    values: dict[int, AtomCompare] = ...
    __slots__: ClassVar[tuple] = ...

class BondCompare(Boost.Python.enum):
    CompareAny: BondCompare = ...
    CompareOrder: BondCompare = ...
    CompareOrderExact: BondCompare = ...
    names: dict[str, BondCompare] = ...
    values: dict[int, BondCompare] = ...
    __slots__: ClassVar[tuple] = ...

class MCSAtomCompare(Boost.Python.instance):
    """
    Base class. Subclass and override MCSAtomCompare.__call__() to define custom atom compare functions, then set MCSParameters.AtomTyper to an instance of the subclass

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def CheckAtomCharge(
        self: MCSAtomCompare,
        parameters: MCSAtomCompareParameters,
        mol1: Mol,
        atom1: int,
        mol2: Mol,
        atom2: int,
    ) -> bool:
        """
        Return True if both atoms have the same formal charge

        C++ signature :bool CheckAtomCharge(RDKit::PyMCSAtomCompare {lvalue},RDKit::MCSAtomCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    def CheckAtomChirality(
        self: MCSAtomCompare,
        parameters: MCSAtomCompareParameters,
        mol1: Mol,
        atom1: int,
        mol2: Mol,
        atom2: int,
    ) -> bool:
        """
        Return True if both atoms have, or have not, a chiral tag

        C++ signature :bool CheckAtomChirality(RDKit::PyMCSAtomCompare {lvalue},RDKit::MCSAtomCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    def CheckAtomRingMatch(
        self: MCSAtomCompare,
        parameters: MCSAtomCompareParameters,
        mol1: Mol,
        atom1: int,
        mol2: Mol,
        atom2: int,
    ) -> bool:
        """
        Return True if both atoms are, or are not, in a ring

        C++ signature :bool CheckAtomRingMatch(RDKit::PyMCSAtomCompare {lvalue},RDKit::MCSAtomCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    def compare(
        self: MCSAtomCompare,
        parameters: MCSAtomCompareParameters,
        mol1: Mol,
        atom1: int,
        mol2: Mol,
        atom2: int,
    ) -> bool:
        """
        override to implement custom atom comparison

        C++ signature :bool compare(RDKit::PyMCSAtomCompare {lvalue},RDKit::MCSAtomCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    @classmethod
    def __call__(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MCSAtomCompareParameters(Boost.Python.instance):
    """
    Parameters controlling how atom-atom matching is done

    C++ signature :void __init__(_object*)

    property CompleteRingsOnly¶
    results cannot include lone ring atoms

    property MatchChiralTag¶
    include atom chirality in the match

    property MatchFormalCharge¶
    include formal charge in the match

    property MatchIsotope¶
    use isotope atom queries in MCSResults

    property MatchValences¶
    include atom valences in the match

    property MaxDistance¶
    Require atoms to be within this many angstroms in 3D

    property RingMatchesRingOnly¶
    ring atoms are only allowed to match other ring atoms"""

    ...
    __instance_size__: ClassVar[int] = ...
    CompleteRingsOnly: Any
    MatchChiralTag: Any
    MatchFormalCharge: Any
    MatchIsotope: Any
    MatchValences: Any
    MaxDistance: Any
    RingMatchesRingOnly: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MCSBondCompare(Boost.Python.instance):
    """
    Base class. Subclass and override MCSBondCompare.__call__() to define custom bond compare functions, then set MCSParameters.BondTyper to an instance of the subclass

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def CheckBondRingMatch(
        self: MCSBondCompare,
        parameters: MCSBondCompareParameters,
        mol1: Mol,
        bond1: int,
        mol2: Mol,
        bond2: int,
    ) -> bool:
        """
        Return True if both bonds are, or are not, part of a ring

        C++ signature :bool CheckBondRingMatch(RDKit::PyMCSBondCompare {lvalue},RDKit::MCSBondCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    def CheckBondStereo(
        self: MCSBondCompare,
        parameters: MCSBondCompareParameters,
        mol1: Mol,
        bond1: int,
        mol2: Mol,
        bond2: int,
    ) -> bool:
        """
        Return True if both bonds have, or have not, a stereo descriptor

        C++ signature :bool CheckBondStereo(RDKit::PyMCSBondCompare {lvalue},RDKit::MCSBondCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    def compare(
        self: MCSBondCompare,
        parameters: MCSBondCompareParameters,
        mol1: Mol,
        bond1: int,
        mol2: Mol,
        bond2: int,
    ) -> bool:
        """
        override to implement custom bond comparison

        C++ signature :bool compare(RDKit::PyMCSBondCompare {lvalue},RDKit::MCSBondCompareParameters,RDKit::ROMol,unsigned int,RDKit::ROMol,unsigned int)
        """
        ...
    @classmethod
    def __call__(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MCSBondCompareParameters(Boost.Python.instance):
    """
    Parameters controlling how bond-bond matching is done

    C++ signature :void __init__(_object*)

    property CompleteRingsOnly¶
    results cannot include partial rings

    property MatchFusedRings¶
    enforce check on ring fusion, i.e. alpha-methylnaphthalene won’t match beta-methylnaphtalene, but decalin will match cyclodecane unless MatchFusedRingsStrict is True

    property MatchFusedRingsStrict¶
    only enforced if MatchFusedRings is True; the ring fusion must be the same in both query and target, i.e. decalin won’t match cyclodecane

    property MatchStereo¶
    include bond stereo in the comparison

    property RingMatchesRingOnly¶
    ring bonds are only allowed to match other ring bonds"""

    ...
    __instance_size__: ClassVar[int] = ...
    CompleteRingsOnly: Any
    MatchFusedRings: Any
    MatchFusedRingsStrict: Any
    MatchStereo: Any
    RingMatchesRingOnly: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MCSParameters(Boost.Python.instance):
    """
    Parameters controlling how the MCS is constructed

    C++ signature :void __init__(_object*)

    property AtomCompareParameters¶
    parameters for comparing atoms

    property AtomTyper¶
    atom typer to be used. Must be one of the members of the rdFMCS.AtomCompare class or an instance of a user-defined subclass of rdFMCS.MCSAtomCompare

    property BondCompareParameters¶
    parameters for comparing bonds

    property BondTyper¶
    bond typer to be used. Must be one of the members of the rdFMCS.BondCompare class or an instance of a user-defined subclass of rdFMCS.MCSBondCompare

    property InitialSeed¶
    SMILES string to be used as the seed of the MCS

    property MaximizeBonds¶
    toggles maximizing the number of bonds (instead of the number of atoms)

    property ProgressCallback¶
    progress callback class. Must be a user-defined subclass of rdFMCS.Progress"""

    __instance_size__: ClassVar[int] = ...
    AtomCompareParameters: Any
    AtomTyper: Any
    BondCompareParameters: Any
    BondTyper: Any
    InitialSeed: Any
    MaximizeBonds: Any
    ProgressCallback: Any
    Threshold: Any
    Timeout: Any
    Verbose: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def SetAtomTyper(self: MCSParameters, comparator: object) -> None:
        """
        DEPRECATED: please use the AtomTyper property instead. Sets the atom typer to be used. The argument must be one of the members of the rdFMCS.MCSAtomCompare class or an instance of a user-defined subclass of rdFMCS.MCSAtomCompare

        C++ signature :void SetAtomTyper(RDKit::PyMCSParameters {lvalue},_object*)"""
        ...
    def SetBondTyper(self: MCSParameters, comparator: object) -> None:
        """
        DEPRECATED: please use the BondTyper property instead. Sets the bond typer to be used. The argument must be one of the members of the rdFMCS.MCSBondCompare class or an instance of a user-defined subclass of rdFMCS.MCSBondCompare

        C++ signature :void SetBondTyper(RDKit::PyMCSParameters {lvalue},_object*)"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MCSProgress(Boost.Python.instance):
    """
    Base class. Subclass and override MCSProgress.__call__() to define a custom callback function

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def callback(self: MCSProgress, stat: object, parameters: object) -> bool:
        """
        DEPRECATED: override __call__ instead.
        override to implement a custom progress callback.

        C++ signature :bool callback(RDKit::PyMCSProgress {lvalue},RDKit::MCSProgressData,RDKit::MCSParameters)
        """
        ...
    @classmethod
    def __call__(cls, *args, **kwargs) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MCSProgressData(Boost.Python.instance):
    """
    Information about the MCS progress

    C++ signature :void __init__(_object*)

    property numAtoms¶
    number of atoms in MCS

    property numBonds¶
    number of bonds in MCS

    property seedProcessed¶
    number of processed seeds"""

    ...
    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def numAtoms(self) -> Any: ...
    @property
    def numBonds(self) -> Any: ...
    @property
    def seedProcessed(self) -> Any: ...

class MCSResult(Boost.Python.instance):
    """
    used to return MCS results
    Raises an exception
    This class cannot be instantiated from Python

    property canceled¶
    if True, the MCS calculation did not finish

    property numAtoms¶
    number of atoms in MCS

    property numBonds¶
    number of bonds in MCS

    property queryMol¶
    query molecule for the MCS

    property smartsString¶
    SMARTS string for the MCS"""

    ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def canceled(self) -> Any: ...
    @property
    def numAtoms(self) -> Any: ...
    @property
    def numBonds(self) -> Any: ...
    @property
    def queryMol(self) -> Any: ...
    @property
    def smartsString(self) -> Any: ...

class RingCompare(Boost.Python.enum):
    IgnoreRingFusion: RingCompare = ...
    PermissiveRingFusion: RingCompare = ...
    StrictRingFusion: RingCompare = ...
    names: dict[str, RingCompare] = ...
    values: dict[int, RingCompare] = ...
    __slots__: ClassVar[tuple] = ...

@overload
def FindMCS(
    self,
    mols: AtomPairsParameters,
    maximizeBonds: bool = True,
    threshold: float = 1.0,
    timeout: int = 3600,
    verbose: bool = False,
    matchValences: bool = False,
    ringMatchesRingOnly: bool = False,
    completeRingsOnly: bool = False,
    matchChiralTag: bool = False,
    atomCompare: AtomCompare = AtomCompare.CompareElements,
    bondCompare: BondCompare = BondCompare.CompareOrder,
    ringCompare: RingCompare = RingCompare.IgnoreRingFusion,
    seedSmarts: str = "",
) -> MCSResult:
    """
    Find the MCS for a set of molecules

    C++ signature :RDKit::MCSResult* FindMCS(boost::python::api::object,RDKit::PyMCSParameters {lvalue})
    """
    ...

@overload
def FindMCS(
    self, mols: AtomPairsParameters, parameters: MCSParameters
) -> MCSResult: ...
