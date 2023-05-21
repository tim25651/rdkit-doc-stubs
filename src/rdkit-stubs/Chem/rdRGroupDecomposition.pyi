"""
Module containing RGroupDecomposition classes and functions.
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

AtomIndexLabels: RGroupLabels
AtomMap: RGroupLabelling
AtomMapLabels: RGroupLabels
AutoDetect: RGroupLabels
DummyAtomLabels: RGroupLabels
Exhaustive: RGroupMatching
FingerprintVariance: RGroupScore
GA: RGroupMatching
Greedy: RGroupMatching
GreedyChunks: RGroupMatching
Isotope: RGroupLabelling
IsotopeLabels: RGroupLabels
MCS: RGroupCoreAlignment
MDLRGroup: RGroupLabelling
MDLRGroupLabels: RGroupLabels
Match: RGroupScore
NoAlignment: RGroupCoreAlignment
NoSymmetrization: RGroupMatching
None_: RGroupCoreAlignment
RelabelDuplicateLabels: RGroupLabels

class RGroupCoreAlignment(Boost.Python.enum):
    MCS: RGroupCoreAlignment = ...
    NoAlignment: RGroupCoreAlignment = ...
    None_: RGroupCoreAlignment = ...
    names: dict[str, RGroupCoreAlignment] = ...
    values: dict[int, RGroupCoreAlignment] = ...
    __slots__: ClassVar[tuple] = ...

class RGroupDecomposition(Boost.Python.instance):
    """
    RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures
    OPTIONS:

    RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or RGroupCoreAlignment.MCSIf set to MCS, cores labels are mapped to each other using their
    Maximum common substructure overlap.

    RGroupLabels: optionally set where the rgroup labels to use are encoded.
    RGroupLabels.IsotopeLabels - labels are stored on isotopes
    RGroupLabels.AtomMapLabels - labels are stored on atommaps
    RGroupLabels.MDLRGroupLabels - labels are stored on MDL R-groups
    RGroupLabels.DummyAtomLabels - labels are stored on dummy atoms
    RGroupLabels.AtomIndexLabels - use the atom index as the label
    RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels
    RGroupLabels.AutoDetect - auto detect the label [default]

    Note: in all cases, any rgroups found on unlabelled atoms will be automaticallylabelled.

    RGroupLabelling: choose where the rlabels are stored on the decomposition
    RGroupLabelling.AtomMap - store rgroups as atom maps (for smiles)
    RGroupLabelling.Isotope - store rgroups on the isotope
    RGroupLabelling.MDLRGroup - store rgroups as mdl rgroups (for molblocks)

    default: AtomMap | MDLRGroup

    onlyMatchAtRGroups: only allow rgroup decomposition at the specified rgroups
    removeAllHydrogenRGroups: remove all user-defined rgroups that only have hydrogens
    removeAllHydrogenRGroupsAndLabels: remove all user-defined rgroups that only have hydrogens, and also remove the corresponding labels from the core
    removeHydrogensPostMatch: remove all hydrogens from the output molecules
    allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or more

    Construct from a molecule or sequence of molecules

    C++ signature :void __init__(_object*,boost::python::api::object)

    __init__( (object)arg1, (AtomPairsParameters)arg2, (RGroupDecompositionParameters)arg3) -> None :Construct from a molecule or sequence of molecules and a parameters object

    C++ signature :void __init__(_object*,boost::python::api::object,RDKit::RGroupDecompositionParameters)
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Add(self, arg1: RGroupDecomposition, arg2: Mol) -> int:
        """
        C++ signature :int Add(RDKit::RGroupDecompositionHelper {lvalue},RDKit::ROMol)
        """
        ...
    def GetRGroupLabels(self, arg1: RGroupDecomposition) -> list:
        """
        Return the current list of found rgroups.
        Note, Process() should be called first

        C++ signature :boost::python::list GetRGroupLabels(RDKit::RGroupDecompositionHelper {lvalue})
        """
        ...
    def GetRGroupsAsColumns(
        self, arg1: RGroupDecomposition, asSmiles: bool = False
    ) -> dict:
        """
        Return the rgroups as columns (note: can be fed directrly into a pandas datatable)
        ARGUMENTS:
        asSmiles: if True return smiles strings, otherwise return molecules [default: False]

        Column structure:columns[rgroup_label] = [ mols_or_smiles ]

        C++ signature :boost::python::dict GetRGroupsAsColumns(RDKit::RGroupDecompositionHelper {lvalue} [,bool=False])
        """
        ...
    def GetRGroupsAsRows(
        self, arg1: RGroupDecomposition, asSmiles: bool = False
    ) -> list:
        """
        Return the rgroups as rows (note: can be fed directrly into a pandas datatable)
        ARGUMENTS:
        asSmiles: if True return smiles strings, otherwise return molecules [default: False]

        Row structure:rows[idx] = {rgroup_label: molecule_or_smiles}

        C++ signature :boost::python::list GetRGroupsAsRows(RDKit::RGroupDecompositionHelper {lvalue} [,bool=False])
        """
        ...
    def Process(self, arg1: RGroupDecomposition) -> bool:
        """
        Process the rgroups (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)

        C++ signature :bool Process(RDKit::RGroupDecompositionHelper {lvalue})"""
        ...
    def ProcessAndScore(self, arg1: RGroupDecomposition) -> tuple:
        """
        Process the rgroups and returns the score (must be done prior to GetRGroupsAsRows/Columns and GetRGroupLabels)

        C++ signature :boost::python::tuple ProcessAndScore(RDKit::RGroupDecompositionHelper {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class RGroupDecompositionParameters(Boost.Python.instance):
    """
    RGroupDecompositionParameters controls how the RGroupDecomposition sets labelling and matches structures
    OPTIONS:

    RGroupCoreAlignment: can be one of RGroupCoreAlignment.None_ or RGroupCoreAlignment.MCSIf set to MCS, cores labels are mapped to each other using their
    Maximum common substructure overlap.

    RGroupLabels: optionally set where the rgroup labels to use are encoded.
    RGroupLabels.IsotopeLabels - labels are stored on isotopes
    RGroupLabels.AtomMapLabels - labels are stored on atommaps
    RGroupLabels.MDLRGroupLabels - labels are stored on MDL R-groups
    RGroupLabels.DummyAtomLabels - labels are stored on dummy atoms
    RGroupLabels.AtomIndexLabels - use the atom index as the label
    RGroupLabels.RelabelDuplicateLabels - fix any duplicate labels
    RGroupLabels.AutoDetect - auto detect the label [default]

    Note: in all cases, any rgroups found on unlabelled atoms will be automaticallylabelled.

    RGroupLabelling: choose where the rlabels are stored on the decomposition
    RGroupLabelling.AtomMap - store rgroups as atom maps (for smiles)
    RGroupLabelling.Isotope - store rgroups on the isotope
    RGroupLabelling.MDLRGroup - store rgroups as mdl rgroups (for molblocks)

    default: AtomMap | MDLRGroup

    onlyMatchAtRGroups: only allow rgroup decomposition at the specified rgroups
    removeAllHydrogenRGroups: remove all user-defined rgroups that only have hydrogens
    removeAllHydrogenRGroupsAndLabels: remove all user-defined rgroups that only have hydrogens, and also remove the corresponding labels from the core
    removeHydrogensPostMatch: remove all hydrogens from the output molecules
    allowNonTerminalRGroups: allow labelled Rgroups of degree 2 or more

    Constructor, takes no arguments

    C++ signature :void __init__(_object*)

    property alignment¶

    property allowMultipleRGroupsOnUnlabelled¶

    property allowNonTerminalRGroups¶

    property chunkSize¶

    property gaMaximumOperations¶

    property gaNumberOperationsWithoutImprovement¶

    property gaNumberRuns¶

    property gaParallelRuns¶

    property gaPopulationSize¶

    property gaRandomSeed¶

    property labels¶

    property matchingStrategy¶

    property onlyMatchAtRGroups¶

    property removeAllHydrogenRGroups¶

    property removeAllHydrogenRGroupsAndLabels¶

    property removeHydrogensPostMatch¶

    property rgroupLabelling¶

    property scoreMethod¶

    property substructMatchParams¶

    property timeout¶"""

    ...
    __instance_size__: ClassVar[int] = ...
    alignment: Any
    allowMultipleRGroupsOnUnlabelled: Any
    allowNonTerminalRGroups: Any
    chunkSize: Any
    gaMaximumOperations: Any
    gaNumberOperationsWithoutImprovement: Any
    gaNumberRuns: Any
    gaParallelRuns: Any
    gaPopulationSize: Any
    gaRandomSeed: Any
    labels: Any
    matchingStrategy: Any
    onlyMatchAtRGroups: Any
    removeAllHydrogenRGroups: Any
    removeAllHydrogenRGroupsAndLabels: Any
    removeHydrogensPostMatch: Any
    rgroupLabelling: Any
    scoreMethod: Any
    timeout: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def substructMatchParams(self) -> Any: ...

class RGroupLabelling(Boost.Python.enum):
    AtomMap: RGroupLabelling = ...
    Isotope: RGroupLabelling = ...
    MDLRGroup: RGroupLabelling = ...
    names: dict[str, RGroupLabelling] = ...
    values: dict[int, RGroupLabelling] = ...
    __slots__: ClassVar[tuple] = ...

class RGroupLabels(Boost.Python.enum):
    AtomIndexLabels: RGroupLabels = ...
    AtomMapLabels: RGroupLabels = ...
    AutoDetect: RGroupLabels = ...
    DummyAtomLabels: RGroupLabels = ...
    IsotopeLabels: RGroupLabels = ...
    MDLRGroupLabels: RGroupLabels = ...
    RelabelDuplicateLabels: RGroupLabels = ...
    names: dict[str, RGroupLabels] = ...
    values: dict[int, RGroupLabels] = ...
    __slots__: ClassVar[tuple] = ...

class RGroupMatching(Boost.Python.enum):
    Exhaustive: RGroupMatching = ...
    GA: RGroupMatching = ...
    Greedy: RGroupMatching = ...
    GreedyChunks: RGroupMatching = ...
    NoSymmetrization: RGroupMatching = ...
    names: dict[str, RGroupMatching] = ...
    values: dict[int, RGroupMatching] = ...
    __slots__: ClassVar[tuple] = ...

class RGroupScore(Boost.Python.enum):
    FingerprintVariance: RGroupScore = ...
    Match: RGroupScore = ...
    names: dict[str, RGroupScore] = ...
    values: dict[int, RGroupScore] = ...
    __slots__: ClassVar[tuple] = ...

def RGroupDecompose(
    self,
    cores: AtomPairsParameters,
    mols: AtomPairsParameters,
    asSmiles: bool = False,
    asRows: bool = True,
    options: RGroupDecompositionParameters = ...,
) -> object:
    """
    Decompose a collecion of molecules into their Rgroups
    ARGUMENTS:

    cores: a set of cores from most to least specific.See RGroupDecompositionParameters for more details
    on how the cores can be labelled

    mols: the molecules to be decomposed
    asSmiles: if True return smiles strings, otherwise return molecules [default: False]
    asRows: return the results as rows (default) otherwise return columns

    RETURNS: row_or_column_results, unmatched

    Row structure:rows[idx] = {rgroup_label: molecule_or_smiles}

    Column structure:columns[rgroup_label] = [ mols_or_smiles ]

    unmatched is a vector of indices in the input mols that were not matched.

    C++ signature :boost::python::api::object RGroupDecompose(boost::python::api::object,boost::python::api::object [,bool=False [,bool=True [,RDKit::RGroupDecompositionParameters=<rdkit.Chem.rdRGroupDecomposition.RGroupDecompositionParameters object at 0x7fd468dcc2e0>]]])
    """
    ...
