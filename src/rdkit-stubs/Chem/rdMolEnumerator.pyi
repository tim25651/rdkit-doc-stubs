"""
rdkit.Chem.rdMolEnumerator module¶
Module containing classes and functions for enumerating molecules
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol, MolBundle

class EnumeratorType(Boost.Python.enum):
    LinkNode: EnumeratorType = ...
    PositionVariation: EnumeratorType = ...
    RepeatUnit: EnumeratorType = ...
    names: dict[str, EnumeratorType] = ...
    values: dict[int, EnumeratorType] = ...
    __slots__: ClassVar[tuple] = ...

class MolEnumeratorParams(Boost.Python.instance):
    """
    Molecular enumerator parameters

    C++ signature :void __init__(_object*)

    __init__( (AtomPairsParameters)arg1, (EnumeratorType)arg2) -> object :

    C++ signature :void* __init__(boost::python::api::object,(anonymous namespace)::EnumeratorTypes)
    """

    __instance_size__: ClassVar[int] = ...
    doRandom: Any
    maxToEnumerate: Any
    randomSeed: Any
    sanitize: Any

    @classmethod
    def __init__(cls, boost, anonymousnamespace) -> Any: ...
    def SetEnumerationOperator(
        self, arg1: MolEnumeratorParams, arg2: EnumeratorType
    ) -> None:
        """
        set the operator to be used for enumeration

        C++ signature :void SetEnumerationOperator(RDKit::MolEnumerator::MolEnumeratorParams*,(anonymous namespace)::EnumeratorTypes)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

@overload
def Enumerate(self, mol: Mol, maxPerOperation: int = 0) -> MolBundle:
    """
    do an enumeration for the supplied parameter type and return a MolBundle
    Limitations:

    Overlapping SRUs, i.e. where one monomer is contained within another, are
    not supported

    C++ signature :RDKit::MolBundle* Enumerate(RDKit::ROMol,RDKit::MolEnumerator::MolEnumeratorParams)
    """
    ...

@overload
def Enumerate(self, mol: Mol, enumParams: MolEnumeratorParams) -> MolBundle: ...
