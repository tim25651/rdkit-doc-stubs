"""
rdkit.Chem.rdTautomerQuery module¶
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol, SubstructMatchParameters
from rdkit.Chem.rdChemReactions import VectSizeT
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

class TautomerQuery(Boost.Python.instance):
    """
    The Tautomer Query Class.
    Creates a query that enables structure search accounting for matching of
    Tautomeric forms

    C++ signature :void* __init__(boost::python::api::object,RDKit::ROMol)

    __init__( (AtomPairsParameters)arg1, (Mol)arg2, (str)arg3) -> object :

    C++ signature :void* __init__(boost::python::api::object,RDKit::ROMol,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @overload
    @classmethod
    def __init__(cls, boost, RDKit) -> Any: ...
    @overload
    @classmethod
    def __init__(cls, boost, RDKit, std) -> Any: ...
    def GetModifiedAtoms(self, arg1: TautomerQuery) -> VectSizeT:
        """
        C++ signature :std::vector<unsigned long, std::allocator<unsigned long> > GetModifiedAtoms(RDKit::TautomerQuery {lvalue})
        """
        ...
    def GetModifiedBonds(self, arg1: TautomerQuery) -> VectSizeT:
        """
        C++ signature :std::vector<unsigned long, std::allocator<unsigned long> > GetModifiedBonds(RDKit::TautomerQuery {lvalue})
        """
        ...
    @overload
    def GetSubstructMatch(
        self: TautomerQuery,
        target: Mol,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatch(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatch(
        self: TautomerQuery, target: Mol, params: SubstructMatchParameters
    ) -> object: ...
    @overload
    def GetSubstructMatches(
        self: TautomerQuery,
        target: Mol,
        uniquify: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatches(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatches(
        self: TautomerQuery, target: Mol, params: SubstructMatchParameters
    ) -> object: ...
    @overload
    def GetSubstructMatchesWithTautomers(
        self: TautomerQuery,
        target: Mol,
        uniquify: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
        maxMatches: int = 1000,
    ) -> object:
        """


        C++ signature :_object* GetSubstructMatchesWithTautomers(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def GetSubstructMatchesWithTautomers(
        self: TautomerQuery, target: Mol, params: SubstructMatchParameters
    ) -> object: ...
    def GetTautomers(self, arg1: TautomerQuery) -> object:
        """
        C++ signature :_object* GetTautomers(RDKit::TautomerQuery)"""
        ...
    def GetTemplateMolecule(self, arg1: TautomerQuery) -> Mol:
        """
        C++ signature :RDKit::ROMol GetTemplateMolecule(RDKit::TautomerQuery {lvalue})
        """
        ...
    @overload
    def IsSubstructOf(
        self: TautomerQuery,
        target: Mol,
        recursionPossible: bool = True,
        useChirality: bool = False,
        useQueryQueryMatches: bool = False,
    ) -> bool:
        """


        C++ signature :bool IsSubstructOf(RDKit::TautomerQuery,RDKit::ROMol,RDKit::SubstructMatchParameters)
        """
        ...
    @overload
    def IsSubstructOf(
        self: TautomerQuery, target: Mol, params: SubstructMatchParameters
    ) -> bool: ...
    def PatternFingerprintTemplate(
        self, arg1: TautomerQuery, fingerprintSize: int = 2048
    ) -> ExplicitBitVect:
        """
        C++ signature :ExplicitBitVect* PatternFingerprintTemplate(RDKit::TautomerQuery {lvalue} [,unsigned int=2048])
        """
        ...
    def ToBinary(self, arg1: TautomerQuery) -> object:
        """
        C++ signature :boost::python::api::object ToBinary(RDKit::TautomerQuery)"""
        ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

def PatternFingerprintTautomerTarget(
    self, target: Mol, fingerprintSize: int = 2048
) -> ExplicitBitVect:
    """
    C++ signature :ExplicitBitVect* PatternFingerprintTautomerTarget(RDKit::ROMol [,unsigned int=2048])
    """
    ...

def TautomerQueryCanSerialize(self) -> bool:
    """
    Returns True if the TautomerQuery is serializable (requires that the RDKit was built with boost::serialization)

    C++ signature :bool TautomerQueryCanSerialize()"""
    ...
