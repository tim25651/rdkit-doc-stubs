"""
rdkit.Chem.rdfragcatalog module¶
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.rdBase import _vectd, _vecti

class FragCatGenerator(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddFragsFromMol(
        self, arg1: FragCatGenerator, arg2: Mol, arg3: FragCatalog
    ) -> int:
        """
        C++ signature :unsigned int AddFragsFromMol(RDKit::FragCatGenerator {lvalue},RDKit::ROMol,RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int>*)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragCatParams(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*,int,int,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,double=1e-08])
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetFuncGroup(self, arg1: FragCatParams, arg2: int) -> Mol:
        """
        C++ signature :RDKit::ROMol const* GetFuncGroup(RDKit::FragCatParams {lvalue},int)
        """
        ...
    def GetLowerFragLength(self, arg1: FragCatParams) -> int:
        """
        C++ signature :unsigned int GetLowerFragLength(RDKit::FragCatParams {lvalue})"""
        ...
    def GetNumFuncGroups(self, arg1: FragCatParams) -> int:
        """
        C++ signature :unsigned int GetNumFuncGroups(RDKit::FragCatParams {lvalue})"""
        ...
    def GetTolerance(self, arg1: FragCatParams) -> float:
        """
        C++ signature :double GetTolerance(RDKit::FragCatParams {lvalue})"""
        ...
    def GetTypeString(self, arg1: FragCatParams) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetTypeString(RDKit::FragCatParams {lvalue})
        """
        ...
    def GetUpperFragLength(self, arg1: FragCatParams) -> int:
        """
        C++ signature :unsigned int GetUpperFragLength(RDKit::FragCatParams {lvalue})"""
        ...
    def Serialize(self, arg1: FragCatParams) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Serialize(RDKit::FragCatParams {lvalue})
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class FragCatalog(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*,RDKit::FragCatParams*)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetBitDescription(self, arg1: FragCatalog, arg2: int) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetBitDescription(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetBitDiscrims(self, arg1: FragCatalog, arg2: int) -> _vectd:
        """
        C++ signature :std::vector<double, std::allocator<double> > GetBitDiscrims(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetBitEntryId(self, arg1: FragCatalog, arg2: int) -> int:
        """
        C++ signature :unsigned int GetBitEntryId(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetBitFuncGroupIds(self, arg1: FragCatalog, arg2: int) -> _vecti:
        """
        C++ signature :std::vector<int, std::allocator<int> > GetBitFuncGroupIds(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetBitOrder(self, arg1: FragCatalog, arg2: int) -> int:
        """
        C++ signature :unsigned int GetBitOrder(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetCatalogParams(self, arg1: FragCatalog) -> FragCatParams:
        """
        C++ signature :RDKit::FragCatParams* GetCatalogParams(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
        ...
    def GetEntryBitId(self, arg1: FragCatalog, arg2: int) -> int:
        """
        C++ signature :unsigned int GetEntryBitId(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetEntryDescription(self, arg1: FragCatalog, arg2: int) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetEntryDescription(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetEntryDownIds(self, arg1: FragCatalog, arg2: int) -> _vecti:
        """
        C++ signature :std::vector<int, std::allocator<int> > GetEntryDownIds(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetEntryFuncGroupIds(self, arg1: FragCatalog, arg2: int) -> _vecti:
        """
        C++ signature :std::vector<int, std::allocator<int> > GetEntryFuncGroupIds(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetEntryOrder(self, arg1: FragCatalog, arg2: int) -> int:
        """
        C++ signature :unsigned int GetEntryOrder(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> const*,unsigned int)
        """
        ...
    def GetFPLength(self, arg1: FragCatalog) -> int:
        """
        C++ signature :unsigned int GetFPLength(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
        ...
    def GetNumEntries(self, arg1: FragCatalog) -> int:
        """
        C++ signature :unsigned int GetNumEntries(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
        ...
    def Serialize(self, arg1: FragCatalog) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Serialize(RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int> {lvalue})
        """
        ...
    @classmethod
    def __getinitargs__(cls) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

class FragFPGenerator(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetFPForMol(
        self, arg1: FragFPGenerator, arg2: Mol, arg3: FragCatalog
    ) -> ExplicitBitVect:
        """
        C++ signature :ExplicitBitVect* GetFPForMol(RDKit::FragFPGenerator {lvalue},RDKit::ROMol,RDCatalog::HierarchCatalog<RDKit::FragCatalogEntry, RDKit::FragCatParams, int>)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
