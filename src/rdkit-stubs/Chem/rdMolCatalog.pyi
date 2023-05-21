"""
rdkit.Chem.rdMolCatalog module¶
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.rdBase import _vecti

class MolCatalog(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AddEdge(self, arg1: MolCatalog, arg2: int, arg3: int) -> None:
        """
        C++ signature :void AddEdge(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue},unsigned int,unsigned int)
        """
        ...
    def AddEntry(self, arg1: MolCatalog, arg2: MolCatalogEntry) -> int:
        """
        C++ signature :unsigned int AddEntry(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int>*,RDKit::MolCatalogEntry*)
        """
        ...
    def GetBitDescription(self, arg1: MolCatalog, arg2: int) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetBitDescription(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
        ...
    def GetBitEntryId(self, arg1: MolCatalog, arg2: int) -> int:
        """
        C++ signature :unsigned int GetBitEntryId(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
        ...
    def GetEntryBitId(self, arg1: MolCatalog, arg2: int) -> int:
        """
        C++ signature :unsigned int GetEntryBitId(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
        ...
    def GetEntryDescription(self, arg1: MolCatalog, arg2: int) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetEntryDescription(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
        ...
    def GetEntryDownIds(self, arg1: MolCatalog, arg2: int) -> _vecti:
        """
        C++ signature :std::vector<int, std::allocator<int> > GetEntryDownIds(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> const*,unsigned int)
        """
        ...
    def GetFPLength(self, arg1: MolCatalog) -> int:
        """
        C++ signature :unsigned int GetFPLength(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue})
        """
        ...
    def GetNumEntries(self, arg1: MolCatalog) -> int:
        """
        C++ signature :unsigned int GetNumEntries(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue})
        """
        ...
    def Serialize(self, arg1: MolCatalog) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > Serialize(RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int> {lvalue})
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

class MolCatalogEntry(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetDescription(self, arg1: MolCatalogEntry) -> str:
        """
        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetDescription(RDKit::MolCatalogEntry {lvalue})
        """
        ...
    def GetMol(self, arg1: MolCatalogEntry) -> Mol:
        """
        C++ signature :RDKit::ROMol GetMol(RDKit::MolCatalogEntry {lvalue})"""
        ...
    def GetOrder(self, arg1: MolCatalogEntry) -> int:
        """
        C++ signature :unsigned int GetOrder(RDKit::MolCatalogEntry {lvalue})"""
        ...
    def SetDescription(self, arg1: MolCatalogEntry, arg2: str) -> None:
        """
        C++ signature :void SetDescription(RDKit::MolCatalogEntry {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetMol(self, arg1: MolCatalogEntry, arg2: Mol) -> None:
        """
        C++ signature :void SetMol(RDKit::MolCatalogEntry*,RDKit::ROMol const*)"""
        ...
    def SetOrder(self, arg1: MolCatalogEntry, arg2: int) -> None:
        """
        C++ signature :void SetOrder(RDKit::MolCatalogEntry {lvalue},unsigned int)"""
        ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

def CreateMolCatalog(self) -> MolCatalog:
    """
    C++ signature :RDCatalog::HierarchCatalog<RDKit::MolCatalogEntry, RDKit::MolCatalogParams, int>* CreateMolCatalog()
    """
    ...
