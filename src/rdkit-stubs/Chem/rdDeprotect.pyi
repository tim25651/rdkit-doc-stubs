"""
rdkit.Chem.rdDeprotect module¶
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class DeprotectData(Boost.Python.instance):
    """
    DeprotectData class, contains a single deprotection reaction and information
    deprotectdata.deprotection_class - functional group being protected
    deprotectdata.reaction_smarts - reaction smarts used for deprotection
    deprotectdata.abbreviation - common abbreviation for the protecting group
    deprotectdata.full_name - full name for the protecting group

    Construct a new DeprotectData instance.>>> reaction_class = "amine"
    >>> reaction_smarts = "[C;R0][C;R0]([C;R0])([O;R0][C;R0](=[O;R0])[NX3;H0,H1:1])C>>[N:1]"
    >>> abbreviation = "Boc"
    >>> full_name = "tert-butyloxycarbonyl"
    >>> data = DeprotectData(reaction_class, reaction_smarts, abbreviation, full_name)
    >>> assert data.isValid()

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    property abbreviation¶

    property deprotection_class¶

    property example¶

    property full_name¶"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def isValid(self, arg1: DeprotectData) -> bool:
        """
        Returns True if the DeprotectData has a valid reaction

        C++ signature :bool isValid(RDKit::Deprotect::DeprotectData {lvalue})"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def abbreviation(self) -> Any: ...
    @property
    def deprotection_class(self) -> Any: ...
    @property
    def example(self) -> Any: ...
    @property
    def full_name(self) -> Any: ...
    @property
    def reaction_smarts(self) -> Any: ...

class DeprotectDataVect(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: DeprotectDataVect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: DeprotectDataVect, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > {lvalue},boost::python::api::object)
        """
        ...
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

def Deprotect(self, mol: Mol, deprotections: AtomPairsParameters = None) -> Mol:
    """
    Return the deprotected version of the molecule.

    C++ signature :boost::shared_ptr<RDKit::ROMol> Deprotect(RDKit::ROMol [,boost::python::api::object=None])
    """
    ...

def DeprotectInPlace(self, mol: Mol, deprotections: AtomPairsParameters = None) -> bool:
    """
    Deprotects the molecule in place.

    C++ signature :bool DeprotectInPlace(RDKit::ROMol {lvalue} [,boost::python::api::object=None])
    """
    ...

def GetDeprotections(self) -> DeprotectDataVect:
    """
    Return the default list of deprotections

    C++ signature :std::vector<RDKit::Deprotect::DeprotectData, std::allocator<RDKit::Deprotect::DeprotectData> > GetDeprotections()
    """
    ...
