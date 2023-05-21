"""
rdkit.Chem.rdAbbreviations module¶
Module containing functions for working with molecular abbreviations
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdAbbreviations import (
    _vectN5RDKit13Abbreviations22AbbreviationDefinitionE,
)
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class AbbreviationDefinition(Boost.Python.instance):
    """
    Abbreviation Definition

    C++ signature :void __init__(_object*)

    property displayLabel¶
    the label in a drawing when the bond comes from the right

    property displayLabelW¶
    the label in a drawing when the bond comes from the west

    property label¶
    the label

    property mol¶
    the query molecule (should have a dummy as the first atom)"""

    ...
    __instance_size__: ClassVar[int] = ...
    displayLabel: Any
    displayLabelW: Any
    label: Any
    mol: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class _vectN5RDKit13Abbreviations22AbbreviationDefinitionE(Boost.Python.instance):
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

def CondenseAbbreviationSubstanceGroups(self, mol: Mol) -> Mol:
    """
    Finds and replaces abbreviation (i.e. “SUP”) substance groups in a molecule. The result is not sanitized.

    C++ signature :RDKit::ROMol* CondenseAbbreviationSubstanceGroups(RDKit::ROMol const*)
    """
    ...

def CondenseMolAbbreviations(
    self,
    mol: Mol,
    abbrevs: AtomPairsParameters,
    maxCoverage: float = 0.4,
    sanitize: bool = True,
) -> Mol:
    """
    Finds and replaces abbreviations in a molecule. The result is not sanitized.

    C++ signature :RDKit::ROMol* CondenseMolAbbreviations(RDKit::ROMol const*,boost::python::api::object [,double=0.4 [,bool=True]])
    """
    ...

def GetDefaultAbbreviations(
    self,
) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    returns a list of the default abbreviation definitions

    C++ signature :std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > GetDefaultAbbreviations()
    """
    ...

def GetDefaultLinkers(self) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    returns a list of the default linker definitions

    C++ signature :std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > GetDefaultLinkers()
    """
    ...

def LabelMolAbbreviations(
    self, mol: Mol, abbrevs: AtomPairsParameters, maxCoverage: float = 0.4
) -> Mol:
    """
    Finds abbreviations and adds to them to a molecule as “SUP” SubstanceGroups

    C++ signature :RDKit::ROMol* LabelMolAbbreviations(RDKit::ROMol const*,boost::python::api::object [,double=0.4])
    """
    ...

def ParseAbbreviations(
    self,
    text: str,
    removeExtraDummies: bool = False,
    allowConnectionToDummies: bool = False,
) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    returns a set of abbreviation definitions from a string

    C++ signature :std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > ParseAbbreviations(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=False [,bool=False]])
    """
    ...

def ParseLinkers(
    self, text: str
) -> _vectN5RDKit13Abbreviations22AbbreviationDefinitionE:
    """
    returns a set of linker definitions from a string

    C++ signature :std::vector<RDKit::Abbreviations::AbbreviationDefinition, std::allocator<RDKit::Abbreviations::AbbreviationDefinition> > ParseLinkers(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...
