"""
Module containing functions for interchange of molecules.
Note that this should be considered beta and that the format

and API will very likely change in future releases.
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class JSONParseParameters(Boost.Python.instance):
    """
    Parameters controlling the JSON parser

    C++ signature :void __init__(_object*)

    property parseConformers¶
    parse conformers in the JSON

    property parseProperties¶
    parse molecular properties in the JSON

    property setAromaticBonds¶
    set bond types to aromatic for bonds flagged aromatic

    property strictValenceCheck¶
    be strict when checking atom valences

    property useHCounts¶
    use atomic H counts from the JSON. You may want to set this to False when parsing queries.
    """

    ...
    __instance_size__: ClassVar[int] = ...
    parseConformers: Any
    parseProperties: Any
    setAromaticBonds: Any
    strictValenceCheck: Any
    useHCounts: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class JSONWriteParameters(Boost.Python.instance):
    """
    Parameters controlling the JSON writer

    C++ signature :void __init__(_object*)

    property useRDKitExtensions¶
    use RDKit extensions to the commonchem format"""

    ...
    __instance_size__: ClassVar[int] = ...
    useRDKitExtensions: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def JSONToMols(self, jsonBlock: str, params: AtomPairsParameters = None) -> tuple:
    """
    Convert JSON to a tuple of molecules

    ARGUMENTS:
    jsonBlock: the molecule to work with
    params: (optional) JSONParseParameters controlling the JSON parsing

    RETURNS:a tuple of Mols

    C++ signature :boost::python::tuple JSONToMols(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,boost::python::api::object=None])
    """
    ...

def MolToJSON(self, mol: Mol, params: AtomPairsParameters = None) -> str:
    """
    Convert a single molecule to JSON

    ARGUMENTS:
    mol: the molecule to work with

    RETURNS:a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToJSON(RDKit::ROMol [,boost::python::api::object=None])
    """
    ...

def MolsToJSON(
    self, mols: AtomPairsParameters, params: AtomPairsParameters = None
) -> str:
    """
    Convert a set of molecules to JSON

    ARGUMENTS:
    mols: the molecules to work with

    RETURNS:a string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolsToJSON(boost::python::api::object [,boost::python::api::object=None])
    """
    ...
