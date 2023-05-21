"""
rdkit.Chem.rdCoordGen module¶
Module containing interface to the CoordGen library.
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class CoordGenParams(Boost.Python.instance):
    """
    Parameters controlling coordinate generation

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...
    coordgenScaling: Any
    dbg_useConstrained: Any
    dbg_useFixed: Any
    minimizerPrecision: Any
    templateFileDir: Any
    treatNonterminalBondsToMetalAsZOBs: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def SetCoordMap(self, arg1: CoordGenParams, arg2: dict) -> None:
        """
        expects a dictionary of Point2D objects with template coordinates

        C++ signature :void SetCoordMap(RDKit::CoordGen::CoordGenParams*,boost::python::dict {lvalue})
        """
        ...
    def SetTemplateMol(self, arg1: CoordGenParams, arg2: Mol) -> None:
        """
        sets a molecule to be used as the template

        C++ signature :void SetTemplateMol(RDKit::CoordGen::CoordGenParams*,RDKit::ROMol const*)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def sketcherBestPrecision(self) -> Any: ...
    @property
    def sketcherCoarsePrecision(self) -> Any: ...
    @property
    def sketcherQuickPrecision(self) -> Any: ...
    @property
    def sketcherStandardPrecision(self) -> Any: ...

def AddCoords(self, mol: Mol, params: AtomPairsParameters = None) -> None:
    """
    Add 2D coordinates.
    ARGUMENTS:

    mol: molecule to modify
    params: (optional) parameters controlling the coordinate generation

    C++ signature :void AddCoords(RDKit::ROMol {lvalue} [,boost::python::api::object {lvalue}=None])
    """
    ...

def SetDefaultTemplateFileDir(self, arg1: str) -> None:
    """
    C++ signature :void SetDefaultTemplateFileDir(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...
