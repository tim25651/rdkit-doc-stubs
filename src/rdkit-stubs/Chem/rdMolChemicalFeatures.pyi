"""
Module containing from chemical feature and functions to generate the
"""
from typing import Any, ClassVar, overload

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.Geometry.rdGeometry import Point3D

class MolChemicalFeature(Boost.Python.instance):
    """
    Class to represent a chemical feature.
    These chemical features may or may not have been derived from molecule object;
    i.e. it is possible to have a chemical feature that was created just from its type
    and location.
    Raises an exception
    This class cannot be instantiated from Python"""

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def ClearCache(self, arg1: MolChemicalFeature) -> None:
        """
        Clears the cache used to store position information.

        C++ signature :void ClearCache(RDKit::MolChemicalFeature {lvalue})"""
        ...
    def GetActiveConformer(self, arg1: MolChemicalFeature) -> int:
        """
        Gets the conformer to use.

        C++ signature :int GetActiveConformer(RDKit::MolChemicalFeature {lvalue})"""
        ...
    def GetAtomIds(self, arg1: MolChemicalFeature) -> object:
        """
        Get the IDs of the atoms that participate in the feature

        C++ signature :_object* GetAtomIds(RDKit::MolChemicalFeature)"""
        ...
    def GetFactory(self, arg1: MolChemicalFeature) -> MolChemicalFeatureFactory:
        """
        Get the factory used to generate this feature

        C++ signature :RDKit::MolChemicalFeatureFactory const* GetFactory(RDKit::MolChemicalFeature {lvalue})
        """
        ...
    def GetFamily(self, arg1: MolChemicalFeature) -> str:
        """
        Get the family to which the feature belongs; donor, acceptor, etc.

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetFamily(RDKit::MolChemicalFeature {lvalue})
        """
        ...
    def GetId(self, arg1: MolChemicalFeature) -> int:
        """
        Returns the identifier of the feature

        C++ signature :int GetId(RDKit::MolChemicalFeature {lvalue})"""
        ...
    def GetMol(self, arg1: MolChemicalFeature) -> Mol:
        """
        Get the molecule used to derive the features

        C++ signature :RDKit::ROMol const* GetMol(RDKit::MolChemicalFeature {lvalue})"""
        ...
    @overload
    def GetPos(self: MolChemicalFeature, confId: int) -> Point3D:
        """
        Get the location of the default chemical feature (first position)

        C++ signature :RDGeom::Point3D GetPos(RDKit::MolChemicalFeature {lvalue})"""
        ...
    @overload
    def GetPos(self: MolChemicalFeature) -> Point3D: ...
    def GetType(self, arg1: MolChemicalFeature) -> str:
        """
        Get the specific type for the feature

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetType(RDKit::MolChemicalFeature {lvalue})
        """
        ...
    def SetActiveConformer(self, arg1: MolChemicalFeature, arg2: int) -> None:
        """
        Sets the conformer to use (must be associated with a molecule).

        C++ signature :void SetActiveConformer(RDKit::MolChemicalFeature {lvalue},int)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MolChemicalFeatureFactory(Boost.Python.instance):
    """
    Class to featurize a molecule
    Raises an exception
    This class cannot be instantiated from Python"""

    def GetFeaturesForMol(self, mol, includeOnly="", confId=-1): ...
    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetFeatureDefs(self, arg1: MolChemicalFeatureFactory) -> dict:
        """
        Get a dictionary with SMARTS definitions for each feature type

        C++ signature :boost::python::dict GetFeatureDefs(RDKit::MolChemicalFeatureFactory)
        """
        ...
    def GetFeatureFamilies(self, arg1: MolChemicalFeatureFactory) -> tuple:
        """
        Get a tuple of feature types

        C++ signature :boost::python::tuple GetFeatureFamilies(RDKit::MolChemicalFeatureFactory)
        """
        ...
    def GetFeaturesForMol(self, mol, includeOnly="", confId=-1): ...
    def GetMolFeature(
        self,
        arg1: MolChemicalFeatureFactory,
        mol: Mol,
        idx: int,
        includeOnly: str = "",
        recompute: bool = True,
        confId: int = -1,
    ) -> MolChemicalFeature:
        """
        returns a particular feature (by index)

        C++ signature :boost::shared_ptr<RDKit::MolChemicalFeature> GetMolFeature(RDKit::MolChemicalFeatureFactory,RDKit::ROMol,int [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’ [,bool=True [,int=-1]]])
        """
        ...
    def GetNumFeatureDefs(self, arg1: MolChemicalFeatureFactory) -> int:
        """
        Get the number of feature definitions

        C++ signature :int GetNumFeatureDefs(RDKit::MolChemicalFeatureFactory {lvalue})
        """
        ...
    def GetNumMolFeatures(
        self, arg1: MolChemicalFeatureFactory, mol: Mol, includeOnly: str = ""
    ) -> int:
        """
        Get the number of features the molecule has

        C++ signature :int GetNumMolFeatures(RDKit::MolChemicalFeatureFactory,RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def BuildFeatureFactory(self, arg1: str) -> MolChemicalFeatureFactory:
    """
    Construct a feature factory given a feature definition in a file

    C++ signature :RDKit::MolChemicalFeatureFactory* BuildFeatureFactory(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def BuildFeatureFactoryFromString(self, arg1: str) -> MolChemicalFeatureFactory:
    """
    Construct a feature factory given a feature definition block

    C++ signature :RDKit::MolChemicalFeatureFactory* BuildFeatureFactoryFromString(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def GetAtomMatch(self, featMatch: AtomPairsParameters, maxAts: int = 1024) -> object:
    """
    Returns an empty list if any of the features passed in share an atom.Otherwise a list of lists of atom indices is returned.

    C++ signature :boost::python::api::object GetAtomMatch(boost::python::api::object [,int=1024])
    """
    ...
