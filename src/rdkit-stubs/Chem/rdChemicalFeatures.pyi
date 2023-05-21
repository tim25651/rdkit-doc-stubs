"""
Module containing free chemical feature functionality
These are feature that are not associated with molecules. They are 
are typically derived from pharmacophores and site-maps.
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Geometry.rdGeometry import Point3D

class FreeChemicalFeature(Boost.Python.instance):
    """
    Class to represent a free chemical features.
    These chemical features are not associated with a molecule, though they can be matched
    to molecular featufres

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1) -> None :Default Constructor

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)family, (str)type, (Point3D)loc [, (int)id=-1]) -> None :Constructor with family, type and location specified

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point3D [,int=-1])

    __init__( (object)arg1, (str)family, (Point3D)loc) -> None :constructor with family and location specified, empty type and id

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point3D)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def GetFamily(self, arg1: FreeChemicalFeature) -> str:
        """
        Get the family of the feature

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetFamily(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
        ...
    def GetId(self, arg1: FreeChemicalFeature) -> int:
        """
        Get the id of the feature

        C++ signature :int GetId(ChemicalFeatures::FreeChemicalFeature {lvalue})"""
        ...
    def GetPos(self, arg1: FreeChemicalFeature) -> Point3D:
        """
        Get the position of the feature

        C++ signature :RDGeom::Point3D GetPos(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
        ...
    def GetType(self, arg1: FreeChemicalFeature) -> str:
        """
        Get the sepcific type for the feature

        C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > GetType(ChemicalFeatures::FreeChemicalFeature {lvalue})
        """
        ...
    def SetFamily(self, arg1: FreeChemicalFeature, arg2: str) -> None:
        """
        Set the family of the feature

        C++ signature :void SetFamily(ChemicalFeatures::FreeChemicalFeature {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    def SetId(self, arg1: FreeChemicalFeature, arg2: int) -> None:
        """
        Set the id of the feature

        C++ signature :void SetId(ChemicalFeatures::FreeChemicalFeature {lvalue},int)"""
        ...
    def SetPos(self, arg1: FreeChemicalFeature, arg2: Point3D) -> None:
        """
        Set the feature position

        C++ signature :void SetPos(ChemicalFeatures::FreeChemicalFeature {lvalue},RDGeom::Point3D)
        """
        ...
    def SetType(self, arg1: FreeChemicalFeature, arg2: str) -> None:
        """
        Set the sepcific type for the feature

        C++ signature :void SetType(ChemicalFeatures::FreeChemicalFeature {lvalue},std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
        """
        ...
    @classmethod
    def __getinitargs__(cls, ChemicalFeatures) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
