"""
rdkit.Geometry.rdGeometry module¶
Module containing geometry objects like points, grids, etc
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.DataStructs.cDataStructs import DiscreteValueType, DiscreteValueVect

class Point2D(Boost.Python.instance):
    """
    A class to represent a two-dimensional point
    Default Constructor

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (float)arg2, (float)arg3) -> None :

    C++ signature :void __init__(_object*,double,double)

    __init__( (object)self, (Point3D)other) -> None :construct from a Point3D (ignoring the z component)

    C++ signature :void __init__(_object*,RDGeom::Point3D)"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...
    x: Any
    y: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AngleTo(self, arg1: Point2D, arg2: Point2D) -> float:
        """
        determines the angle between a vector to this point (between 0 and PI)

        C++ signature :double AngleTo(RDGeom::Point2D {lvalue},RDGeom::Point2D)"""
        ...
    def DirectionVector(self, arg1: Point2D, arg2: Point2D) -> Point2D:
        """
        return a normalized direction vector from this point to another

        C++ signature :RDGeom::Point2D DirectionVector(RDGeom::Point2D {lvalue},RDGeom::Point2D)
        """
        ...
    def DotProduct(self, arg1: Point2D, arg2: Point2D) -> float:
        """
        Dot product with another point

        C++ signature :double DotProduct(RDGeom::Point2D {lvalue},RDGeom::Point2D)"""
        ...
    def Length(self, arg1: Point2D) -> float:
        """
        Length of the vector

        C++ signature :double Length(RDGeom::Point2D {lvalue})"""
        ...
    def LengthSq(self, arg1: Point2D) -> float:
        """
        Square of the length

        C++ signature :double LengthSq(RDGeom::Point2D {lvalue})"""
        ...
    def Normalize(self, arg1: Point2D) -> None:
        """
        Normalize the vector (using L2 norm)

        C++ signature :void Normalize(RDGeom::Point2D {lvalue})"""
        ...
    def SignedAngleTo(self, arg1: Point2D, arg2: Point2D) -> float:
        """
        determines the signed angle between a vector to this point (between 0 and 2*PI)

        C++ signature :double SignedAngleTo(RDGeom::Point2D {lvalue},RDGeom::Point2D)"""
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDGeom) -> Any: ...
    @classmethod
    def __getitem__(cls, RDGeom, int) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __iadd__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __idiv__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __imul__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __isub__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __len__(cls, RDGeom) -> Any: ...
    @classmethod
    def __mul__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...
    @classmethod
    def __truediv__(cls, RDGeom, double) -> Any: ...

class Point3D(Boost.Python.instance):
    """
    A class to represent a three-dimensional point
    The x, y, and z coordinates can be read and written using either attributes
    (i.e. pt.x = 4) or indexing (i.e. pt[0] = 4).
    Default Constructor

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (float)arg2, (float)arg3, (float)arg4) -> None :

    C++ signature :void __init__(_object*,double,double,double)"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...
    x: Any
    y: Any
    z: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AngleTo(self, arg1: Point3D, arg2: Point3D) -> float:
        """
        determines the angle between a vector to this point (between 0 and PI)

        C++ signature :double AngleTo(RDGeom::Point3D {lvalue},RDGeom::Point3D)"""
        ...
    def CrossProduct(self, arg1: Point3D, arg2: Point3D) -> Point3D:
        """
        Get the cross product between two points

        C++ signature :RDGeom::Point3D CrossProduct(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
        ...
    def DirectionVector(self, arg1: Point3D, arg2: Point3D) -> Point3D:
        """
        return a normalized direction vector from this point to another

        C++ signature :RDGeom::Point3D DirectionVector(RDGeom::Point3D {lvalue},RDGeom::Point3D)
        """
        ...
    def Distance(self, arg1: Point3D, arg2: Point3D) -> float:
        """
        Distance from this point to another point

        C++ signature :double Distance(RDGeom::Point3D,RDGeom::Point3D)"""
        ...
    def DotProduct(self, arg1: Point3D, arg2: Point3D) -> float:
        """
        Dot product with another point

        C++ signature :double DotProduct(RDGeom::Point3D {lvalue},RDGeom::Point3D)"""
        ...
    def Length(self, arg1: Point3D) -> float:
        """
        Length of the vector

        C++ signature :double Length(RDGeom::Point3D {lvalue})"""
        ...
    def LengthSq(self, arg1: Point3D) -> float:
        """
        Square of the length

        C++ signature :double LengthSq(RDGeom::Point3D {lvalue})"""
        ...
    def Normalize(self, arg1: Point3D) -> None:
        """
        Normalize the vector (using L2 norm)

        C++ signature :void Normalize(RDGeom::Point3D {lvalue})"""
        ...
    def SignedAngleTo(self, arg1: Point3D, arg2: Point3D) -> float:
        """
        determines the signed angle between a vector to this point (between 0 and 2*PI)

        C++ signature :double SignedAngleTo(RDGeom::Point3D {lvalue},RDGeom::Point3D)"""
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDGeom) -> Any: ...
    @classmethod
    def __getitem__(cls, RDGeom, int) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __iadd__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __idiv__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __imul__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __isub__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __len__(cls, RDGeom) -> Any: ...
    @classmethod
    def __mul__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...
    @classmethod
    def __truediv__(cls, RDGeom, double) -> Any: ...

class PointND(Boost.Python.instance):
    """
    A class to represent an N-dimensional point

    C++ signature :void __init__(_object*,unsigned int)"""

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def AngleTo(self, arg1: PointND, arg2: PointND) -> float:
        """
        determines the angle between a vector to this point (between 0 and PI)

        C++ signature :double AngleTo(RDGeom::PointND {lvalue},RDGeom::PointND)"""
        ...
    def DirectionVector(self, arg1: PointND, arg2: PointND) -> PointND:
        """
        return a normalized direction vector from this point to another

        C++ signature :RDGeom::PointND DirectionVector(RDGeom::PointND {lvalue},RDGeom::PointND)
        """
        ...
    def Distance(self, arg1: Point3D, arg2: Point3D) -> float:
        """
        Distance from this point to another point

        C++ signature :double Distance(RDGeom::Point3D,RDGeom::Point3D)"""
        ...
    def DotProduct(self, arg1: PointND, arg2: PointND) -> float:
        """
        Dot product with another point

        C++ signature :double DotProduct(RDGeom::PointND {lvalue},RDGeom::PointND)"""
        ...
    def Length(self, arg1: PointND) -> float:
        """
        Length of the vector

        C++ signature :double Length(RDGeom::PointND {lvalue})"""
        ...
    def LengthSq(self, arg1: PointND) -> float:
        """
        Square of the length

        C++ signature :double LengthSq(RDGeom::PointND {lvalue})"""
        ...
    def Normalize(self, arg1: PointND) -> None:
        """
        Normalize the vector (using L2 norm)

        C++ signature :void Normalize(RDGeom::PointND {lvalue})"""
        ...
    @classmethod
    def __add__(cls, other) -> Any: ...
    @classmethod
    def __getinitargs__(cls, RDGeom) -> Any: ...
    @classmethod
    def __getitem__(cls, RDGeom, int) -> Any: ...
    @classmethod
    def __getstate__(cls, RDGeom) -> Any: ...
    @classmethod
    def __iadd__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __idiv__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __imul__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __isub__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __len__(cls, RDGeom) -> Any: ...
    @classmethod
    def __mul__(cls, RDGeom, double) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setitem__(cls, RDGeom, int, double) -> Any: ...
    @classmethod
    def __setstate__(cls, RDGeom, boost) -> Any: ...
    @classmethod
    def __sub__(cls, other) -> Any: ...
    @classmethod
    def __truediv__(cls, RDGeom, double) -> Any: ...

class UniformGrid3D_(Boost.Python.instance):
    """
    Class to represent a uniform three-dimensional
    cubic grid. Each grid point can store a poisitive integer value. For the sake
    of efficiency these value can either be binary, fit in 2, 4, 8 or 16 bits
    pickle constructor

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def CompareParams(self, arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> bool:
        """
        Compare the parameters between two grid object

        C++ signature :bool CompareParams(RDGeom::UniformGrid3D {lvalue},RDGeom::UniformGrid3D)
        """
        ...
    def GetGridIndex(
        self, arg1: UniformGrid3D_, arg2: int, arg3: int, arg4: int
    ) -> int:
        """
        Get the index to the grid point with the three integer indices provided

        C++ signature :int GetGridIndex(RDGeom::UniformGrid3D {lvalue},unsigned int,unsigned int,unsigned int)
        """
        ...
    def GetGridIndices(self, arg1: UniformGrid3D_, arg2: int) -> tuple:
        """
        Returns the integer indices of the grid index provided.

        C++ signature :boost::python::tuple GetGridIndices(RDGeom::UniformGrid3D,unsigned int)
        """
        ...
    def GetGridPointIndex(self, arg1: UniformGrid3D_, arg2: Point3D) -> int:
        """
        Get the index to the grid point closest to the specified point

        C++ signature :int GetGridPointIndex(RDGeom::UniformGrid3D {lvalue},RDGeom::Point3D)
        """
        ...
    def GetGridPointLoc(self, arg1: UniformGrid3D_, arg2: int) -> Point3D:
        """
        Get the location of the specified grid point

        C++ signature :RDGeom::Point3D GetGridPointLoc(RDGeom::UniformGrid3D {lvalue},unsigned int)
        """
        ...
    def GetNumX(self, arg1: UniformGrid3D_) -> int:
        """
        Get the number of grid points along x-axis

        C++ signature :unsigned int GetNumX(RDGeom::UniformGrid3D {lvalue})"""
        ...
    def GetNumY(self, arg1: UniformGrid3D_) -> int:
        """
        Get the number of grid points along y-axis

        C++ signature :unsigned int GetNumY(RDGeom::UniformGrid3D {lvalue})"""
        ...
    def GetNumZ(self, arg1: UniformGrid3D_) -> int:
        """
        Get the number of grid points along z-axis

        C++ signature :unsigned int GetNumZ(RDGeom::UniformGrid3D {lvalue})"""
        ...
    def GetOccupancyVect(self, arg1: UniformGrid3D_) -> DiscreteValueVect:
        """
        Get the occupancy vector for the grid

        C++ signature :RDKit::DiscreteValueVect const* GetOccupancyVect(RDGeom::UniformGrid3D {lvalue})
        """
        ...
    def GetOffset(self, arg1: UniformGrid3D_) -> Point3D:
        """
        Get the location of the center of the grid

        C++ signature :RDGeom::Point3D GetOffset(RDGeom::UniformGrid3D {lvalue})"""
        ...
    def GetSize(self, arg1: UniformGrid3D_) -> int:
        """
        Get the size of the grid (number of grid points)

        C++ signature :unsigned int GetSize(RDGeom::UniformGrid3D {lvalue})"""
        ...
    def GetSpacing(self, arg1: UniformGrid3D_) -> float:
        """
        Get the grid spacing

        C++ signature :double GetSpacing(RDGeom::UniformGrid3D {lvalue})"""
        ...
    def GetVal(self, arg1: UniformGrid3D_, arg2: int) -> int:
        """
        Get the value at the specified grid point

        C++ signature :int GetVal(RDGeom::UniformGrid3D,unsigned int)"""
        ...
    def GetValPoint(self, arg1: UniformGrid3D_, arg2: Point3D) -> int:
        """
        Get the value at the closest grid point

        C++ signature :int GetValPoint(RDGeom::UniformGrid3D,RDGeom::Point3D)"""
        ...
    def SetSphereOccupancy(
        self: UniformGrid3D_,
        center: Point3D,
        radius: float,
        stepSize: float,
        maxLayers: int = -1,
        ignoreOutOfBound: bool = True,
    ) -> None:
        """
        Set the occupancy on the grid for a sphere or specified radiusand multiple layers around this sphere, with decreasing values of

        occupancy

        C++ signature :void SetSphereOccupancy(RDGeom::UniformGrid3D {lvalue},RDGeom::Point3D,double,double [,int=-1 [,bool=True]])
        """
        ...
    def SetVal(self, arg1: UniformGrid3D_, arg2: int, arg3: int) -> None:
        """
        Set the value at the specified grid point

        C++ signature :void SetVal(RDGeom::UniformGrid3D {lvalue},unsigned int,unsigned int)
        """
        ...
    def SetValPoint(self, arg1: UniformGrid3D_, arg2: Point3D, arg3: int) -> None:
        """
        Set the value at grid point closest to the specified point

        C++ signature :void SetValPoint(RDGeom::UniformGrid3D {lvalue},RDGeom::Point3D,unsigned int)
        """
        ...
    @classmethod
    def __getinitargs__(cls, RDGeom) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __iadd__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __iand__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __ior__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __isub__(cls, boost, RDGeom) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...

def ComputeDihedralAngle(
    self, arg1: Point3D, arg2: Point3D, arg3: Point3D, arg4: Point3D
) -> float:
    """
    calculates the dihedral angle determined by four Point3D objects

    C++ signature :double ComputeDihedralAngle(RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D)
    """
    ...

def ComputeGridCentroid(
    self, arg1: UniformGrid3D_, arg2: Point3D, arg3: float
) -> tuple:
    """
    Compute the grid point at the center of sphere around a Point3D

    C++ signature :boost::python::tuple ComputeGridCentroid(RDGeom::UniformGrid3D,RDGeom::Point3D,double)
    """
    ...

def ComputeSignedDihedralAngle(
    self, arg1: Point3D, arg2: Point3D, arg3: Point3D, arg4: Point3D
) -> float:
    """
    calculates the signed dihedral angle determined by four Point3D objects

    C++ signature :double ComputeSignedDihedralAngle(RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D,RDGeom::Point3D)
    """
    ...

def FindGridTerminalPoints(
    self, arg1: UniformGrid3D_, arg2: float, arg3: float
) -> tuple:
    """
    Find a grid’s terminal points (defined in the subshape algorithm).

    C++ signature :boost::python::tuple FindGridTerminalPoints(RDGeom::UniformGrid3D,double,double)
    """
    ...

def ProtrudeDistance(self, arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> float:
    """
    Compute the protrude distance between two grid objects

    C++ signature :double ProtrudeDistance(RDGeom::UniformGrid3D,RDGeom::UniformGrid3D)
    """
    ...

def TanimotoDistance(self, arg1: UniformGrid3D_, arg2: UniformGrid3D_) -> float:
    """
    Compute the tanimoto distance between two grid objects

    C++ signature :double TanimotoDistance(RDGeom::UniformGrid3D,RDGeom::UniformGrid3D)
    """
    ...

def TverskyIndex(
    self, arg1: UniformGrid3D_, arg2: UniformGrid3D_, arg3: float, arg4: float
) -> float:
    """
    Compute the tversky index between two grid objects

    C++ signature :double TverskyIndex(RDGeom::UniformGrid3D,RDGeom::UniformGrid3D,double,double)
    """
    ...

def UniformGrid3D(
    self,
    dimX: float,
    dimY: float,
    dimZ: float,
    spacing: float = 0.5,
    valType: DiscreteValueType = DiscreteValueType.TWOBITVALUE,
    offSet: Point3D = None,
) -> UniformGrid3D_:
    """
    Faking the constructor

    C++ signature :RDGeom::UniformGrid3D* UniformGrid3D(double,double,double [,double=0.5 [,RDKit::DiscreteValueVect::DiscreteValueType=rdkit.DataStructs.cDataStructs.DiscreteValueType.TWOBITVALUE [,RDGeom::Point3D const*=None]]])
    """
    ...

def WriteGridToFile(self, arg1: UniformGrid3D_, arg2: str) -> None:
    """
    Write the grid to a grid file

    C++ signature :void WriteGridToFile(RDGeom::UniformGrid3D,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...
