"""
rdkit.Chem.FeatMaps.FeatMapPoint module¶
"""
from _typeshed import Incomplete
from rdkit.Chem import ChemicalFeatures as ChemicalFeatures
from rdkit.Chem import rdChemicalFeatures

class FeatMapPoint(rdChemicalFeatures.FreeChemicalFeature):
    """
    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    __init__( (object)arg1) -> None :Default Constructor

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)family, (str)type, (Point3D)loc [, (int)id=-1]) -> None :Constructor with family, type and location specified

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point3D [,int=-1])

    __init__( (object)arg1, (str)family, (Point3D)loc) -> None :constructor with family and location specified, empty type and id

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDGeom::Point3D)
    """

    featDirs: None = ...
    weight: float = ...

    def GetDirMatch(self, other, useBest=True):
        """
        >>> from rdkit import Geometry
        >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
        >>> fmp = FeatMapPoint()
        >>> fmp.initFromFeat(sfeat)
        >>> fmp.GetDirMatch(sfeat)
        1.0

        >>> sfeat.featDirs=[Geometry.Point3D(0,0,1),Geometry.Point3D(0,0,-1)]
        >>> fmp.featDirs=[Geometry.Point3D(0,0,1),Geometry.Point3D(1,0,0)]
        >>> fmp.GetDirMatch(sfeat)
        1.0
        >>> fmp.GetDirMatch(sfeat,useBest=True)
        1.0
        >>> fmp.GetDirMatch(sfeat,useBest=False)
        0.0

        >>> sfeat.featDirs=[Geometry.Point3D(0,0,1)]
        >>> fmp.GetDirMatch(sfeat,useBest=False)
        0.5

        >>> sfeat.featDirs=[Geometry.Point3D(0,0,1)]
        >>> fmp.featDirs=[Geometry.Point3D(0,0,-1)]
        >>> fmp.GetDirMatch(sfeat)
        -1.0
        >>> fmp.GetDirMatch(sfeat,useBest=False)
        -1.0"""
        ...
    def GetDist2(self, other):
        """
        >>> from rdkit import Geometry
        >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
        >>> fmp = FeatMapPoint()
        >>> fmp.initFromFeat(sfeat)
        >>> fmp.GetDist2(sfeat)
        0.0
        >>> sfeat.SetPos(Geometry.Point3D(2,0,0))
        >>> fmp.GetDist2(sfeat)
        4.0"""
        ...
    featDirs: None = ...

    def __init__(self, *args, **kwargs) -> None: ...
    def initFromFeat(self, feat):
        """
        >>> from rdkit import Geometry
        >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
        >>> fmp = FeatMapPoint()
        >>> fmp.initFromFeat(sfeat)
        >>> fmp.GetFamily()==sfeat.GetFamily()
        True
        >>> fmp.GetType()==sfeat.GetType()
        True
        >>> list(fmp.GetPos())
        [0.0, 0.0, 0.0]
        >>> fmp.featDirs == []
        True

        >>> sfeat.featDirs = [Geometry.Point3D(1.0,0,0)]
        >>> fmp.initFromFeat(sfeat)
        >>> len(fmp.featDirs)
        1"""
        ...
    def GetDist2(self, other):
        """
        >>> from rdkit import Geometry
        >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
        >>> fmp = FeatMapPoint()
        >>> fmp.initFromFeat(sfeat)
        >>> fmp.GetDist2(sfeat)
        0.0
        >>> sfeat.SetPos(Geometry.Point3D(2,0,0))
        >>> fmp.GetDist2(sfeat)
        4.0"""
        ...
    def GetDirMatch(self, other, useBest=True):
        """
        >>> from rdkit import Geometry
        >>> sfeat = ChemicalFeatures.FreeChemicalFeature('Aromatic','Foo',Geometry.Point3D(0,0,0))
        >>> fmp = FeatMapPoint()
        >>> fmp.initFromFeat(sfeat)
        >>> fmp.GetDirMatch(sfeat)
        1.0

        >>> sfeat.featDirs=[Geometry.Point3D(0,0,1),Geometry.Point3D(0,0,-1)]
        >>> fmp.featDirs=[Geometry.Point3D(0,0,1),Geometry.Point3D(1,0,0)]
        >>> fmp.GetDirMatch(sfeat)
        1.0
        >>> fmp.GetDirMatch(sfeat,useBest=True)
        1.0
        >>> fmp.GetDirMatch(sfeat,useBest=False)
        0.0

        >>> sfeat.featDirs=[Geometry.Point3D(0,0,1)]
        >>> fmp.GetDirMatch(sfeat,useBest=False)
        0.5

        >>> sfeat.featDirs=[Geometry.Point3D(0,0,1)]
        >>> fmp.featDirs=[Geometry.Point3D(0,0,-1)]
        >>> fmp.GetDirMatch(sfeat)
        -1.0
        >>> fmp.GetDirMatch(sfeat,useBest=False)
        -1.0"""
        ...
