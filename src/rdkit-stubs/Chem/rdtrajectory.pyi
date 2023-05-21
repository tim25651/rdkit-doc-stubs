"""
rdkit.Chem.rdtrajectory module¶
Module containing Trajectory and Snapshot objects
"""
from typing import Any

import Boost.Python
from rdkit.Chem.rdchem import Mol
from rdkit.Geometry.rdGeometry import Point2D, Point3D

class Snapshot(Boost.Python.instance):
    """
    A class which allows storing coordinates from a trajectory

    Constructor;coordList: list of floats containing the coordinates for this Snapshot;
    energy:    the energy for this Snapshot.

    C++ signature :void* __init__(boost::python::api::object,boost::python::list {lvalue} [,double=0.0])

    __init__( (AtomPairsParameters)arg1, (Snapshot)other) -> object :Copy constructor

    C++ signature :void* __init__(boost::python::api::object,RDKit::Snapshot*)"""

    @classmethod
    def __init__(cls, boost, RDKit) -> Any: ...
    def GetEnergy(self: Snapshot) -> float:
        """
        returns the energy for this Snapshot

        C++ signature :double GetEnergy(RDKit::Snapshot {lvalue})"""
        ...
    def GetPoint2D(self: Snapshot, pointNum: int) -> Point2D:
        """
        return the coordinates at pointNum as a Point2D object; requires the Trajectory dimension to be == 2

        C++ signature :RDGeom::Point2D GetPoint2D(RDKit::Snapshot {lvalue},unsigned int)
        """
        ...
    def GetPoint3D(self: Snapshot, pointNum: int) -> Point3D:
        """
        return the coordinates at pointNum as a Point3D object; requires the Trajectory dimension to be >= 2

        C++ signature :RDGeom::Point3D GetPoint3D(RDKit::Snapshot {lvalue},unsigned int)
        """
        ...
    def SetEnergy(self: Snapshot, energy: float) -> None:
        """
        sets the energy for this Snapshot

        C++ signature :void SetEnergy(RDKit::Snapshot {lvalue},double)"""
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class Trajectory(Boost.Python.instance):
    """
    A class which allows storing Snapshots from a trajectory

    Constructor;dimension:    dimensionality of this Trajectory’s coordinate tuples;
    numPoints:    number of coordinate tuples associated to each Snapshot;
    snapshotList: list of Snapshot objects used to initialize the Trajectory (optional; defaults to []).

    C++ signature :void* __init__(boost::python::api::object,unsigned int,unsigned int [,boost::python::list=[]])

    __init__( (AtomPairsParameters)arg1, (Trajectory)other) -> object :Copy constructor

    C++ signature :void* __init__(boost::python::api::object,RDKit::Trajectory*)"""

    @classmethod
    def __init__(cls, boost, RDKit) -> Any: ...
    def AddConformersToMol(
        self: Trajectory, mol: Mol, from_: int = -1, to: int = -1
    ) -> int:
        """
        adds conformations from the Trajectory to mol
        from is the first Snapshot that will be added as a Conformer; defaults to -1 (first available)
        to is the last Snapshot that will be added as a Conformer; defaults to -1 (all)

        C++ signature :unsigned int AddConformersToMol(RDKit::Trajectory {lvalue},RDKit::ROMol {lvalue} [,int=-1 [,int=-1]])
        """
        ...
    def AddSnapshot(self: Trajectory, s: Snapshot) -> int:
        """
        appends Snapshot s to this Trajectory; returns the zero-based index position of the added snapshot

        C++ signature :unsigned int AddSnapshot(RDKit::Trajectory {lvalue},RDKit::Snapshot)
        """
        ...
    def Clear(self: Trajectory) -> None:
        """
        removes all Snapshots from the Trajectory

        C++ signature :void Clear(RDKit::Trajectory {lvalue})"""
        ...
    def Dimension(self: Trajectory) -> int:
        """
        returns the dimensionality of this Trajectory’s coordinate tuples

        C++ signature :unsigned int Dimension(RDKit::Trajectory {lvalue})"""
        ...
    def GetSnapshot(self: Trajectory, snapshotNum: int) -> Snapshot:
        """
        returns the Snapshot snapshotNum, where the latter is the zero-based index of the retrieved Snapshot

        C++ signature :RDKit::Snapshot* GetSnapshot(RDKit::Trajectory*,unsigned int)"""
        ...
    def InsertSnapshot(self: Trajectory, snapshotNum: int, s: Snapshot) -> int:
        """
        inserts Snapshot s into the Trajectory at the position snapshotNum, where the latter is the zero-based index of the Trajectory’s Snapshot before which the Snapshot s will be inserted; returns the zero-based index position of the inserted snapshot

        C++ signature :unsigned int InsertSnapshot(RDKit::Trajectory {lvalue},unsigned int,RDKit::Snapshot)
        """
        ...
    def NumPoints(self: Trajectory) -> int:
        """
        returns the number of coordinate tuples associated to each Snapshot

        C++ signature :unsigned int NumPoints(RDKit::Trajectory {lvalue})"""
        ...
    def RemoveSnapshot(self: Trajectory, snapshotNum: int) -> int:
        """
        removes Snapshot snapshotNum from the Trajectory, where snapshotNum is the zero-based index of Snapshot to be removed

        C++ signature :unsigned int RemoveSnapshot(RDKit::Trajectory {lvalue},unsigned int)
        """
        ...
    @classmethod
    def __len__(cls, RDKit) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def ReadAmberTrajectory(self, fName: str, traj: Trajectory) -> int:
    """
    reads coordinates from an AMBER trajectory file into the Trajectory object; returns the number of Snapshot objects read in

    C++ signature :unsigned int ReadAmberTrajectory(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::Trajectory {lvalue})
    """
    ...

def ReadGromosTrajectory(self, fName: str, traj: Trajectory) -> int:
    """
    reads coordinates from a GROMOS trajectory file into the Trajectory object; returns the number of Snapshot objects read in

    C++ signature :unsigned int ReadGromosTrajectory(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,RDKit::Trajectory {lvalue})
    """
    ...
