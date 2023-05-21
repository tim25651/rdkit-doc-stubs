"""
rdkit.SimDivFilters.rdSimDivPickers module¶
Module containing the diversity and similarity pickers
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters
from rdkit.rdBase import _vecti, _vectSt6vectorIiSaIiEE

CENTROID: ClusterMethod
CLINK: ClusterMethod
GOWER: ClusterMethod
MCQUITTY: ClusterMethod
SLINK: ClusterMethod
UPGMA: ClusterMethod
WARD: ClusterMethod

class ClusterMethod(Boost.Python.enum):
    CENTROID: ClusterMethod = ...
    CLINK: ClusterMethod = ...
    GOWER: ClusterMethod = ...
    MCQUITTY: ClusterMethod = ...
    SLINK: ClusterMethod = ...
    UPGMA: ClusterMethod = ...
    WARD: ClusterMethod = ...
    names: dict[str, ClusterMethod] = ...
    values: dict[int, ClusterMethod] = ...
    __slots__: ClassVar[tuple] = ...

class HierarchicalClusterPicker(Boost.Python.instance):
    """
    A class for diversity picking of items using Hierarchical Clustering

    C++ signature :void __init__(_object*,RDPickers::HierarchicalClusterPicker::ClusterMethod)
    """

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def Cluster(
        self,
        arg1: HierarchicalClusterPicker,
        arg2: AtomPairsParameters,
        arg3: int,
        arg4: int,
    ) -> _vectSt6vectorIiSaIiEE:
        """
        Return a list of clusters of item from the pool using hierarchical clustering

        ARGUMENTS:
        distMat: 1D distance matrix (only the lower triangle elements)
        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool

        C++ signature :std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > Cluster(RDPickers::HierarchicalClusterPicker*,boost::python::api::object {lvalue},int,int)
        """
        ...
    def Pick(
        self,
        arg1: HierarchicalClusterPicker,
        arg2: AtomPairsParameters,
        arg3: int,
        arg4: int,
    ) -> _vecti:
        """
        Pick a diverse subset of items from a pool of items using hierarchical clustering

        ARGUMENTS:
        distMat: 1D distance matrix (only the lower triangle elements)
        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool

        C++ signature :std::vector<int, std::allocator<int> > Pick(RDPickers::HierarchicalClusterPicker*,boost::python::api::object {lvalue},int,int)
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class LeaderPicker(Boost.Python.instance):
    """
    A class for diversity picking of items using Roger Sayle’s Leader algorithm (analogous to sphere exclusion). The algorithm is currently unpublished, but a description is available in this presentation from the 2019 RDKit UGM: https://github.com/rdkit/UGM_2019/raw/master/Presentations/Sayle_Clustering.pdf

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def LazyBitVectorPick(
        self: LeaderPicker,
        objects: AtomPairsParameters,
        poolSize: int,
        threshold: float,
        pickSize: int = 0,
        firstPicks: AtomPairsParameters = (),
        numThreads: int = 1,
    ) -> _vecti:
        """
        Pick a subset of items from a collection of bit vectors using Tanimoto distance. The threshold value is a distance (i.e. 1-similarity). Note that the numThreads argument is currently ignored.

        C++ signature :std::vector<int, std::allocator<int> > LazyBitVectorPick(RDPickers::LeaderPicker*,boost::python::api::object,int,double [,int=0 [,boost::python::api::object=() [,int=1]]])
        """
        ...
    def LazyPick(
        self: LeaderPicker,
        distFunc: AtomPairsParameters,
        poolSize: int,
        threshold: float,
        pickSize: int = 0,
        firstPicks: AtomPairsParameters = (),
        numThreads: int = 1,
    ) -> _vecti:
        """
        Pick a subset of items from a pool of items using the user-provided function to determine distances. Note that the numThreads argument is currently ignored.

        C++ signature :std::vector<int, std::allocator<int> > LazyPick(RDPickers::LeaderPicker*,boost::python::api::object,int,double [,int=0 [,boost::python::api::object=() [,int=1]]])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...

class MaxMinPicker(Boost.Python.instance):
    """
    A class for diversity picking of items using the MaxMin Algorithm

    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def LazyBitVectorPick(
        self: MaxMinPicker,
        objects: AtomPairsParameters,
        poolSize: int,
        pickSize: int,
        firstPicks: AtomPairsParameters = (),
        seed: int = -1,
        useCache: AtomPairsParameters = None,
    ) -> _vecti:
        """
        Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

        vectors: a sequence of the bit vectors that should be picked from.
        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool
        firstPicks: (optional) the first items to be picked (seeds the list)
        seed: (optional) seed for the random number generator
        useCache: IGNORED.

        C++ signature :std::vector<int, std::allocator<int> > LazyBitVectorPick(RDPickers::MaxMinPicker*,boost::python::api::object,int,int [,boost::python::api::object=() [,int=-1 [,boost::python::api::object=None]]])
        """
        ...
    def LazyBitVectorPickWithThreshold(
        self: MaxMinPicker,
        objects: AtomPairsParameters,
        poolSize: int,
        pickSize: int,
        threshold: float,
        firstPicks: AtomPairsParameters = (),
        seed: int = -1,
    ) -> tuple:
        """
        Pick a subset of items from a pool of bit vectors using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

        vectors: a sequence of the bit vectors that should be picked from.
        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool
        threshold: stop picking when the distance goes below this value
        firstPicks: (optional) the first items to be picked (seeds the list)
        seed: (optional) seed for the random number generator

        C++ signature :boost::python::tuple LazyBitVectorPickWithThreshold(RDPickers::MaxMinPicker*,boost::python::api::object,int,int,double [,boost::python::api::object=() [,int=-1]])
        """
        ...
    def LazyPick(
        self: MaxMinPicker,
        distFunc: AtomPairsParameters,
        poolSize: int,
        pickSize: int,
        firstPicks: AtomPairsParameters = (),
        seed: int = -1,
        useCache: AtomPairsParameters = None,
    ) -> _vecti:
        """
        Pick a subset of items from a pool of items using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

        distFunc: a function that should take two indices and return thedistance between those two points.
        NOTE: the implementation caches distance values, so the
        client code does not need to do so; indeed, it should not.

        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool
        firstPicks: (optional) the first items to be picked (seeds the list)
        seed: (optional) seed for the random number generator
        useCache: IGNORED

        C++ signature :std::vector<int, std::allocator<int> > LazyPick(RDPickers::MaxMinPicker*,boost::python::api::object,int,int [,boost::python::api::object=() [,int=-1 [,boost::python::api::object=None]]])
        """
        ...
    def LazyPickWithThreshold(
        self: MaxMinPicker,
        distFunc: AtomPairsParameters,
        poolSize: int,
        pickSize: int,
        threshold: float,
        firstPicks: AtomPairsParameters = (),
        seed: int = -1,
    ) -> tuple:
        """
        Pick a subset of items from a pool of items using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604
        ARGUMENTS:

        distFunc: a function that should take two indices and return thedistance between those two points.
        NOTE: the implementation caches distance values, so the
        client code does not need to do so; indeed, it should not.

        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool
        threshold: stop picking when the distance goes below this value
        firstPicks: (optional) the first items to be picked (seeds the list)
        seed: (optional) seed for the random number generator

        C++ signature :boost::python::tuple LazyPickWithThreshold(RDPickers::MaxMinPicker*,boost::python::api::object,int,int,double [,boost::python::api::object=() [,int=-1]])
        """
        ...
    def Pick(
        self: MaxMinPicker,
        distMat: AtomPairsParameters,
        poolSize: int,
        pickSize: int,
        firstPicks: AtomPairsParameters = (),
        seed: int = -1,
    ) -> _vecti:
        """
        Pick a subset of items from a pool of items using the MaxMin Algorithm
        Ashton, M. et. al., Quant. Struct.-Act. Relat., 21 (2002), 598-604

        ARGUMENTS:
        distMat: 1D distance matrix (only the lower triangle elements)
        poolSize: number of items in the pool
        pickSize: number of items to pick from the pool
        firstPicks: (optional) the first items to be picked (seeds the list)
        seed: (optional) seed for the random number generator

        C++ signature :std::vector<int, std::allocator<int> > Pick(RDPickers::MaxMinPicker*,boost::python::api::object,int,int [,boost::python::api::object=() [,int=-1]])
        """
        ...
    @classmethod
    def __reduce__(cls) -> Any: ...
