"""
rdkit.Chem.Scaffolds.rdScaffoldNetwork module¶
Module containing functions for creating a Scaffold Network
"""
from typing import Any, ClassVar

import Boost.Python
from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

class EdgeType(Boost.Python.enum):
    Fragment: EdgeType = ...
    Generic: EdgeType = ...
    GenericBond: EdgeType = ...
    Initialize: EdgeType = ...
    RemoveAttachment: EdgeType = ...
    names: dict[str, EdgeType] = ...
    values: dict[int, EdgeType] = ...
    __slots__: ClassVar[tuple] = ...

class NetworkEdge(Boost.Python.instance):
    """
    A scaffold network edge
    Raises an exception
    This class cannot be instantiated from Python

    property beginIdx¶
    index of the begin node in node list

    property endIdx¶
    index of the end node in node list

    property type¶
    type of the edge"""

    ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @property
    def beginIdx(self) -> Any: ...
    @property
    def endIdx(self) -> Any: ...
    @property
    def type(self) -> Any: ...

class NetworkEdge_VECT(Boost.Python.instance):
    """
    C++ signature :void __init__(_object*)"""

    __instance_size__: ClassVar[int] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    def append(self, arg1: NetworkEdge_VECT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void append(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},boost::python::api::object)
        """
        ...
    def extend(self, arg1: NetworkEdge_VECT, arg2: AtomPairsParameters) -> None:
        """
        C++ signature :void extend(std::vector<RDKit::ScaffoldNetwork::NetworkEdge, std::allocator<RDKit::ScaffoldNetwork::NetworkEdge> > {lvalue},boost::python::api::object)
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

class ScaffoldNetwork(Boost.Python.instance):
    """
    A scaffold network

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (str)arg2) -> None :

    C++ signature :void __init__(_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)

    property counts¶
    the number of times each node was encountered while building the network.

    property edges¶
    the sequence of network edges

    property molCounts¶
    the number of moleclues each node was found in.

    property nodes¶
    the sequence of SMILES defining the nodes"""

    ...
    __getstate_manages_dict__: ClassVar[bool] = ...
    __instance_size__: ClassVar[int] = ...
    __safe_for_unpickling__: ClassVar[bool] = ...

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __getinitargs__(cls, RDKit) -> Any: ...
    @classmethod
    def __getstate__(cls, boost) -> Any: ...
    @classmethod
    def __reduce__(cls) -> Any: ...
    @classmethod
    def __setstate__(cls, state) -> Any: ...
    @property
    def counts(self) -> Any: ...
    @property
    def edges(self) -> Any: ...
    @property
    def molCounts(self) -> Any: ...
    @property
    def nodes(self) -> Any: ...

class ScaffoldNetworkParams(Boost.Python.instance):
    """
    Scaffold network parameters

    C++ signature :void __init__(_object*)

    __init__( (object)arg1, (_vectNSt7__cxx1112basic_stringIcSt11char_traitsIcESaIcEEE)bondBreakerSmartsList) -> None :Constructor taking a list of Reaction SMARTS for the fragmentation reactions

    C++ signature :void __init__(_object*,std::vector<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >, std::allocator<std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > > >)

    property collectMolCounts¶
    keep track of the number of molecules each scaffold was found in

    property flattenChirality¶
    remove chirality and bond stereo when flattening

    property flattenIsotopes¶
    remove isotopes when flattening

    property flattenKeepLargest¶
    keep only the largest fragment when doing flattening

    property includeGenericBondScaffolds¶
    include scaffolds with all bonds replaced by single bonds

    property includeGenericScaffolds¶
    include scaffolds with all atoms replaced by dummies

    property includeScaffoldsWithAttachments¶
    Include the version of the scaffold with attachment points

    property includeScaffoldsWithoutAttachments¶
    remove attachment points from scaffolds and include the result

    property keepOnlyFirstFragment¶
    keep only the first fragment from the bond breaking rule

    property pruneBeforeFragmenting¶
    Do a pruning/flattening step before starting fragmenting"""

    ...
    __instance_size__: ClassVar[int] = ...
    collectMolCounts: Any
    flattenChirality: Any
    flattenIsotopes: Any
    flattenKeepLargest: Any
    includeGenericBondScaffolds: Any
    includeGenericScaffolds: Any
    includeScaffoldsWithAttachments: Any
    includeScaffoldsWithoutAttachments: Any
    keepOnlyFirstFragment: Any
    pruneBeforeFragmenting: Any

    @classmethod
    def __init__(cls, *args, **kwargs) -> None: ...
    @classmethod
    def __reduce__(cls) -> Any: ...

def BRICSScaffoldParams(self) -> ScaffoldNetworkParams:
    """
    Returns parameters for generating scaffolds using BRICS fragmentation rules

    C++ signature :RDKit::ScaffoldNetwork::ScaffoldNetworkParams* BRICSScaffoldParams()
    """
    ...

def CreateScaffoldNetwork(
    self, mols: AtomPairsParameters, params: ScaffoldNetworkParams
) -> ScaffoldNetwork:
    """
    create (and return) a new network from a sequence of molecules

    C++ signature :RDKit::ScaffoldNetwork::ScaffoldNetwork* CreateScaffoldNetwork(boost::python::api::object,RDKit::ScaffoldNetwork::ScaffoldNetworkParams)
    """
    ...

def UpdateScaffoldNetwork(
    self,
    mols: AtomPairsParameters,
    network: ScaffoldNetwork,
    params: ScaffoldNetworkParams,
) -> None:
    """
    update an existing network by adding molecules

    C++ signature :void UpdateScaffoldNetwork(boost::python::api::object,RDKit::ScaffoldNetwork::ScaffoldNetwork {lvalue},RDKit::ScaffoldNetwork::ScaffoldNetworkParams)
    """
    ...
