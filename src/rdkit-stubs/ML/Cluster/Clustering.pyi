"""
rdkit.ML.Cluster.Clustering module¶
"""
from typing import Any

from rdkit.Chem.rdMolDescriptors import AtomPairsParameters

def MurtaghCluster(
    self, data: AtomPairsParameters, nPts: int, sz: int, option: int
) -> object:
    """
    TODO: provide docstring

    C++ signature :_object* MurtaghCluster(boost::python::api::object,int,int,int)"""
    ...

def MurtaghDistCluster(
    self, data: AtomPairsParameters, nPts: int, option: int
) -> object:
    """
    TODO: provide docstring

    C++ signature :_object* MurtaghDistCluster(boost::python::api::object,int,int)"""
    ...
