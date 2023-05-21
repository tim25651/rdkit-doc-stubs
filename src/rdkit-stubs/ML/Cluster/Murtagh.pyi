"""
rdkit.ML.Cluster.Murtagh module¶
Interface to the C++ Murtagh hierarchic clustering code
"""
from _typeshed import Incomplete
from rdkit.ML.Cluster import Clusters as Clusters
from rdkit.ML.Cluster.Clustering import MurtaghCluster as MurtaghCluster
from rdkit.ML.Cluster.Clustering import MurtaghDistCluster as MurtaghDistCluster

WARDS: int
SLINK: int
CLINK: int
UPGMA: int
MCQUITTY: int
GOWER: int
CENTROID: int
methods: Incomplete

def ClusterData(self, data, nPts, method, isDistData=0):
    """
    clusters the data points passed in and returns the cluster tree
    Arguments

    data: a list of lists (or array, or whatever) with the input
    data (see discussion of _isDistData_ argument for the exception)
    nPts: the number of points to be used

    method: determines which clustering algorithm should be used.The defined constants for these are:
    ‘WARDS, SLINK, CLINK, UPGMA’

    isDistData: set this toggle when the data passed in is adistance matrix.  The distance matrix should be stored
    symmetrically so that _LookupDist (above) can retrieve
    the results:

    for i<j: d_ij = dists[j*(j-1)//2 + i]

    Returns

    a single entry list with the cluster tree"""
    ...
