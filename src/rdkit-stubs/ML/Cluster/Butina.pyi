"""
rdkit.ML.Cluster.Butina module¶
Implementation of the clustering algorithm published in:
Butina JCICS 39 747-750 (1999)
"""
from _typeshed import Incomplete
from rdkit import RDLogger as RDLogger
from rdkit.ML.KNN import DistFunctions

def ClusterData(
    self,
    data,
    nPts,
    distThresh,
    isDistData=False,
    distFunc=DistFunctions.EuclideanDist,
    reordering=False,
):
    """
    clusters the data points passed in and returns the list of clusters
    Arguments

    data: a list of items with the input data
    (see discussion of _isDistData_ argument for the exception)
    nPts: the number of points to be used
    distThresh: elements within this range of each other are considered
    to be neighbors

    isDistData: set this toggle when the data passed in is adistance matrix.  The distance matrix should be stored
    symmetrically. An example of how to do this:

    dists = []
    for i in range(nPts):

    for j in range(i):dists.append( distfunc(i,j) )

    distFunc: a function to calculate distances between points.Receives 2 points as arguments, should return a float

    reordering: if this toggle is set, the number of neighbors is updatedfor the unassigned molecules after a new cluster is created such
    that always the molecule with the largest number of unassigned
    neighbors is selected as the next cluster center.

    Returns

    a tuple of tuples containing information about the clusters:
    ( (cluster1_elem1, cluster1_elem2, …),(cluster2_elem1, cluster2_elem2, …),
    …

    )
    The first element for each cluster is its centroid."""
    ...

logger: Incomplete

def EuclideanDist(self, pi, pj): ...
def ClusterData(
    self,
    data,
    nPts,
    distThresh,
    isDistData=False,
    distFunc=DistFunctions.EuclideanDist,
    reordering=False,
):
    """
    clusters the data points passed in and returns the list of clusters
    Arguments

    data: a list of items with the input data
    (see discussion of _isDistData_ argument for the exception)
    nPts: the number of points to be used
    distThresh: elements within this range of each other are considered
    to be neighbors

    isDistData: set this toggle when the data passed in is adistance matrix.  The distance matrix should be stored
    symmetrically. An example of how to do this:

    dists = []
    for i in range(nPts):

    for j in range(i):dists.append( distfunc(i,j) )

    distFunc: a function to calculate distances between points.Receives 2 points as arguments, should return a float

    reordering: if this toggle is set, the number of neighbors is updatedfor the unassigned molecules after a new cluster is created such
    that always the molecule with the largest number of unassigned
    neighbors is selected as the next cluster center.

    Returns

    a tuple of tuples containing information about the clusters:
    ( (cluster1_elem1, cluster1_elem2, …),(cluster2_elem1, cluster2_elem2, …),
    …

    )
    The first element for each cluster is its centroid."""
    ...
