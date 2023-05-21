"""
rdkit.ML.Cluster.ClusterUtils module¶
utility functions for clustering
"""

def FindClusterCentroidFromDists(self, cluster, dists):
    """
    find the point in a cluster which has the smallest summedEuclidean distance to all others

    Arguments

    cluster: the cluster to work with
    dists: the distance matrix to use for the points

    Returns

    the index of the centroid point"""
    ...

def GetNodeList(self, cluster):
    """
    returns an ordered list of all nodes below cluster
    the ordering is done using the lengths of the child nodes

    Arguments

    cluster: the cluster in question

    Returns

    a list of the leaves below this cluster"""
    ...

def GetNodesDownToCentroids(self, cluster, above=1):
    """
    returns an ordered list of all nodes below cluster"""
    ...

def FindClusterCentroidFromDists(self, cluster, dists):
    """
    find the point in a cluster which has the smallest summedEuclidean distance to all others

    Arguments

    cluster: the cluster to work with
    dists: the distance matrix to use for the points

    Returns

    the index of the centroid point"""
    ...

def SplitIntoNClusters(self, cluster, n, breadthFirst=True):
    """
    splits a cluster tree into a set of branches
    Arguments

    cluster: the root of the cluster tree
    n: the number of clusters to include in the split
    breadthFirst: toggles breadth first (vs depth first) cleavage
    of the cluster tree.

    Returns

    a list of sub clusters"""
    ...
