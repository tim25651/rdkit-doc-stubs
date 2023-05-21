"""
rdkit.ML.InfoTheory.BitClusterer module¶
"""
from rdkit import DataStructs as DataStructs

class BitClusterer(object):
    """
    Class to cluster a set of bits based on their correllation
    The correlation matrix is first built using by reading the fingerprints
    from a database or a list of fingerprints"""

    def __init__(self, idList, nCluster, type=...) -> None: ...
    def ClusterBits(self, corrMat): ...
    def SetClusters(self, clusters): ...
    def GetClusters(self): ...
    def MapToClusterFP(self, fp):
        """
        Map the fingerprint to a smaller sized (= number of clusters) fingerprint
        Each cluster get a bit in the new fingerprint and is turned on if any of the bits in
        the cluster are turned on in the original fingerprint"""
        ...
    def MapToClusterScores(self, fp):
        """
        Map the fingerprint to a real valued vector of score based on the bit clusters
        The dimension of the vector is same as the number of clusters. Each value in the
        vector corresponds to the number of bits in the corresponding cluster
        that are turned on in the fingerprint

        ARGUMENTS:
        fp : the fingerprint"""
        ...
    def SetClusters(self, clusters): ...
    def MapToClusterFP(self, fp):
        """
        Map the fingerprint to a smaller sized (= number of clusters) fingerprint
        Each cluster get a bit in the new fingerprint and is turned on if any of the bits in
        the cluster are turned on in the original fingerprint"""
        ...
