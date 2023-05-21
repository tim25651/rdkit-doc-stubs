"""
rdkit.Chem.Fingerprints.ClusterMols module¶

utility functionality for clustering molecules using fingerprintsincludes a command line app for clustering

Sample Usage:python ClusterMols.py  -d data.gdb -t daylight_sig     –idName=”CAS_TF” -o clust1.pkl     –actTable=”dop_test” –actName=”moa_quant”
"""
from rdkit import DataStructs as DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols as FingerprintMols
from rdkit.Chem.Fingerprints import MolSimilarity as MolSimilarity
from rdkit.ML.Cluster import Murtagh as Murtagh

def ClusterFromDetails(self, details):
    """
    Returns the cluster tree"""
    ...

def ClusterPoints(
    self,
    data,
    metric,
    algorithmId,
    haveLabels=False,
    haveActs=True,
    returnDistances=False,
): ...

message = FingerprintMols.message
error = FingerprintMols.error

def GetDistanceMatrix(self, data, metric, isSimilarity=1):
    """
    data should be a list of tuples with fingerprints in position 1
    (the rest of the elements of the tuple are not important)

    Returns the symmetric distance matrix
    (see ML.Cluster.Resemblance for layout documentation)"""
    ...

def ClusterPoints(
    self,
    data,
    metric,
    algorithmId,
    haveLabels=False,
    haveActs=True,
    returnDistances=False,
): ...
def ClusterFromDetails(self, details):
    """
    Returns the cluster tree"""
    ...
