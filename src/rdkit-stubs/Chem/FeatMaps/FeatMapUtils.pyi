"""
rdkit.Chem.FeatMaps.FeatMapUtils module¶
"""
from rdkit.Chem.FeatMaps import FeatMaps as FeatMaps

class DirMergeMode(object):
    NoMerge: int = ...
    Sum: int = ...

    def valid(self, dirMergeMode):
        """
        Check that dirMergeMode is valid"""
        ...

class MergeMethod(object):
    WeightedAverage: int = ...
    Average: int = ...
    UseLarger: int = ...
    WeightedAverage: int = ...

    def valid(self, mergeMethod):
        """
        Check that mergeMethod is valid"""
        ...

class MergeMetric(object):
    Distance: int = ...
    NoMerge: int = ...
    Distance: int = ...
    Overlap: int = ...

    def valid(self, mergeMetric):
        """
        Check that mergeMetric is valid"""
        ...

def CombineFeatMaps(self, fm1, fm2, mergeMetric=0, mergeTol=1.5, dirMergeMode=0):
    """
    the parameters will be taken from fm1"""
    ...

class DirMergeMode(object):
    NoMerge: int = ...
    Sum: int = ...

    def valid(self, dirMergeMode):
        """
        Check that dirMergeMode is valid"""
        ...

def GetFeatFeatDistMatrix(self, fm, mergeMetric, mergeTol, dirMergeMode, compatFunc):
    """
    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging."""
    ...

def MergeFeatPoints(
    self,
    fm,
    mergeMetric=0,
    mergeTol=1.5,
    dirMergeMode=0,
    mergeMethod=0,
    compatFunc=familiesMatch,
):
    """
    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging.
    returns whether or not any points were actually merged"""
    ...

def familiesMatch(self, f1, f2): ...
def feq(self, v1, v2, tol=0.0001): ...
def MergeFeatPoints(
    self,
    fm,
    mergeMetric=0,
    mergeTol=1.5,
    dirMergeMode=0,
    mergeMethod=0,
    compatFunc=familiesMatch,
):
    """
    NOTE that mergeTol is a max value for merging when using distance-based
    merging and a min value when using score-based merging.
    returns whether or not any points were actually merged"""
    ...

def CombineFeatMaps(self, fm1, fm2, mergeMetric=0, mergeTol=1.5, dirMergeMode=0):
    """
    the parameters will be taken from fm1"""
    ...
