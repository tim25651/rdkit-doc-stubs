"""
rdkit.Chem.Subshape.SubshapeAligner module¶
"""
from collections.abc import Generator

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import Geometry as Geometry
from rdkit import RDLogger as RDLogger
from rdkit.Chem.Subshape import SubshapeObjects as SubshapeObjects
from rdkit.Numerics import Alignment as Alignment

class SubshapeAligner(object):
    coarseGridToleranceMult: float = ...
    dirThresh: float = ...
    triangleRMSTol: float = ...
    distMetric: int = ...
    shapeDistTol: float = ...
    numFeatThresh: int = ...
    dirThresh: float = ...
    edgeTol: float = ...
    coarseGridToleranceMult: float = ...
    medGridToleranceMult: float = ...
    numFeatThresh: int = ...
    shapeDistTol: float = ...
    triangleRMSTol: float = ...

    def GetSubshapeAlignments(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        tgtConf=-1,
        queryConf=-1,
        pruneStats=None,
    ): ...
    def GetTriangleMatches(self, target, query):
        """
        this is a generator function returning the possible triangle
        matches between the two shapes"""
        ...
    def PruneMatchesUsingDirection(
        self, target, query, alignments, pruneStats=None
    ): ...
    def PruneMatchesUsingFeatures(self, target, query, alignments, pruneStats=None): ...
    def PruneMatchesUsingDirection(
        self, target, query, alignments, pruneStats=None
    ): ...
    def PruneMatchesUsingShape(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        alignments,
        tgtConf=-1,
        queryConf=-1,
        pruneStats=None,
    ): ...
    def GetSubshapeAlignments(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        tgtConf=-1,
        queryConf=-1,
        pruneStats=None,
    ): ...
    def __call__(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        tgtConf: int = ...,
        queryConf: int = ...,
        pruneStats: Incomplete | None = ...,
    ) -> Generator[Incomplete, None, None]: ...

logger: Incomplete

class SubshapeAlignment(object):
    transform: None = ...
    triangleSSD: None = ...
    targetTri: None = ...
    queryTri: None = ...
    alignedConfId: int = ...
    dirMatch: float = ...
    queryTri: None = ...
    shapeDist: float = ...
    targetTri: None = ...
    transform: None = ...
    triangleSSD: None = ...

class SubshapeDistanceMetric(object):
    PROTRUDE: int = ...
    TANIMOTO: int = ...
    PROTRUDE: int = ...

def ClusterAlignments(
    self, mol, alignments, builder, neighborTol=0.1, distMetric=1, tempConfId=1001
):
    """
    clusters a set of alignments and returns the cluster centroid"""
    ...

def GetShapeShapeDistance(self, s1, s2, distMetric):
    """
    returns the distance between two shapes according to the provided metric"""
    ...

def ClusterAlignments(
    self, mol, alignments, builder, neighborTol=0.1, distMetric=1, tempConfId=1001
):
    """
    clusters a set of alignments and returns the cluster centroid"""
    ...

def TransformMol(self, mol, tform, confId=-1, newConfId=100):
    """
    Applies the transformation to a molecule and sets it up with a single conformer"""
    ...

class SubshapeAligner(object):
    coarseGridToleranceMult: float = ...
    dirThresh: float = ...
    triangleRMSTol: float = ...
    distMetric: int = ...
    shapeDistTol: float = ...
    numFeatThresh: int = ...
    dirThresh: float = ...
    edgeTol: float = ...
    coarseGridToleranceMult: float = ...
    medGridToleranceMult: float = ...
    numFeatThresh: int = ...
    shapeDistTol: float = ...
    triangleRMSTol: float = ...

    def GetSubshapeAlignments(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        tgtConf=-1,
        queryConf=-1,
        pruneStats=None,
    ): ...
    def GetTriangleMatches(self, target, query):
        """
        this is a generator function returning the possible triangle
        matches between the two shapes"""
        ...
    def PruneMatchesUsingDirection(
        self, target, query, alignments, pruneStats=None
    ): ...
    def PruneMatchesUsingFeatures(self, target, query, alignments, pruneStats=None): ...
    def PruneMatchesUsingDirection(
        self, target, query, alignments, pruneStats=None
    ): ...
    def PruneMatchesUsingShape(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        alignments,
        tgtConf=-1,
        queryConf=-1,
        pruneStats=None,
    ): ...
    def GetSubshapeAlignments(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        tgtConf=-1,
        queryConf=-1,
        pruneStats=None,
    ): ...
    def __call__(
        self,
        targetMol,
        target,
        queryMol,
        query,
        builder,
        tgtConf: int = ...,
        queryConf: int = ...,
        pruneStats: Incomplete | None = ...,
    ) -> Generator[Incomplete, None, None]: ...
