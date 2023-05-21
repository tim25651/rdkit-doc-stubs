"""
rdkit.Chem.Subshape.BuilderUtils module¶
"""
from rdkit import Geometry as Geometry
from rdkit.Chem.Subshape import SubshapeObjects as SubshapeObjects

def AppendSkeletonPoints(
    self,
    shapeGrid,
    termPts,
    winRad,
    stepDist,
    maxGridVal=3,
    maxDistC=15.0,
    distTol=1.5,
    symFactor=1.5,
    verbose=False,
): ...
def AssignMolFeatsToPoints(self, pts, mol, featFactory, winRad): ...
def CalculateDirectionsAtPoint(self, pt, shapeGrid, winRad): ...
def ClusterTerminalPts(self, pts, winRad, scale): ...
def ComputeGridIndices(self, shapeGrid, winRad): ...
def ComputeShapeGridCentroid(self, pt, shapeGrid, winRad): ...
def ExpandTerminalPts(self, shape, pts, winRad, maxGridVal=3.0, targetNumPts=5):
    """
    find additional terminal points until a target number is reached"""
    ...

def FindFarthestGridPoint(self, shape, loc, winRad, maxGridVal):
    """
    find the grid point with max occupancy that is furthest from a
    given location"""
    ...

def FindGridPointBetweenPoints(self, pt1, pt2, shapeGrid, winRad): ...
def FindTerminalPtsFromConformer(self, conf, winRad, nbrCount): ...
def FindTerminalPtsFromShape(self, shape, winRad, fraction, maxGridVal=3): ...
def FindTerminalPtsFromConformer(self, conf, winRad, nbrCount): ...
def FindGridPointBetweenPoints(self, pt1, pt2, shapeGrid, winRad): ...
def ClusterTerminalPts(self, pts, winRad, scale): ...
def GetMoreTerminalPoints(self, shape, pts, winRad, maxGridVal, targetNumber=5):
    """
    adds a set of new terminal points using a max-min algorithm"""
    ...

def FindFarthestGridPoint(self, shape, loc, winRad, maxGridVal):
    """
    find the grid point with max occupancy that is furthest from a
    given location"""
    ...

def ExpandTerminalPts(self, shape, pts, winRad, maxGridVal=3.0, targetNumPts=5):
    """
    find additional terminal points until a target number is reached"""
    ...

def AppendSkeletonPoints(
    self,
    shapeGrid,
    termPts,
    winRad,
    stepDist,
    maxGridVal=3,
    maxDistC=15.0,
    distTol=1.5,
    symFactor=1.5,
    verbose=False,
): ...
def CalculateDirectionsAtPoint(self, pt, shapeGrid, winRad): ...
def AssignMolFeatsToPoints(self, pts, mol, featFactory, winRad): ...
