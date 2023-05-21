"""
rdkit.Chem.Pharm3D.Pharmacophore module¶
"""
from _typeshed import Incomplete
from rdkit import Geometry as Geometry
from rdkit.Chem import ChemicalFeatures as ChemicalFeatures
from rdkit.RDLogger import logger as logger

class ExplicitPharmacophore(object):
    """
    this is a pharmacophore with explicit point locations and radii"""

    def getFeature(self, i): ...
    def __init__(
        self, feats: Incomplete | None = ..., radii: Incomplete | None = ...
    ) -> None: ...
    def getFeatures(self): ...
    def getRadii(self): ...
    def getFeature(self, i): ...
    def getRadius(self, i): ...
    def setRadius(self, i, rad): ...
    def initFromString(self, text): ...
    def initFromFile(self, inF): ...
    def initFromLines(self, lines): ...
    def initFromString(self, text): ...
    def setRadius(self, i, rad): ...

class Pharmacophore(object):
    def getFeature(self, i): ...
    def __init__(self, feats, initMats: bool = ...) -> None: ...
    def getFeatures(self): ...
    def getLowerBound(self, i, j): ...
    def getLowerBound2D(self, i, j): ...
    def getFeature(self, i): ...
    def getUpperBound(self, i, j): ...
    def getUpperBound2D(self, i, j): ...
    def setLowerBound(self, i, j, val, checkBounds=False): ...
    def setLowerBound2D(self, i, j, val, checkBounds=False): ...
    def getLowerBound(self, i, j): ...
    def setUpperBound(self, i, j, val, checkBounds=False): ...
    def setLowerBound(self, i, j, val, checkBounds=False): ...
    def getUpperBound2D(self, i, j): ...
    def getLowerBound2D(self, i, j): ...
    def setUpperBound2D(self, i, j, val, checkBounds=False): ...
    def setLowerBound2D(self, i, j, val, checkBounds=False): ...

class ExplicitPharmacophore(object):
    """
    this is a pharmacophore with explicit point locations and radii"""

    def getFeature(self, i): ...
    def __init__(
        self, feats: Incomplete | None = ..., radii: Incomplete | None = ...
    ) -> None: ...
    def getFeatures(self): ...
    def getRadii(self): ...
    def getFeature(self, i): ...
    def getRadius(self, i): ...
    def setRadius(self, i, rad): ...
    def initFromString(self, text): ...
    def initFromFile(self, inF): ...
    def initFromLines(self, lines): ...
    def initFromString(self, text): ...
    def setRadius(self, i, rad): ...