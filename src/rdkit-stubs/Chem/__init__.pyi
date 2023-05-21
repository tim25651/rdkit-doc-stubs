"""
Module contents¶
A module for molecules and stuff
see Chem/index.html in the doc tree for documentation
"""
from _typeshed import Incomplete
from rdkit import DataStructs as DataStructs
from rdkit import RDConfig as RDConfig
from rdkit import rdBase as rdBase
from rdkit.Chem import rdchem as rdchem
from rdkit.Chem import rdCoordGen as rdCoordGen
from rdkit.Chem.inchi import *
from rdkit.Chem.rdchem import *
from rdkit.Chem.rdCIPLabeler import *
from rdkit.Chem.rdmolfiles import *
from rdkit.Chem.rdMolInterchange import *
from rdkit.Chem.rdmolops import *
from rdkit.Geometry import rdGeometry as rdGeometry

templDir: Incomplete

def QuickSmartsMatch(smi, sma, unique: bool = ..., display: bool = ...): ...
def CanonSmiles(smi, useChiral: int = ...): ...
def SupplierFromFilename(fileN, delim: str = ..., **kwargs): ...
def FindMolChiralCenters(
    mol,
    force: bool = ...,
    includeUnassigned: bool = ...,
    includeCIP: bool = ...,
    useLegacyImplementation: bool = ...,
): ...
