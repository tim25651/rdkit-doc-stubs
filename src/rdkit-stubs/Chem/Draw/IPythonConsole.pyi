"""
rdkit.Chem.Draw.IPythonConsole module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Draw as Draw
from rdkit.Chem import rdchem as rdchem
from rdkit.Chem import rdChemReactions as rdChemReactions
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D

from . import InteractiveRenderer as InteractiveRenderer

def DisableSubstructMatchRendering(self): ...

molSize: Incomplete
highlightSubstructs: bool
kekulizeStructures: bool
highlightByReactant: bool
ipython_useSVG: bool
ipython_showProperties: bool
ipython_maxProperties: int
ipython_3d: bool
molSize_3d: Incomplete
drawing_type_3d: str
bgcolor_3d: str
drawOptions: Incomplete

def addMolToView(self, mol, view, confId=-1, drawAs=None): ...
def drawMol3D(self, m, view=None, confId=-1, drawAs=None, bgColor=None, size=None): ...
def display_pil_image(self, img):
    """
    displayhook function for PIL Images, rendered as PNG"""
    ...

def ShowMols(self, mols, maxMols=50, **kwargs): ...
def DrawMorganBit(
    self, mol, bitId, bitInfo, drawOptions: rdMolDraw2D.MolDrawOptions = ..., **kwargs
): ...
def DrawMorganBits(
    self, *args, drawOptions: rdMolDraw2D.MolDrawOptions = ..., **kwargs
): ...
def DrawRDKitBit(
    self, mol, bitId, bitInfo, drawOptions: rdMolDraw2D.MolDrawOptions = ..., **kwargs
): ...
def DrawRDKitBits(
    self, *args, drawOptions: rdMolDraw2D.MolDrawOptions = ..., **kwargs
): ...
def EnableSubstructMatchRendering(self): ...
def InstallIPythonRenderer(self): ...
def ShowMols(self, mols, maxMols=50, **kwargs): ...
def DisableSubstructMatchRendering(self): ...
def UninstallIPythonRenderer(self): ...
def addMolToView(self, mol, view, confId=-1, drawAs=None): ...
def display_pil_image(self, img):
    """
    displayhook function for PIL Images, rendered as PNG"""
    ...

def drawMol3D(self, m, view=None, confId=-1, drawAs=None, bgColor=None, size=None): ...
