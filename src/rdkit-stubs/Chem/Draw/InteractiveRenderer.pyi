from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Draw as Draw

from . import rdMolDraw2D as rdMolDraw2D

log: Incomplete
rdkitStructureRendererJsUrl: str
minimalLibJsUrl: str
parentNodeQuery: str

def filterDefaultDrawOpts(molDrawOptions): ...
def setEnabled(shouldEnable: bool = ..., quiet: bool = ...): ...
def isEnabled(mol: Incomplete | None = ...): ...
def getOpts(mol): ...
def setOpts(mol, opts) -> None: ...
def setOpt(mol, key, value) -> None: ...
def clearOpts(mol) -> None: ...
def clearOpt(mol, key) -> None: ...
def setNoInteractive(mol, shouldSet: bool = ...) -> None: ...
def isNoInteractive(mol): ...
def injectHTMLFooterAfterTable(html): ...
def generateHTMLFooter(doc, element): ...
def newlineToXml(molblock): ...
def xmlToNewline(xmlblock): ...
def toDataMol(mol): ...
def camelCaseOptToDataTag(opt): ...
def generateHTMLBody(mol, size, **kwargs): ...
def MolsToHTMLTable(
    mols,
    molsPerRow: int = ...,
    subImgSize=...,
    legends: Incomplete | None = ...,
    highlightAtomLists: Incomplete | None = ...,
    highlightBondLists: Incomplete | None = ...,
    useSVG: bool = ...,
    drawOptions: Incomplete | None = ...,
    **kwargs,
): ...
