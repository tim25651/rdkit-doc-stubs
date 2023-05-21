"""
Module contents¶
"""
from typing import NamedTuple

from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit import rdBase as rdBase
from rdkit.Chem import rdDepictor as rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D as rdMolDraw2D
from rdkit.Chem.Draw.MolDrawing import DrawingOptions as DrawingOptions
from rdkit.Chem.Draw.MolDrawing import MolDrawing as MolDrawing
from rdkit.Chem.Draw.rdMolDraw2D import *

def MolDraw2DFromQPainter(
    qpainter,
    width: int = ...,
    height: int = ...,
    panelWidth: int = ...,
    panelHeight: int = ...,
): ...
def MolToImage(
    mol,
    size=...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    fitImage: bool = ...,
    options: Incomplete | None = ...,
    canvas: Incomplete | None = ...,
    **kwargs,
): ...
def MolToFile(
    mol,
    filename,
    size=...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    imageType: Incomplete | None = ...,
    fitImage: bool = ...,
    options: Incomplete | None = ...,
    **kwargs,
) -> None: ...
def MolToImageFile(
    mol, filename, size=..., kekulize: bool = ..., wedgeBonds: bool = ..., **kwargs
) -> None: ...
def ShowMol(
    mol,
    size=...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    title: str = ...,
    stayInFront: bool = ...,
    **kwargs,
) -> None: ...
def MolToMPL(
    mol,
    size=...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    imageType: Incomplete | None = ...,
    fitImage: bool = ...,
    options: Incomplete | None = ...,
    **kwargs,
): ...
def calcAtomGaussians(
    mol, a: float = ..., step: float = ..., weights: Incomplete | None = ...
): ...
def MolsToImage(mols, subImgSize=..., legends: Incomplete | None = ..., **kwargs): ...
def shouldKekulize(mol, kekulize): ...
def MolsToGridImage(
    mols,
    molsPerRow: int = ...,
    subImgSize=...,
    legends: Incomplete | None = ...,
    highlightAtomLists: Incomplete | None = ...,
    highlightBondLists: Incomplete | None = ...,
    useSVG: bool = ...,
    returnPNG: bool = ...,
    **kwargs,
): ...
def ReactionToImage(
    rxn,
    subImgSize=...,
    useSVG: bool = ...,
    drawOptions: Incomplete | None = ...,
    returnPNG: bool = ...,
    **kwargs,
): ...
def MolToQPixmap(
    mol,
    size=...,
    kekulize: bool = ...,
    wedgeBonds: bool = ...,
    fitImage: bool = ...,
    options: Incomplete | None = ...,
    **kwargs,
): ...
def DrawMorganBit(mol, bitId, bitInfo, whichExample: int = ..., **kwargs): ...
def DrawMorganBits(tpls, **kwargs): ...

class FingerprintEnv(NamedTuple):
    submol: Incomplete
    highlightAtoms: Incomplete
    atomColors: Incomplete
    highlightBonds: Incomplete
    bondColors: Incomplete
    highlightRadii: Incomplete

def DrawMorganEnvs(
    envs,
    molsPerRow: int = ...,
    subImgSize=...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor=...,
    ringColor=...,
    centerColor=...,
    extraColor=...,
    legends: Incomplete | None = ...,
    drawOptions: Incomplete | None = ...,
    **kwargs,
): ...
def DrawMorganEnv(
    mol,
    atomId,
    radius,
    molSize=...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor=...,
    ringColor=...,
    centerColor=...,
    extraColor=...,
    drawOptions: Incomplete | None = ...,
    **kwargs,
): ...
def DrawRDKitBits(tpls, **kwargs): ...
def DrawRDKitBit(mol, bitId, bitInfo, whichExample: int = ..., **kwargs): ...
def DrawRDKitEnvs(
    envs,
    molsPerRow: int = ...,
    subImgSize=...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor=...,
    extraColor=...,
    nonAromaticColor: Incomplete | None = ...,
    legends: Incomplete | None = ...,
    drawOptions: Incomplete | None = ...,
    **kwargs,
): ...
def DrawRDKitEnv(
    mol,
    bondPath,
    molSize=...,
    baseRad: float = ...,
    useSVG: bool = ...,
    aromaticColor=...,
    extraColor=...,
    nonAromaticColor: Incomplete | None = ...,
    drawOptions: Incomplete | None = ...,
    **kwargs,
): ...
def SetComicMode(opts) -> None: ...
