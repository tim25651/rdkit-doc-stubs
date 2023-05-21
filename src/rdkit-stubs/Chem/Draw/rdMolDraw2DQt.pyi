"""
rdkit.Chem.Draw.rdMolDraw2DQt module¶
Module containing a C++ implementation of 2D molecule drawing using Qt
"""
from rdkit.Chem.Draw import rdMolDraw2D

class MolDraw2DQt(rdMolDraw2D.MolDraw2D):
    """
    Qt molecule drawer
    Raises an exception
    This class cannot be instantiated from Python"""

    ...

def MolDraw2DFromQPainter_(
    self,
    width: int,
    height: int,
    pointer_to_QPainter: int,
    panelWidth: int = -1,
    panelHeight: int = -1,
) -> MolDraw2DQt:
    """
    Returns a MolDraw2DQt instance set to use a QPainter.
    Use sip.unwrapinstance(qptr) to get the required pointer information. Please note that this is somewhat fragile.

    C++ signature :RDKit::MolDraw2DQt* MolDraw2DFromQPainter_(int,int,unsigned long [,int=-1 [,int=-1]])
    """
    ...
