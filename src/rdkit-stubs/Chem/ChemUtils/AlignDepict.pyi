"""
rdkit.Chem.ChemUtils.AlignDepict module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import Geometry as Geometry
from rdkit.Chem import rdDepictor as rdDepictor

def AlignDepict(self, mol, core, corePattern=None, acceptFailure=False):
    """
    Arguments:

    mol:          the molecule to be aligned, this will come backwith a single conformer.

    core:         a molecule with the core atoms to align to;this should have a depiction.

    corePattern:  (optional) an optional molecule to be used togenerate the atom mapping between the molecule
    and the core."""
    ...

def initParser(self):
    """
    Initialize the parser"""
    ...

def main(self):
    """
    Main application"""
    ...

def processArgs(self, args): ...
def main(self):
    """
    Main application"""
    ...
