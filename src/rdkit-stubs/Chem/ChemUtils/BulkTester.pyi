"""
rdkit.Chem.ChemUtils.BulkTester module¶
"""
import _io
from rdkit import Chem as Chem
from rdkit.Chem import Randomize as Randomize

def TestMolecule(self, mol): ...
def TestSupplier(
    self,
    suppl,
    stopAfter=-1,
    reportInterval=100,
    reportTo: _io.TextIOWrapper = ...,
    nameProp="_Name",
): ...
