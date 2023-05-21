"""
rdkit.Chem.ChemUtils.TemplateExpand module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem
from rdkit.Chem import Crippen as Crippen
from rdkit.Chem.ChemUtils.AlignDepict import AlignDepict as AlignDepict

def ConstructSidechains(self, suppl, sma=None, replace=True, useAll=False): ...

logger: Incomplete

def Usage(self): ...

nDumped: int

def Explode(
    self, template, sidechains, outF, autoNames=True, do3D=False, useTethers=False
): ...
def MoveDummyNeighborsToBeginning(self, mol, useAll=False): ...
def Usage(self): ...
def ConstructSidechains(self, suppl, sma=None, replace=True, useAll=False): ...
