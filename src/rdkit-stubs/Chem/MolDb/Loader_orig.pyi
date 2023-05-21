"""
rdkit.Chem.MolDb.Loader_orig module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem
from rdkit.Chem import Crippen as Crippen
from rdkit.Chem import Descriptors as Descriptors
from rdkit.Chem import Lipinski as Lipinski
from rdkit.Dbase import DbModule as DbModule
from rdkit.Dbase.DbConnection import DbConnect as DbConnect

logger: Incomplete

def ProcessMol(
    self,
    mol,
    typeConversions,
    globalProps,
    nDone,
    nameProp="_Name",
    nameCol="compound_id",
    redraw=False,
    keepHs=False,
    skipProps=False,
    addComputedProps=False,
    skipSmiles=False,
    uniqNames=None,
    namesSeen=None,
): ...
def ConvertRows(self, rows, globalProps, defaultVal, skipSmiles): ...
def LoadDb(
    self,
    suppl,
    dbName,
    nameProp="_Name",
    nameCol="compound_id",
    silent=False,
    redraw=False,
    errorsTo=None,
    keepHs=False,
    defaultVal="N/A",
    skipProps=False,
    regName="molecules",
    skipSmiles=False,
    maxRowsCached=-1,
    uniqNames=False,
    addComputedProps=False,
    lazySupplier=False,
    startAnew=True,
): ...
def ProcessMol(
    self,
    mol,
    typeConversions,
    globalProps,
    nDone,
    nameProp="_Name",
    nameCol="compound_id",
    redraw=False,
    keepHs=False,
    skipProps=False,
    addComputedProps=False,
    skipSmiles=False,
    uniqNames=None,
    namesSeen=None,
): ...
