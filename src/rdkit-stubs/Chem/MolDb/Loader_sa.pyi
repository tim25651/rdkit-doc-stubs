"""
rdkit.Chem.MolDb.Loader_sa module¶
"""
from typing import Any

import sqlalchemy
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import AllChem as AllChem
from rdkit.Chem import Crippen as Crippen
from rdkit.Chem import Descriptors as Descriptors
from rdkit.Chem import Lipinski as Lipinski

decBase: Incomplete

class Compound(sqlalchemy.orm.decl_api.Base):
    """
    A simple constructor that allows initialization from kwargs.
    Sets attributes on the constructed instance using the names and
    values in kwargs.
    Only keys that are present as
    attributes of the instance’s class are allowed. These could be,
    for example, any mapped columns or relationships."""

    __tablename__: str
    guid: Any = ...
    molpkl: Any = ...

def RegisterSchema(self, dbUrl, echo=False): ...

ConnectToSchema = RegisterSchema

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
    numForPropScan=10,
    startAnew=True,
): ...

logger: Incomplete

def ProcessMol(
    self,
    session,
    mol,
    globalProps,
    nDone,
    nameProp="_Name",
    nameCol="compound_id",
    redraw=False,
    keepHs=False,
    skipProps=False,
    addComputedProps=False,
    skipSmiles=False,
): ...
def RegisterSchema(self, dbUrl, echo=False): ...
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
    numForPropScan=10,
    startAnew=True,
): ...
