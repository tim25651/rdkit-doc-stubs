"""
rdkit.VLib.NodeLib.DbMolSupply module¶
"""
from rdkit import Chem as Chem
from rdkit import RDConfig as RDConfig
from rdkit.Chem.Suppliers import DbMolSupplier as DbMolSupplier
from rdkit.VLib import Supply
from rdkit.VLib.Supply import SupplyNode as SupplyNode

class DbMolSupplyNode(Supply.SupplyNode):
    """
    Supplies molecules from a db result set:

    Sample Usage:>>> from rdkit.Dbase.DbConnection import DbConnect
    >>> dbName = os.path.join(RDConfig.RDCodeDir,'Chem','Fingerprints',                             'test_data','data.gdb')
    >>> conn = DbConnect(dbName,'simple_mols')
    >>> dataset = conn.GetData()
    >>> suppl = DbMolSupplyNode(dataset)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    12
    >>> ms[0].GetProp("ID")
    'ether-1'
    >>> ms[10].GetProp("ID")
    'acid-4'
    >>> suppl.reset()
    >>> suppl.next().GetProp("ID")
    'ether-1'
    >>> suppl.next().GetProp("ID")
    'acid-1'
    >>> suppl.reset()"""

    def next(self): ...
    def __init__(self, dbResults, **kwargs) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def next(self): ...

def GetNode(self, dbName, tableName): ...
