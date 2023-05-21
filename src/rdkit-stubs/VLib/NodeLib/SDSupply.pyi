"""
rdkit.VLib.NodeLib.SDSupply module¶
"""
from rdkit import Chem as Chem
from rdkit.VLib import Supply
from rdkit.VLib.Supply import SupplyNode as SupplyNode

class SDSupplyNode(Supply.SupplyNode):
    """
    SD supplier

    Sample Usage:>>> import os
    >>> from rdkit import RDConfig
    >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',                               'test_data','NCI_aids.10.sdf')
    >>> suppl = SDSupplyNode(fileN)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    10
    >>> ms[0].GetProp("_Name")
    '48'
    >>> ms[1].GetProp("_Name")
    '78'
    >>> suppl.reset()
    >>> suppl.next().GetProp("_Name")
    '48'
    >>> suppl.next().GetProp("_Name")
    '78'"""

    def next(self): ...
    def __init__(self, fileName, **kwargs) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def next(self): ...
