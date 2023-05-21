"""
rdkit.VLib.NodeLib.SmilesSupply module¶
"""
from rdkit import Chem as Chem
from rdkit.VLib import Supply
from rdkit.VLib.Supply import SupplyNode as SupplyNode

class SmilesSupplyNode(Supply.SupplyNode):
    """
    Smiles supplier

    Sample Usage:>>> import os
    >>> from rdkit import RDConfig
    >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',                               'test_data','pgp_20.txt')
    >>> suppl = SmilesSupplyNode(fileN,delim="	",smilesColumn=2,nameColumn=1,titleLine=1)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    20
    >>> ms[0].GetProp("_Name")
    'ALDOSTERONE'
    >>> ms[0].GetProp("ID")
    'RD-PGP-0001'
    >>> ms[1].GetProp("_Name")
    'AMIODARONE'
    >>> ms[3].GetProp("ID")
    'RD-PGP-0004'
    >>> suppl.reset()
    >>> suppl.next().GetProp("_Name")
    'ALDOSTERONE'
    >>> suppl.next().GetProp("_Name")
    'AMIODARONE'
    >>> suppl.reset()"""

    def next(self): ...
    def __init__(
        self,
        fileName,
        delim: str = ...,
        nameColumn: int = ...,
        smilesColumn: int = ...,
        titleLine: int = ...,
        **kwargs,
    ) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def next(self): ...
