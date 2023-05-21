"""
rdkit.VLib.NodeLib.SmilesDupeFilter module¶
"""
from rdkit import Chem as Chem
from rdkit.VLib import Filter
from rdkit.VLib.Filter import FilterNode as FilterNode

class DupeFilter(Filter.FilterNode):
    """
    canonical-smiles based duplicate filter
    Assumptions:

    inputs are molecules

    Sample Usage:>>> import os
    >>> from rdkit import RDConfig
    >>> from rdkit.VLib.NodeLib.SDSupply import SDSupplyNode
    >>> fileN = os.path.join(RDConfig.RDCodeDir,'VLib','NodeLib',                             'test_data','NCI_aids.10.sdf')
    >>> suppl = SDSupplyNode(fileN)
    >>> filt = DupeFilter()
    >>> filt.AddParent(suppl)
    >>> ms = [x for x in filt]
    >>> len(ms)
    10
    >>> ms[0].GetProp("_Name")
    '48'
    >>> ms[1].GetProp("_Name")
    '78'
    >>> filt.reset()
    >>> filt.next().GetProp("_Name")
    '48'"""

    def filter(self, cmpd): ...
    def __init__(self, **kwargs) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def filter(self, cmpd): ...
