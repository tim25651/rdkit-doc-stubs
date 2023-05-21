"""
rdkit.VLib.NodeLib.SmartsMolFilter module¶
"""
from rdkit import Chem as Chem
from rdkit.VLib import Filter
from rdkit.VLib.Filter import FilterNode as FilterNode

class SmartsFilter(Filter.FilterNode):
    """
    filter out molecules matching one or more SMARTS patterns
    There is a count associated with each pattern. Molecules are
    allowed to match the pattern up to this number of times.
    Assumptions:

    inputs are molecules

    Sample Usage:>>> smis = ['C1CCC1','C1CCC1C=O','CCCC','CCC=O','CC(=O)C','CCN','NCCN','NCC=O']
    >>> mols = [Chem.MolFromSmiles(x) for x in smis]
    >>> from rdkit.VLib.Supply import SupplyNode
    >>> suppl = SupplyNode(contents=mols)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    8

    We can pass in SMARTS strings:
    >>> smas = [‘C=O’,’CN’]
    >>> counts = [1,2]
    >>> filt = SmartsFilter(patterns=smas,counts=counts)
    >>> filt.AddParent(suppl)
    >>> ms = [x for x in filt]
    >>> len(ms)
    5
    Alternatively, we can pass in molecule objects:
    >>> mols =[Chem.MolFromSmarts(x) for x in smas]
    >>> counts = [1,2]
    >>> filt.Destroy()
    >>> filt = SmartsFilter(patterns=mols,counts=counts)
    >>> filt.AddParent(suppl)
    >>> ms = [x for x in filt]
    >>> len(ms)
    5
    Negation does what you’d expect:
    >>> filt.SetNegate(1)
    >>> ms = [x for x in filt]
    >>> len(ms)
    3"""

    def __init__(self, patterns=..., counts=..., **kwargs) -> None: ...
    def filter(self, cmpd): ...
