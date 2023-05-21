"""
rdkit.VLib.NodeLib.SmartsRemover module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.VLib import Transform
from rdkit.VLib.Transform import TransformNode as TransformNode

class SmartsRemover(Transform.TransformNode):
    """
    transforms molecules by removing atoms matching smarts patterns
    Assumptions:

    inputs are molecules

    Sample Usage:>>> smis = ['C1CCC1.C=O','C1CCC1C=O','CCC=O.C=O','NCC=O.C=O.CN']
    >>> mols = [Chem.MolFromSmiles(x) for x in smis]
    >>> from rdkit.VLib.Supply import SupplyNode
    >>> suppl = SupplyNode(contents=mols)
    >>> ms = [x for x in suppl]
    >>> len(ms)
    4

    We can pass in SMARTS strings:
    >>> smas = [‘C=O’,’CN’]
    >>> tform = SmartsRemover(patterns=smas)
    >>> tform.AddParent(suppl)
    >>> ms = [x for x in tform]
    >>> len(ms)
    4
    >>> Chem.MolToSmiles(ms[0])
    ‘C1CCC1’
    >>> Chem.MolToSmiles(ms[1])
    ‘O=CC1CCC1’
    >>> Chem.MolToSmiles(ms[2])
    ‘CCC=O’
    >>> Chem.MolToSmiles(ms[3])
    ‘NCC=O’
    We can also remove pieces of the molecule that are not complete
    fragments:
    >>> tform.Destroy()
    >>> smas = [‘C=O’,’CN’]
    >>> smas = [Chem.MolFromSmarts(x) for x in smas]
    >>> tform = SmartsRemover(patterns=smas,wholeFragments=0)
    >>> tform.AddParent(suppl)
    >>> ms = [x for x in tform]
    >>> len(ms)
    4
    >>> Chem.MolToSmiles(ms[0])
    ‘C1CCC1’
    >>> Chem.MolToSmiles(ms[1])
    ‘C1CCC1’
    >>> Chem.MolToSmiles(ms[3])
    ‘’
    Or patterns themselves:
    >>> tform.Destroy()
    >>> smas = [‘C=O’,’CN’]
    >>> smas = [Chem.MolFromSmarts(x) for x in smas]
    >>> tform = SmartsRemover(patterns=smas)
    >>> tform.AddParent(suppl)
    >>> ms = [x for x in tform]
    >>> len(ms)
    4
    >>> Chem.MolToSmiles(ms[0])
    ‘C1CCC1’
    >>> Chem.MolToSmiles(ms[3])
    ‘NCC=O’"""

    def __init__(self, patterns=..., wholeFragments: int = ..., **kwargs) -> None: ...
    def transform(self, cmpd): ...

biggerTest: str
__test__: Incomplete
