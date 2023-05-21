"""
rdkit.VLib.Transform module¶
"""
from _typeshed import Incomplete
from rdkit.VLib import Node
from rdkit.VLib.Node import VLibNode as VLibNode

class TransformNode(Node.VLibNode):
    """
    base class for nodes which filter their input
    Assumptions:

    transform function takes a number of arguments equal to the
    number of inputs we have.  We return whatever it returns
    inputs (parents) can be stepped through in lockstep

    Usage Example:
    >>> from rdkit.VLib.Supply import SupplyNode
    >>> def func(a,b):
    ...   return a+b
    >>> tform = TransformNode(func)
    >>> suppl1 = SupplyNode(contents=[1,2,3,3])
    >>> suppl2 = SupplyNode(contents=[1,2,3,1])
    >>> tform.AddParent(suppl1)
    >>> tform.AddParent(suppl2)
    >>> v = [x for x in tform]
    >>> v
    [2, 4, 6, 4]
    >>> tform.reset()
    >>> v = [x for x in tform]
    >>> v
    [2, 4, 6, 4]

    If we don’t provide a function, just return the inputs:
    >>> tform = TransformNode()
    >>> suppl1 = SupplyNode(contents=[1,2,3,3])
    >>> suppl2 = SupplyNode(contents=[1,2,3,1])
    >>> tform.AddParent(suppl1)
    >>> tform.AddParent(suppl2)
    >>> v = [x for x in tform]
    >>> v
    [(1, 1), (2, 2), (3, 3), (3, 1)]"""

    def __init__(self, func: Incomplete | None = ..., **kwargs) -> None: ...
    def next(self):
        """
        part of the iterator interface
        raises StopIteration on failure"""
        ...
