"""
rdkit.VLib.Supply module¶
"""
from _typeshed import Incomplete
from rdkit.VLib import Node
from rdkit.VLib.Node import VLibNode as VLibNode

class SupplyNode(Node.VLibNode):
    """
    base class for nodes which supply things

    Assumptions:
    no parents

    Usage Example:
    >>> supplier = SupplyNode(contents=[1,2,3])
    >>> supplier.next()
    1
    >>> supplier.next()
    2
    >>> supplier.next()
    3
    >>> supplier.next()
    Traceback (most recent call last):
        ...
    StopIteration
    >>> supplier.reset()
    >>> supplier.next()
    1
    >>> [x for x in supplier]
    [1, 2, 3]"""

    def AddParent(self, parent, notify=1):
        """
        >>> p1 = VLibNode()
        >>> p2 = VLibNode()
        >>> c1 = VLibNode()
        >>> c1.AddParent(p1)
        >>> len(c1.GetParents())
        1
        >>> len(p1.GetChildren())
        1
        >>> c1.AddParent(p2,notify=0)
        >>> len(c1.GetParents())
        2
        >>> len(p2.GetChildren())
        0
        >>> p2.AddChild(c1,notify=0)
        >>> len(c1.GetParents())
        2
        >>> len(p2.GetChildren())
        1"""
        ...
    def next(self):
        """
        part of the iterator interface
        raises StopIteration on failure"""
        ...
    def __init__(self, contents: Incomplete | None = ..., **kwargs) -> None: ...
    def reset(self):
        """
        resets our iteration state"""
        ...
    def next(self):
        """
        part of the iterator interface
        raises StopIteration on failure"""
        ...
    def AddParent(self, parent, notify=1):
        """
        >>> p1 = VLibNode()
        >>> p2 = VLibNode()
        >>> c1 = VLibNode()
        >>> c1.AddParent(p1)
        >>> len(c1.GetParents())
        1
        >>> len(p1.GetChildren())
        1
        >>> c1.AddParent(p2,notify=0)
        >>> len(c1.GetParents())
        2
        >>> len(p2.GetChildren())
        0
        >>> p2.AddChild(c1,notify=0)
        >>> len(c1.GetParents())
        2
        >>> len(p2.GetChildren())
        1"""
        ...
