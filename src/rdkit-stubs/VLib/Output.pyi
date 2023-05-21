"""
rdkit.VLib.Output module¶
"""
from _typeshed import Incomplete
from rdkit.VLib import Node
from rdkit.VLib.Node import VLibNode as VLibNode

class OutputNode(Node.VLibNode):
    """
    base class for nodes which dump output
    Assumptions:

    destination supports a write() method
    strFunc, if provided, returns a string representation of
    the input
    inputs (parents) can be stepped through in lockstep

    Usage Example:
    >>> from rdkit.VLib.Supply import SupplyNode
    >>> supplier = SupplyNode(contents=[1,2,3])
    >>> from io import StringIO
    >>> sio = StringIO()
    >>> node = OutputNode(dest=sio,strFunc=lambda x:'%s '%(str(x)))
    >>> node.AddParent(supplier)
    >>> node.next()
    1
    >>> sio.getvalue()
    '1 '
    >>> node.next()
    2
    >>> sio.getvalue()
    '1 2 '"""

    def __init__(
        self, dest: Incomplete | None = ..., strFunc: Incomplete | None = ..., **kwargs
    ) -> None: ...
    def next(self):
        """
        part of the iterator interface
        raises StopIteration on failure"""
        ...
