"""
rdkit.ML.DecTree.QuantTree module¶
Defines the class _QuantTreeNode_, used to represent decision trees with automatic
quantization bounds

_QuantTreeNode_ is derived from _DecTree.DecTreeNode_
"""
from _typeshed import Incomplete
from rdkit.ML.DecTree import DecTree
from rdkit.ML.DecTree import DecTree as DecTree
from rdkit.ML.DecTree import Tree as Tree

class QuantTreeNode(DecTree.DecTreeNode):
    """
    constructor
    Arguments

    parent: the parent of this node in the tree
    name: the name of the node
    label: the node’s label (should be an integer)
    data: an optional data field
    level: an integer indicating the level of this node in the hierarchy
    (used for printing)
    isTerminal: flags a node as being terminal.  This is useful for those
    times when it’s useful to know such things."""

    qBounds: Incomplete
    nBounds: int

    def __init__(self, *args, **kwargs) -> None: ...
    def ClassifyExample(self, example, appendExamples=0):
        """
        Recursively classify an example by running it through the tree
        Arguments

        example: the example to be classified
        appendExamples: if this is nonzero then this node (and all children)
        will store the example

        Returns

        the classification of _example_

        NOTE:In the interest of speed, I don’t use accessor functions
        here.  So if you subclass DecTreeNode for your own trees, you’ll
        have to either include ClassifyExample or avoid changing the names
        of the instance variables this needs."""
        ...
    def GetQuantBounds(self): ...
    def SetQuantBounds(self, qBounds): ...
    def GetQuantBounds(self): ...
    def __cmp__(self, other): ...
    def __lt__(self, other): ...
    def __eq__(self, other): ...
