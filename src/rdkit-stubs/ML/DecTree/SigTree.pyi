"""
rdkit.ML.DecTree.SigTree module¶
Defines the class SigTreeNode, used to represent trees that
use signatures (bit vectors) to represent data.  As inputs (examples),
SigTreeNode’s expect 3-sequences: (label,sig,act)

_SigTreeNode_ is derived from _DecTree.DecTreeNode_
"""
from rdkit.DataStructs.VectCollection import VectCollection as VectCollection
from rdkit.ML.DecTree import DecTree
from rdkit.ML.DecTree import DecTree as DecTree

class SigTreeNode(DecTree.DecTreeNode):
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

    def ClassifyExample(self, example, appendExamples=0):
        """
        Recursively classify an example by running it through the tree
        Arguments

        example: the example to be classified, a sequence at least
        2 long:

        ( id, sig )

        where sig is a BitVector (or something supporting __getitem__)
        additional fields will be ignored.

        appendExamples: if this is nonzero then this node (and all children)
        will store the example

        Returns

        the classification of _example_"""
        ...
    def NameModel(self, *args, **kwargs):
        """
        Set the names of each node in the tree from a list of variable names.
        Arguments

        varNames: a list of names to be assigned

        Notes

        this works its magic by recursively traversing all children
        The assumption is made here that the varNames list can be indexed
        by the labels of tree nodes"""
        ...
    def ClassifyExample(self, example, appendExamples=0):
        """
        Recursively classify an example by running it through the tree
        Arguments

        example: the example to be classified, a sequence at least
        2 long:

        ( id, sig )

        where sig is a BitVector (or something supporting __getitem__)
        additional fields will be ignored.

        appendExamples: if this is nonzero then this node (and all children)
        will store the example

        Returns

        the classification of _example_"""
        ...
    def __init__(self, *args, **kwargs) -> None: ...
