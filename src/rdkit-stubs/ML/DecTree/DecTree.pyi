"""
rdkit.ML.DecTree.DecTree module¶
Defines the class _DecTreeNode_, used to represent decision trees
_DecTreeNode_ is derived from _Tree.TreeNode_
"""
from _typeshed import Incomplete
from rdkit.ML.DecTree import Tree
from rdkit.ML.DecTree import Tree as Tree

class DecTreeNode(Tree.TreeNode):
    """
    This is used to represent decision trees
    _DecTreeNode_s are simultaneously the roots and branches of decision trees.
    Everything is nice and recursive.

    _DecTreeNode_s can save the following pieces of internal state, accessible viastandard setter/getter functions:

    _Examples_: a list of examples which have been classified
    _BadExamples_: a list of examples which have been misclassified
    _TrainingExamples_: the list of examples used to train the tree
    _TestExamples_: the list of examples used to test the tree

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

    def AddChild(self, name, label=None, data=None, isTerminal=0):
        """
        Constructs and adds a child with the specified data to our list
        Arguments

        name: the name of the new node
        label: the label of the new node (should be an integer)
        data: the data to be stored in the new node
        isTerminal: a toggle to indicate whether or not the new node is
        a terminal (leaf) node.

        **Returns*

        the _DecTreeNode_ which is constructed"""
        ...
    examples: Incomplete
    badExamples: Incomplete
    trainingExamples: Incomplete
    testExamples: Incomplete

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
    def ClearExamples(self): ...
    def GetBadExamples(self): ...
    def AddChild(self, name, label=None, data=None, isTerminal=0):
        """
        Constructs and adds a child with the specified data to our list
        Arguments

        name: the name of the new node
        label: the label of the new node (should be an integer)
        data: the data to be stored in the new node
        isTerminal: a toggle to indicate whether or not the new node is
        a terminal (leaf) node.

        **Returns*

        the _DecTreeNode_ which is constructed"""
        ...
    def GetExamples(self): ...
    def GetTestExamples(self): ...
    def GetTrainingExamples(self): ...
    def SetBadExamples(self, examples): ...
    def SetExamples(self, examples): ...
    def SetTestExamples(self, examples): ...
    def GetBadExamples(self): ...
    def SetBadExamples(self, examples): ...
    def GetTrainingExamples(self): ...
    def SetTrainingExamples(self, examples): ...
    def GetTestExamples(self): ...
    def SetTestExamples(self, examples): ...
    def ClearExamples(self): ...
