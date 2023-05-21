"""
rdkit.ML.DecTree.Tree module¶
Implements a class used to represent N-ary trees
"""
from _typeshed import Incomplete

class TreeNode(object):
    """
    This is your bog standard Tree class.
    the root of the tree is just a TreeNode like all other members.
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
        Creates a new TreeNode and adds a child to the tree
        Arguments

        name: the name of the new node
        label: the label of the new node (should be an integer)
        data: the data to be stored in the new node
        isTerminal: a toggle to indicate whether or not the new node is
        a terminal (leaf) node.

        **Returns*

        the _TreeNode_ which is constructed"""
        ...
    def AddChildNode(self, node):
        """
        Adds a TreeNode to the local list of children
        Arguments

        node: the node to be added

        Note

        the level of the node (used in printing) is set as well"""
        ...
    def Destroy(self):
        """
        Destroys this node and all of its children"""
        ...
    def GetChildren(self):
        """
        Returns a python list of the children of this node"""
        ...
    def GetData(self):
        """
        Returns the data stored at this node"""
        ...
    def GetLabel(self):
        """
        Returns the label of this node"""
        ...
    def GetLevel(self):
        """
        Returns the level of this node"""
        ...
    def GetName(self):
        """
        Returns the name of this node"""
        ...
    def GetParent(self):
        """
        Returns the parent of this node"""
        ...
    def GetTerminal(self):
        """
        Returns whether or not this node is terminal"""
        ...
    NameModel = NameTree
    children: Incomplete
    parent: Incomplete
    name: Incomplete
    data: Incomplete
    terminalNode: Incomplete
    label: Incomplete
    level: Incomplete
    examples: Incomplete

    def __init__(
        self,
        parent,
        name,
        label: Incomplete | None = ...,
        data: Incomplete | None = ...,
        level: int = ...,
        isTerminal: int = ...,
    ) -> None: ...
    def NameTree(self, varNames):
        """
        Set the names of each node in the tree from a list of variable names.
        Arguments

        varNames: a list of names to be assigned

        Notes

        this works its magic by recursively traversing all children
        The assumption is made here that the varNames list can be indexed
        by the labels of tree nodes"""
        ...
    def Pickle(self, fileName="foo.pkl"):
        """
        Pickles the tree and writes it to disk"""
        ...
    def Print(self, level=0, showData=0):
        """
        Pretty prints the tree
        Arguments

        level: sets the number of spaces to be added at the beginning of the output
        showData: if this is nonzero, the node’s _data_ value will be printed as well

        Note

        this works recursively"""
        ...
    NameModel = NameTree

    def AddChildNode(self, node):
        """
        Adds a TreeNode to the local list of children
        Arguments

        node: the node to be added

        Note

        the level of the node (used in printing) is set as well"""
        ...
    def AddChild(self, name, label=None, data=None, isTerminal=0):
        """
        Creates a new TreeNode and adds a child to the tree
        Arguments

        name: the name of the new node
        label: the label of the new node (should be an integer)
        data: the data to be stored in the new node
        isTerminal: a toggle to indicate whether or not the new node is
        a terminal (leaf) node.

        **Returns*

        the _TreeNode_ which is constructed"""
        ...
    def PruneChild(self, child):
        """
        Removes the child node
        Arguments

        child: a TreeNode"""
        ...
    def ReplaceChildIndex(self, index, newChild):
        """
        Replaces a given child with a new one
        Arguments

        index: an integer
        child: a TreeNode"""
        ...
    def GetChildren(self):
        """
        Returns a python list of the children of this node"""
        ...
    def Destroy(self):
        """
        Destroys this node and all of its children"""
        ...
    def GetName(self):
        """
        Returns the name of this node"""
        ...
    def SetName(self, name):
        """
        Sets the name of this node"""
        ...
    def GetData(self):
        """
        Returns the data stored at this node"""
        ...
    def SetData(self, data):
        """
        Sets the data stored at this node"""
        ...
    def GetTerminal(self):
        """
        Returns whether or not this node is terminal"""
        ...
    def SetTerminal(self, isTerminal):
        """
        Sets whether or not this node is terminal"""
        ...
    def GetLabel(self):
        """
        Returns the label of this node"""
        ...
    def SetLabel(self, label):
        """
        Sets the label of this node (should be an integer)"""
        ...
    def GetLevel(self):
        """
        Returns the level of this node"""
        ...
    def SetLevel(self, level):
        """
        Sets the level of this node"""
        ...
    def SetName(self, name):
        """
        Sets the name of this node"""
        ...
    def GetParent(self):
        """
        Returns the parent of this node"""
        ...
    def SetParent(self, parent):
        """
        Sets the parent of this node"""
        ...
    def SetTerminal(self, isTerminal):
        """
        Sets whether or not this node is terminal"""
        ...
    def Print(self, level=0, showData=0):
        """
        Pretty prints the tree
        Arguments

        level: sets the number of spaces to be added at the beginning of the output
        showData: if this is nonzero, the node’s _data_ value will be printed as well

        Note

        this works recursively"""
        ...
    def Pickle(self, fileName="foo.pkl"):
        """
        Pickles the tree and writes it to disk"""
        ...
    def __cmp__(self, other): ...
    def __lt__(self, other): ...
    def __eq__(self, other): ...
