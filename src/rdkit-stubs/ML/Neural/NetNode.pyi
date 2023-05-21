"""
rdkit.ML.Neural.NetNode module¶
Contains the class _NetNode_ which is used to represent nodes in neural nets
Network Architecture:

A tacit assumption in all of this stuff is that we’re dealing with
feedforward networks.
The network itself is stored as a list of _NetNode_ objects.  The list
is ordered in the sense that nodes in earlier/later layers than a
given node are guaranteed to come before/after that node in the list.
This way we can easily generate the values of each node by moving
sequentially through the list, we’re guaranteed that every input for a
node has already been filled in.
Each node stores a list (_inputNodes_) of indices of its inputs in the
main node list.
"""
from _typeshed import Incomplete

from . import ActFuncs as ActFuncs

class NetNode(object):
    """
    a node in a neural network
    Constructor
    Arguments

    nodeIndex: the integer index of this node in _nodeList_
    nodeList: the list of other _NetNodes_ already in the network
    inputNodes: a list of this node’s inputs
    weights: a list of this node’s weights

    actFunc: the activation function to be used here.  Must support the APIof _ActFuncs.ActFunc_.

    actFuncParms: a tuple of extra arguments to be passed to the activation functionconstructor.

    NoteThere should be only one copy of _inputNodes_, every _NetNode_ just has a pointer
    to it so that changes made at one node propagate automatically to the others."""

    def Eval(self, valVect):
        """
        Given a set of inputs (valVect), returns the output of this node
        Arguments

        valVect: a list of inputs

        Returns

        the result of running the values in valVect through this node"""
        ...
    def GetInputs(self):
        """
        returns the input list"""
        ...
    def GetWeights(self):
        """
        returns the weight list"""
        ...
    inputNodes: Incomplete

    def SetInputs(self, inputNodes):
        """
        Sets the input list
        Arguments

        inputNodes: a list of _NetNode_s which are to be used as inputs

        Note

        If this _NetNode_ already has weights set and _inputNodes_ is a different length,
        this will bomb out with an assertion."""
        ...
    def GetInputs(self):
        """
        returns the input list"""
        ...
    weights: Incomplete

    def SetWeights(self, weights):
        """
        Sets the weight list
        Arguments

        weights: a list of values which are to be used as weights

        Note

        If this _NetNode_ already has _inputNodes_  and _weights_ is a different length,
        this will bomb out with an assertion."""
        ...
    def GetWeights(self):
        """
        returns the weight list"""
        ...
    nodeIndex: Incomplete
    nodeList: Incomplete
    actFunc: Incomplete

    def __init__(
        self,
        nodeIndex,
        nodeList,
        inputNodes: Incomplete | None = ...,
        weights: Incomplete | None = ...,
        actFunc=...,
        actFuncParms=...,
    ) -> None: ...
