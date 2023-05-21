"""
rdkit.ML.Neural.Network module¶
Contains the class _Network_ which is used to represent neural nets
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
from rdkit.ML.Neural import ActFuncs as ActFuncs
from rdkit.ML.Neural import NetNode as NetNode

class Network(object):
    """
    a neural network
    Constructor
    This constructs and initializes the network based upon the specified
    node counts.
    A fully connected network with random weights is constructed.
    Arguments

    nodeCounts: a list containing the number of nodes to be in each layer.
    the ordering is:(nInput,nHidden1,nHidden2, … , nHiddenN, nOutput)

    nodeConnections: I don’t know why this is here, but it’s optional.  ;-)

    actFunc: the activation function to be used here.  Must support the APIof _ActFuncs.ActFunc_.

    actFuncParms: a tuple of extra arguments to be passed to the activation functionconstructor.

    weightBounds:  a float which provides the boundary on the random initial weights"""

    def ClassifyExample(self, example, appendExamples=0):
        """
        classifies a given example and returns the results of the output layer.
        Arguments

        example: the example to be classified

        NOTE:

        if the output layer is only one element long,
        a scalar (not a list) will be returned.  This is why a lot of the other
        network code claims to only support single valued outputs."""
        ...
    def ConstructNodes(self, nodeCounts, actFunc, actFuncParms):
        """
        build an unconnected network and set node counts
        Arguments

        nodeCounts: a list containing the number of nodes to be in each layer.
        the ordering is:(nInput,nHidden1,nHidden2, … , nHiddenN, nOutput)"""
        ...
    def ConstructRandomWeights(self, minWeight=-1, maxWeight=1):
        """
        initialize all the weights in the network to random numbers
        Arguments

        minWeight: the minimum value a weight can take
        maxWeight: the maximum value a weight can take"""
        ...
    nConnections: Incomplete

    def FullyConnectNodes(self):
        """
        Fully connects each layer in the network to the one above it

        Notethis sets the connections, but does not assign weights"""
        ...
    def GetAllNodes(self):
        """
        returns a list of all nodes"""
        ...
    def GetHiddenLayerNodeList(self, which):
        """
        returns a list of hidden nodes in the specified layer"""
        ...
    nodeCounts: Incomplete
    numInputNodes: Incomplete
    numOutputNodes: Incomplete
    numHiddenLayers: Incomplete
    numInHidden: Incomplete
    nodeList: Incomplete
    layerIndices: Incomplete

    def ConstructNodes(self, nodeCounts, actFunc, actFuncParms):
        """
        build an unconnected network and set node counts
        Arguments

        nodeCounts: a list containing the number of nodes to be in each layer.
        the ordering is:(nInput,nHidden1,nHidden2, … , nHiddenN, nOutput)"""
        ...
    def GetInputNodeList(self):
        """
        returns a list of input node indices"""
        ...
    def GetLastOutputs(self):
        """
        returns the complete list of output layer values from the last time this node
        classified anything"""
        ...
    def GetNode(self, which):
        """
        returns a particular node"""
        ...
    def GetNumHidden(self):
        """
        returns the number of hidden layers"""
        ...
    def GetNumNodes(self):
        """
        returns the total number of nodes"""
        ...
    def GetOutputNodeList(self):
        """
        returns a list of output node indices"""
        ...
    def GetHiddenLayerNodeList(self, which):
        """
        returns a list of hidden nodes in the specified layer"""
        ...
    def GetNumNodes(self):
        """
        returns the total number of nodes"""
        ...
    def GetNumHidden(self):
        """
        returns the number of hidden layers"""
        ...
    def GetNode(self, which):
        """
        returns a particular node"""
        ...
    def GetAllNodes(self):
        """
        returns a list of all nodes"""
        ...
    lastResults: Incomplete

    def ClassifyExample(self, example, appendExamples=0):
        """
        classifies a given example and returns the results of the output layer.
        Arguments

        example: the example to be classified

        NOTE:

        if the output layer is only one element long,
        a scalar (not a list) will be returned.  This is why a lot of the other
        network code claims to only support single valued outputs."""
        ...
    def GetLastOutputs(self):
        """
        returns the complete list of output layer values from the last time this node
        classified anything"""
        ...
    def __init__(
        self,
        nodeCounts,
        nodeConnections: Incomplete | None = ...,
        actFunc=...,
        actFuncParms=...,
        weightBounds: int = ...,
    ) -> None: ...
