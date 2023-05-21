"""
rdkit.ML.Cluster.Clusters module¶
contains the Cluster class for representing hierarchical cluster trees
"""
from _typeshed import Incomplete

class Cluster(object):
    """
    a class for storing clusters/data
    General Remarks

    It is assumed that the bottom of any cluster hierarchy tree is composed of
    the individual data points which were clustered.

    Clusters objects store the following pieces of data, most areaccessible via standard Setters/Getters:

    Children: Not Settable, the list of children.  You can add children
    with the _AddChild()_ and _AddChildren()_ methods.
    Note this can be of arbitrary length,
    but the current algorithms I have only produce trees with two children
    per cluster

    Metric: the metric for this cluster (i.e. how far apart its children are)
    Index: the order in which this cluster was generated

    Points: Not Settable, the list of original points in this cluster(calculated recursively from the children)

    PointsPositions: Not Settable, the list of positions of the originalpoints in this cluster (calculated recursively from the children)

    Position: the location of the cluster Note for a cluster this
    probably means the location of the average of all the Points which are
    its children.
    Data: a data field.  This is used with the original points to store their
    data value (i.e. the value we’re using to classify)
    Name: the name of this cluster

    Constructor
    Arguments

    see the class documentation for the meanings of these arguments
    my wrists are tired"""

    metric: Incomplete
    children: Incomplete
    pos: Incomplete
    index: Incomplete
    name: Incomplete
    data: Incomplete

    def __init__(
        self,
        metric: float = ...,
        children: Incomplete | None = ...,
        position: Incomplete | None = ...,
        index: int = ...,
        name: Incomplete | None = ...,
        data: Incomplete | None = ...,
    ) -> None: ...
    def SetMetric(self, metric): ...
    def GetMetric(self): ...
    def SetIndex(self, index): ...
    def GetIndex(self): ...
    def SetPosition(self, pos): ...
    def GetPosition(self): ...
    def GetPointsPositions(self): ...
    def GetPoints(self): ...
    def FindSubtree(self, index):
        """
        finds and returns the subtree with a particular index"""
        ...
    def AddChild(self, child):
        """
        Adds a child to our list
        Arguments

        child: a Cluster"""
        ...
    def AddChildren(self, children):
        """
        Adds a bunch of children to our list
        Arguments

        children: a list of Clusters"""
        ...
    def Compare(self, other, ignoreExtras=1):
        """
        not as choosy as self==other"""
        ...
    def FindSubtree(self, index):
        """
        finds and returns the subtree with a particular index"""
        ...
    def RemoveChild(self, child):
        """
        Removes a child from our list
        Arguments

        child: a Cluster"""
        ...
    def GetChildren(self): ...
    def SetData(self, data): ...
    def GetData(self): ...
    def GetIndex(self): ...
    def GetMetric(self): ...
    def SetName(self, name): ...
    def GetName(self): ...
    def GetPoints(self): ...
    def GetPointsPositions(self): ...
    def GetPosition(self): ...
    def IsTerminal(self): ...
    def Print(self, level=0, showData=0, offset="\t"): ...
    def RemoveChild(self, child):
        """
        Removes a child from our list
        Arguments

        child: a Cluster"""
        ...
    def SetData(self, data): ...
    def SetIndex(self, index): ...
    def SetMetric(self, metric): ...
    def SetName(self, name): ...
    def SetPosition(self, pos): ...
    def Compare(self, other, ignoreExtras=1):
        """
        not as choosy as self==other"""
        ...
    def IsTerminal(self): ...
    def __len__(self) -> int: ...
    def __cmp__(self, other): ...

def cmp(self, t1, t2): ...

CMPTOL: float

class Cluster(object):
    """
    a class for storing clusters/data
    General Remarks

    It is assumed that the bottom of any cluster hierarchy tree is composed of
    the individual data points which were clustered.

    Clusters objects store the following pieces of data, most areaccessible via standard Setters/Getters:

    Children: Not Settable, the list of children.  You can add children
    with the _AddChild()_ and _AddChildren()_ methods.
    Note this can be of arbitrary length,
    but the current algorithms I have only produce trees with two children
    per cluster

    Metric: the metric for this cluster (i.e. how far apart its children are)
    Index: the order in which this cluster was generated

    Points: Not Settable, the list of original points in this cluster(calculated recursively from the children)

    PointsPositions: Not Settable, the list of positions of the originalpoints in this cluster (calculated recursively from the children)

    Position: the location of the cluster Note for a cluster this
    probably means the location of the average of all the Points which are
    its children.
    Data: a data field.  This is used with the original points to store their
    data value (i.e. the value we’re using to classify)
    Name: the name of this cluster

    Constructor
    Arguments

    see the class documentation for the meanings of these arguments
    my wrists are tired"""

    metric: Incomplete
    children: Incomplete
    pos: Incomplete
    index: Incomplete
    name: Incomplete
    data: Incomplete

    def __init__(
        self,
        metric: float = ...,
        children: Incomplete | None = ...,
        position: Incomplete | None = ...,
        index: int = ...,
        name: Incomplete | None = ...,
        data: Incomplete | None = ...,
    ) -> None: ...
    def SetMetric(self, metric): ...
    def GetMetric(self): ...
    def SetIndex(self, index): ...
    def GetIndex(self): ...
    def SetPosition(self, pos): ...
    def GetPosition(self): ...
    def GetPointsPositions(self): ...
    def GetPoints(self): ...
    def FindSubtree(self, index):
        """
        finds and returns the subtree with a particular index"""
        ...
    def AddChild(self, child):
        """
        Adds a child to our list
        Arguments

        child: a Cluster"""
        ...
    def AddChildren(self, children):
        """
        Adds a bunch of children to our list
        Arguments

        children: a list of Clusters"""
        ...
    def Compare(self, other, ignoreExtras=1):
        """
        not as choosy as self==other"""
        ...
    def FindSubtree(self, index):
        """
        finds and returns the subtree with a particular index"""
        ...
    def RemoveChild(self, child):
        """
        Removes a child from our list
        Arguments

        child: a Cluster"""
        ...
    def GetChildren(self): ...
    def SetData(self, data): ...
    def GetData(self): ...
    def GetIndex(self): ...
    def GetMetric(self): ...
    def SetName(self, name): ...
    def GetName(self): ...
    def GetPoints(self): ...
    def GetPointsPositions(self): ...
    def GetPosition(self): ...
    def IsTerminal(self): ...
    def Print(self, level=0, showData=0, offset="\t"): ...
    def RemoveChild(self, child):
        """
        Removes a child from our list
        Arguments

        child: a Cluster"""
        ...
    def SetData(self, data): ...
    def SetIndex(self, index): ...
    def SetMetric(self, metric): ...
    def SetName(self, name): ...
    def SetPosition(self, pos): ...
    def Compare(self, other, ignoreExtras=1):
        """
        not as choosy as self==other"""
        ...
    def IsTerminal(self): ...
    def __len__(self) -> int: ...
    def __cmp__(self, other): ...
