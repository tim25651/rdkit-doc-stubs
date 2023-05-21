"""
rdkit.ML.Cluster.ClusterVis module¶
Cluster tree visualization using Sping
"""
from _typeshed import Incomplete
from rdkit.piddle import piddle as piddle
from rdkit.sping import pid as pid
from rdkit.sping.colors import Color

from . import ClusterUtils as ClusterUtils

class ClusterRenderer(object):
    canvas: Incomplete
    size: Incomplete
    ptColors: Incomplete
    lineWidth: Incomplete
    showIndices: Incomplete
    showNodes: Incomplete
    stopAtCentroids: Incomplete
    logScale: Incomplete
    tooClose: Incomplete

    def __init__(
        self,
        canvas,
        size,
        ptColors=...,
        lineWidth: Incomplete | None = ...,
        showIndices: int = ...,
        showNodes: int = ...,
        stopAtCentroids: int = ...,
        logScale: int = ...,
        tooClose: int = ...,
    ) -> None: ...
    ySpace: Incomplete

    def DrawTree(self, cluster, minHeight=2.0): ...

piddle = pid

class VisOpts(object):
    """
    stores visualization options for cluster viewing
    Instance variables

    x/yOffset: amount by which the drawing is offset from the edges of the canvas
    lineColor: default color for drawing the cluster tree
    lineWidth: the width of the lines used to draw the tree"""

    xOffset: int = ...
    yOffset: int = ...
    lineColor: Color = ...
    hideColor: Color = ...
    terminalColors: list[Color] = ...
    lineWidth: int = ...
    hideWidth: float = ...
    nodeRad: int = ...
    nodeColor: Color = ...
    highlightColor: Color = ...
    highlightRad: int = ...
    lineColor: Color = ...
    lineWidth: int = ...
    nodeColor: Color = ...
    nodeRad: int = ...
    terminalColors: list[Color] = ...
    xOffset: int = ...
    yOffset: int = ...

def ClusterToImg(
    self,
    cluster,
    fileName,
    size=(300, 300),
    ptColors=[],
    lineWidth=None,
    showIndices=0,
    stopAtCentroids=0,
    logScale=0,
):
    """
    handles the work of drawing a cluster tree to an image file

    Arguments

    cluster: the cluster tree to be drawn
    fileName: the name of the file to be created
    size: the size of output canvas
    ptColors: if this is specified, the _colors_ will be used to color
    the terminal nodes of the cluster tree.  (color == _pid.Color_)
    lineWidth: if specified, it will be used for the widths of the lines
    used to draw the tree

    Notes

    The extension on  _fileName_ determines the type of image file created.
    All formats supported by PIL can be used.
    if _ptColors_ is the wrong length for the number of possible terminal
    node types, this will throw an IndexError
    terminal node types are determined using their _GetData()_ methods"""
    ...

class ClusterRenderer(object):
    canvas: Incomplete
    size: Incomplete
    ptColors: Incomplete
    lineWidth: Incomplete
    showIndices: Incomplete
    showNodes: Incomplete
    stopAtCentroids: Incomplete
    logScale: Incomplete
    tooClose: Incomplete

    def __init__(
        self,
        canvas,
        size,
        ptColors=...,
        lineWidth: Incomplete | None = ...,
        showIndices: int = ...,
        showNodes: int = ...,
        stopAtCentroids: int = ...,
        logScale: int = ...,
        tooClose: int = ...,
    ) -> None: ...
    ySpace: Incomplete

    def DrawTree(self, cluster, minHeight=2.0): ...

def DrawClusterTree(
    self,
    cluster,
    canvas,
    size,
    ptColors=[],
    lineWidth=None,
    showIndices=0,
    showNodes=1,
    stopAtCentroids=0,
    logScale=0,
    tooClose=-1,
):
    """
    handles the work of drawing a cluster tree on a Sping canvas

    Arguments

    cluster: the cluster tree to be drawn
    canvas:  the Sping canvas on which to draw
    size: the size of _canvas_
    ptColors: if this is specified, the _colors_ will be used to color
    the terminal nodes of the cluster tree.  (color == _pid.Color_)
    lineWidth: if specified, it will be used for the widths of the lines
    used to draw the tree

    Notes

    _Canvas_ is neither _save_d nor _flush_ed at the end of this
    if _ptColors_ is the wrong length for the number of possible terminal
    node types, this will throw an IndexError
    terminal node types are determined using their _GetData()_ methods"""
    ...

def ClusterToPDF(
    self,
    cluster,
    fileName,
    size=(300, 300),
    ptColors=[],
    lineWidth=None,
    showIndices=0,
    stopAtCentroids=0,
    logScale=0,
):
    """
    handles the work of drawing a cluster tree to an PDF file

    Arguments

    cluster: the cluster tree to be drawn
    fileName: the name of the file to be created
    size: the size of output canvas
    ptColors: if this is specified, the _colors_ will be used to color
    the terminal nodes of the cluster tree.  (color == _pid.Color_)
    lineWidth: if specified, it will be used for the widths of the lines
    used to draw the tree

    Notes

    if _ptColors_ is the wrong length for the number of possible terminal
    node types, this will throw an IndexError
    terminal node types are determined using their _GetData()_ methods"""
    ...

def ClusterToSVG(
    self,
    cluster,
    fileName,
    size=(300, 300),
    ptColors=[],
    lineWidth=None,
    showIndices=0,
    stopAtCentroids=0,
    logScale=0,
):
    """
    handles the work of drawing a cluster tree to an SVG file

    Arguments

    cluster: the cluster tree to be drawn
    fileName: the name of the file to be created
    size: the size of output canvas
    ptColors: if this is specified, the _colors_ will be used to color
    the terminal nodes of the cluster tree.  (color == _pid.Color_)
    lineWidth: if specified, it will be used for the widths of the lines
    used to draw the tree

    Notes

    if _ptColors_ is the wrong length for the number of possible terminal
    node types, this will throw an IndexError
    terminal node types are determined using their _GetData()_ methods"""
    ...

def DrawClusterTree(
    self,
    cluster,
    canvas,
    size,
    ptColors=[],
    lineWidth=None,
    showIndices=0,
    showNodes=1,
    stopAtCentroids=0,
    logScale=0,
    tooClose=-1,
):
    """
    handles the work of drawing a cluster tree on a Sping canvas

    Arguments

    cluster: the cluster tree to be drawn
    canvas:  the Sping canvas on which to draw
    size: the size of _canvas_
    ptColors: if this is specified, the _colors_ will be used to color
    the terminal nodes of the cluster tree.  (color == _pid.Color_)
    lineWidth: if specified, it will be used for the widths of the lines
    used to draw the tree

    Notes

    _Canvas_ is neither _save_d nor _flush_ed at the end of this
    if _ptColors_ is the wrong length for the number of possible terminal
    node types, this will throw an IndexError
    terminal node types are determined using their _GetData()_ methods"""
    ...

def ClusterToImg(
    self,
    cluster,
    fileName,
    size=(300, 300),
    ptColors=[],
    lineWidth=None,
    showIndices=0,
    stopAtCentroids=0,
    logScale=0,
):
    """
    handles the work of drawing a cluster tree to an image file

    Arguments

    cluster: the cluster tree to be drawn
    fileName: the name of the file to be created
    size: the size of output canvas
    ptColors: if this is specified, the _colors_ will be used to color
    the terminal nodes of the cluster tree.  (color == _pid.Color_)
    lineWidth: if specified, it will be used for the widths of the lines
    used to draw the tree

    Notes

    The extension on  _fileName_ determines the type of image file created.
    All formats supported by PIL can be used.
    if _ptColors_ is the wrong length for the number of possible terminal
    node types, this will throw an IndexError
    terminal node types are determined using their _GetData()_ methods"""
    ...
