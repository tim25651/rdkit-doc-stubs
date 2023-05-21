"""
rdkit.Chem.Draw.qtCanvas module¶
"""
from _typeshed import Incomplete
from rdkit.Chem.Draw import canvasbase
from rdkit.Chem.Draw.canvasbase import CanvasBase as CanvasBase
from rdkit.Chem.Draw.rdMolDraw2DQt import rdkitQtVersion as rdkitQtVersion

QPainter_Antialiasing: Incomplete
QPainter_SmoothPixmapTransform: Incomplete
GlobalColor_white: Incomplete
PenStile_SolidLine: Incomplete
PenStile_DashLine: Incomplete

class Canvas(canvasbase.CanvasBase):
    def addCanvasDashedWedge(
        self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs
    ):
        """
        Draw a dashed wedge
        The wedge is identified by the three points p1, p2, and p3.
        It will be drawn using the given color; if color2 is specified
        it will be used for the second half of the wedge
        TODO: fix comment, I’m not sure what dash does"""
        ...
    size: Incomplete
    qsize: Incomplete
    pixmap: Incomplete
    painter: Incomplete

    def __init__(self, size) -> None: ...
    def addCanvasLine(self, p1, p2, color=(0, 0, 0), color2=None, **kwargs):
        """
        Draw a single line on the canvas
        This function will draw a line between p1 and p2 with the
        given color.
        If color2 is specified, it will be used to draw the second half
        of the segment"""
        ...
    def addCanvasPolygon(self, ps, color=(0, 0, 0), fill=True, stroke=False, **kwargs):
        """
        Draw a polygon
        Draw a polygon identified by vertexes given in ps using
        the given color"""
        ...
    def addCanvasText(self, text, pos, font, color=(0, 0, 0), **kwargs):
        """
        Draw some text
        The provided text is drawn at position pos using the given
        font and the chosen color."""
        ...
    def addCanvasPolygon(self, ps, color=(0, 0, 0), fill=True, stroke=False, **kwargs):
        """
        Draw a polygon
        Draw a polygon identified by vertexes given in ps using
        the given color"""
        ...
    def addCanvasDashedWedge(
        self, p1, p2, p3, dash=(2, 2), color=(0, 0, 0), color2=None, **kwargs
    ):
        """
        Draw a dashed wedge
        The wedge is identified by the three points p1, p2, and p3.
        It will be drawn using the given color; if color2 is specified
        it will be used for the second half of the wedge
        TODO: fix comment, I’m not sure what dash does"""
        ...
    def flush(self):
        """
        Complete any remaining draw operation
        This is supposed to be the last operation on the canvas before
        saving it"""
        ...
