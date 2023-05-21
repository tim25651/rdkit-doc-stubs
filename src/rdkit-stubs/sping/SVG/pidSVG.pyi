from math import *

from _typeshed import Incomplete
from rdkit.sping.PDF import pdfmetrics as pdfmetrics
from rdkit.sping.pid import *

SVG_HEADER: str

class SVGCanvas(Canvas):
    size: Incomplete
    def __init__(
        self,
        size=...,
        name: str = ...,
        includeXMLHeader: bool = ...,
        extraHeaderText: str = ...,
    ) -> None: ...
    def clear(self) -> None: ...
    def flush(self) -> None: ...
    def save(
        self, file: Incomplete | None = ..., format: Incomplete | None = ...
    ) -> None: ...
    def text(self): ...
    def drawLine(
        self,
        x1,
        y1,
        x2,
        y2,
        color: Incomplete | None = ...,
        width: Incomplete | None = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawPolygon(
        self,
        pointlist,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor=...,
        closed: int = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawEllipse(
        self,
        x1,
        y1,
        x2,
        y2,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor=...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawArc(
        self,
        x1,
        y1,
        x2,
        y2,
        theta1: int = ...,
        extent: int = ...,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor: Incomplete | None = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawCurve(
        self,
        x1,
        y1,
        x2,
        y2,
        x3,
        y3,
        x4,
        y4,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor=...,
        closed: int = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawString(
        self,
        s,
        x,
        y,
        font: Incomplete | None = ...,
        color: Incomplete | None = ...,
        angle: int = ...,
        **kwargs,
    ) -> None: ...
    def drawFigure(
        self,
        partList,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor: Incomplete | None = ...,
        closed: int = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawImage(
        self,
        image,
        x1,
        y1,
        x2: Incomplete | None = ...,
        y2: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def stringWidth(self, s, font: Incomplete | None = ...): ...
    def fontAscent(self, font: Incomplete | None = ...): ...
    def fontDescent(self, font: Incomplete | None = ...): ...

def test(): ...
def dashtest(): ...
def testit(canvas, s, x, y, font: Incomplete | None = ...) -> None: ...
def test2() -> None: ...
