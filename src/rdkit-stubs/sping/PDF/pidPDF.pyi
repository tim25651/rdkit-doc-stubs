from _typeshed import Incomplete
from rdkit.sping import pagesizes as pagesizes
from rdkit.sping.pid import *

from . import pdfgen as pdfgen
from . import pdfgeom as pdfgeom
from . import pdfmetrics as pdfmetrics

DEFAULT_PAGE_SIZE: Incomplete
font_face_map: Incomplete
ps_font_map: Incomplete

class PDFCanvas(Canvas):
    pdf: Incomplete
    pagesize: Incomplete
    filename: Incomplete
    drawingsize: Incomplete
    pageTransitionString: str
    pageNumber: int
    def __init__(
        self, size: Incomplete | None = ..., name: str = ..., pagesize=...
    ) -> None: ...
    defaultFont: Incomplete
    defaultLineColor: Incomplete
    defaultFillColor: Incomplete
    defaultLineWidth: Incomplete
    def showPage(self) -> None: ...
    def isInteractive(self): ...
    def canUpdate(self): ...
    def clear(self) -> None: ...
    def flush(self) -> None: ...
    def save(
        self, file: Incomplete | None = ..., format: Incomplete | None = ...
    ) -> None: ...
    def setInfoLine(self, s) -> None: ...
    def __setattr__(self, key, value) -> None: ...
    def resetDefaults(self) -> None: ...
    def stringWidth(self, s, font: Incomplete | None = ...): ...
    def fontHeight(self, font: Incomplete | None = ...): ...
    def fontAscent(self, font: Incomplete | None = ...): ...
    def fontDescent(self, font: Incomplete | None = ...): ...
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
    def drawLines(
        self,
        lineList,
        color: Incomplete | None = ...,
        width: Incomplete | None = ...,
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
        fillColor: Incomplete | None = ...,
        closed: int = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawRect(
        self,
        x1,
        y1,
        x2,
        y2,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor: Incomplete | None = ...,
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
        fillColor: Incomplete | None = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawArc(
        self,
        x1,
        y1,
        x2,
        y2,
        startAng: int = ...,
        extent: int = ...,
        edgeColor: Incomplete | None = ...,
        edgeWidth: Incomplete | None = ...,
        fillColor: Incomplete | None = ...,
        dash: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...
    def drawPolygon(
        self,
        pointlist,
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
    def drawLiteral(self, literal) -> None: ...

def test(): ...
