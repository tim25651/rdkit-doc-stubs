from _typeshed import Incomplete
from rdkit.sping.pid import *

from . import psmetrics as psmetrics

class PostScriptLevelException(ValueError): ...

linesep: str
PiddleLegalFonts: Incomplete
Roman: str
Bold: str
Italic: str
PSFontMapStdEnc: Incomplete
PSFontMapLatin1Enc: Incomplete

def latin1FontEncoding(fontname): ...
def dashLineDefinition(): ...

class PsDSC:
    def __init__(self) -> None: ...
    def documentHeader(self): ...
    def boundingBoxStr(self, x0, y0, x1, y1): ...
    inPageFlag: int
    def BeginPageStr(self, pageSetupStr, pageName: Incomplete | None = ...): ...
    def EndPageStr(self): ...

class EpsDSC(PsDSC):
    def __init__(self) -> None: ...
    def documentHeader(self): ...

class PSCanvas(Canvas):
    filename: Incomplete
    drawImage: Incomplete
    code: Incomplete
    dsc: Incomplete
    defaultFont: Incomplete
    fontMapEncoding: Incomplete
    pageNum: int
    def __init__(
        self, size=..., name: str = ..., PostScriptLevel: int = ..., fontMapEncoding=...
    ) -> None: ...
    def psBeginDocument(self) -> None: ...
    def psEndDocument(self) -> None: ...
    def psBeginPage(self, pageName: Incomplete | None = ...) -> None: ...
    def psEndPage(self) -> None: ...
    def nextPage(self) -> None: ...
    def clear(self) -> None: ...
    def resetToDefaults(self) -> None: ...
    def flush(self) -> None: ...
    def save(
        self, file: Incomplete | None = ..., format: Incomplete | None = ...
    ) -> None: ...
    def stringWidth(self, s, font: Incomplete | None = ...): ...
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
    def drawRoundRect(
        self,
        x1,
        y1,
        x2,
        y2,
        rx: int = ...,
        ry: int = ...,
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
