from _typeshed import Incomplete
from rdkit.sping.pid import *

f: Incomplete

class PILCanvas(Canvas):
    def __init__(self, size=..., name: str = ...) -> None: ...
    def __setattr__(self, attribute, value) -> None: ...
    def getImage(self): ...
    def save(
        self, file: Incomplete | None = ..., format: Incomplete | None = ...
    ) -> None: ...
    def clear(self) -> None: ...
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
    def drawString(
        self,
        s,
        x,
        y,
        font: Incomplete | None = ...,
        color: Incomplete | None = ...,
        angle: int = ...,
        **kwargs,
    ): ...
    def drawImage(
        self,
        image,
        x1,
        y1,
        x2: Incomplete | None = ...,
        y2: Incomplete | None = ...,
        **kwargs,
    ) -> None: ...

def test(): ...
def testit(canvas, s, x, y, font: Incomplete | None = ...) -> None: ...
def test2() -> None: ...
