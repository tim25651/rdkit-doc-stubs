from _typeshed import Incomplete

StandardEnglishFonts: Incomplete
widths: Incomplete
ascent_descent: Incomplete

def parseAFMfile(filename): ...

class FontCache:
    def __init__(self) -> None: ...
    def loadfont(self, fontname) -> None: ...
    def getfont(self, fontname): ...
    def stringwidth(self, text, font): ...
    def status(self): ...

TheFontCache: Incomplete
stringwidth: Incomplete