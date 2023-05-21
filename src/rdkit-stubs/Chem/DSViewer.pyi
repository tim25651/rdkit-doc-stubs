from _typeshed import Incomplete
from rdkit import Chem as Chem

class Displayable:
    doc: Incomplete
    id: Incomplete
    visible: bool
    children: Incomplete
    def __init__(self, doc, id: int = ...) -> None: ...
    def Select(self, atoms=..., state: bool = ..., recurse: bool = ...): ...
    def Hide(self, recurse: bool = ...) -> None: ...
    def Show(self, recurse: bool = ...) -> None: ...
    def ShowOnly(self, recurse: bool = ...) -> None: ...
    def __del__(self) -> None: ...

class MolViewer:
    app: Incomplete
    doc: Incomplete
    displayables: Incomplete
    def __init__(self, force: int = ..., title: str = ..., **kwargs) -> None: ...
    def DeleteAll(self) -> None: ...
    def DeleteAllExcept(self, excludes) -> None: ...
    def ShowMol(
        self,
        mol,
        name: str = ...,
        showOnly: bool = ...,
        highlightFeatures=...,
        molB: str = ...,
        confId: int = ...,
        zoom: bool = ...,
    ) -> None: ...
    def LoadFile(self, filename, name, showOnly: bool = ...): ...
    def GetSelectedAtoms(self, whichSelection: str = ...): ...
    def HighlightAtoms(self, indices, where, extraHighlight: bool = ...) -> None: ...
    def SelectAtoms(self, itemId, atomIndices, selName: str = ...) -> None: ...
    def SetDisplayUpdate(self, val) -> None: ...
    def GetAtomCoords(self, sels): ...
    def AddPharmacophore(self, locs, colors, label, sphereRad: float = ...) -> None: ...
    def SetDisplayStyle(self, obj, style: str = ...) -> None: ...
    def HideAll(self) -> None: ...
    def HideObject(self, objName) -> None: ...
    def DisplayObject(self, objName) -> None: ...
    def Zoom(self, objName) -> None: ...
    def SelectProteinNeighborhood(
        self,
        aroundObj,
        inObj,
        distance: float = ...,
        name: str = ...,
        showSurface: bool = ...,
    ) -> None: ...
    def Redraw(self) -> None: ...