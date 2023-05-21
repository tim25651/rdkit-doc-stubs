"""
rdkit.Dbase.DbReport module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import rdDepictor as rdDepictor
from rdkit.Chem.Draw import DrawUtils as DrawUtils
from rdkit.Chem.Draw.MolDrawing import DrawingOptions as DrawingOptions
from rdkit.Chem.Draw.MolDrawing import MolDrawing as MolDrawing
from rdkit.Dbase import DbInfo as DbInfo
from rdkit.Dbase.DbConnection import DbConnect as DbConnect
from rdkit.Reports.PDFReport import PDFReport as PDFReport
from rdkit.Reports.PDFReport import ReportUtils as ReportUtils
from rdkit.utils import cactvs as cactvs

GetReportlabTable: Incomplete
QuickReport: Incomplete
hasCDX: int

class CDXImageTransformer:
    smiCol: Incomplete
    tempHandler: Incomplete
    width: Incomplete
    verbose: Incomplete

    def __init__(
        self,
        smiCol,
        width: int = ...,
        verbose: int = ...,
        tempHandler: Incomplete | None = ...,
    ) -> None: ...
    def __call__(self, arg): ...

class CactvsImageTransformer:
    smiCol: Incomplete
    tempHandler: Incomplete
    width: Incomplete
    verbose: Incomplete

    def __init__(
        self,
        smiCol,
        width: float = ...,
        verbose: int = ...,
        tempHandler: Incomplete | None = ...,
    ) -> None: ...
    def __call__(self, arg): ...

class ReportLabImageTransformer:
    smiCol: Incomplete
    width: Incomplete
    verbose: Incomplete

    def __init__(
        self,
        smiCol,
        width: float = ...,
        verbose: int = ...,
        tempHandler: Incomplete | None = ...,
    ) -> None: ...
    def __call__(self, arg): ...

class RDImageTransformer:
    smiCol: Incomplete
    tempHandler: Incomplete
    width: Incomplete
    verbose: Incomplete

    def __init__(
        self,
        smiCol,
        width: float = ...,
        verbose: int = ...,
        tempHandler: Incomplete | None = ...,
    ) -> None: ...
    def __call__(self, arg): ...
