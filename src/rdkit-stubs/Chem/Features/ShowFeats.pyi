"""
rdkit.Chem.Features.ShowFeats module¶
"""
from _typeshed import Incomplete
from rdkit import Geometry as Geometry
from rdkit import RDConfig as RDConfig

logger: Incomplete
BEGIN: int
END: int
TRIANGLE_FAN: int
COLOR: int
VERTEX: int
NORMAL: int
SPHERE: int
CYLINDER: int
ALPHA: int

def ShowArrow(
    self,
    viewer,
    tail,
    head,
    radius,
    color,
    label,
    transparency=0,
    includeArrowhead=True,
): ...
def ShowMolFeats(
    self,
    mol,
    factory,
    viewer,
    radius=0.5,
    confId=-1,
    showOnly=True,
    name="",
    transparency=0.0,
    colors=None,
    excludeTypes=[],
    useFeatDirs=True,
    featLabel=None,
    dirLabel=None,
    includeArrowheads=True,
    writeFeats=False,
    showMol=True,
    featMapFile=False,
): ...

parser: Incomplete
