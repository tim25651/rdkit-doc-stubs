from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors

from .utils import memoized_property as memoized_property

log: Incomplete

class FragmentPattern:
    name: Incomplete
    smarts_str: Incomplete
    def __init__(self, name, smarts) -> None: ...
    def smarts(self): ...

REMOVE_FRAGMENTS: Incomplete
LEAVE_LAST: bool
PREFER_ORGANIC: bool

def is_organic(fragment): ...

class FragmentRemover:
    fragments: Incomplete
    leave_last: Incomplete
    def __init__(self, fragments=..., leave_last=...) -> None: ...
    def __call__(self, mol): ...
    def remove(self, mol): ...

class LargestFragmentChooser:
    prefer_organic: Incomplete
    def __init__(self, prefer_organic=...) -> None: ...
    def __call__(self, mol): ...
    def choose(self, mol): ...