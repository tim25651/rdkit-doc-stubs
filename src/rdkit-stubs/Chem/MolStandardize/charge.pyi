from _typeshed import Incomplete
from rdkit import Chem as Chem

from .utils import memoized_property as memoized_property

log: Incomplete

class AcidBasePair:
    name: Incomplete
    acid_str: Incomplete
    base_str: Incomplete
    def __init__(self, name, acid, base) -> None: ...
    def acid(self): ...
    def base(self): ...

ACID_BASE_PAIRS: Incomplete

class ChargeCorrection:
    name: Incomplete
    smarts_str: Incomplete
    charge: Incomplete
    def __init__(self, name, smarts, charge) -> None: ...
    def smarts(self): ...

CHARGE_CORRECTIONS: Incomplete

class Reionizer:
    acid_base_pairs: Incomplete
    charge_corrections: Incomplete
    def __init__(self, acid_base_pairs=..., charge_corrections=...) -> None: ...
    def __call__(self, mol): ...
    def reionize(self, mol): ...

class Uncharger:
    def __init__(self) -> None: ...
    def __call__(self, mol): ...
    def uncharge(self, mol): ...
