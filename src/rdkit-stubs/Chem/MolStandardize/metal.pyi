from _typeshed import Incomplete
from rdkit import Chem as Chem

log: Incomplete

class MetalDisconnector:
    def __init__(self) -> None: ...
    def __call__(self, mol): ...
    def disconnect(self, mol): ...