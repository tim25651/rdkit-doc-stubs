"""
rdkit.Chem.Suppliers.MolSupplier module¶
Supplies an abstract class for working with sequences of molecules
"""

class MolSupplier(object):
    """
    we must, at minimum, support forward iteration"""

    def NextMol(self):
        """
        Must be implemented in child class"""
        ...
    def __init__(self) -> None: ...
    def Reset(self): ...
    def __iter__(self): ...
    def next(self): ...
    def NextMol(self):
        """
        Must be implemented in child class"""
        ...
    __next__ = next
