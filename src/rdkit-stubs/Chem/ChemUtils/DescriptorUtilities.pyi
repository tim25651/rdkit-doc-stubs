"""
rdkit.Chem.ChemUtils.DescriptorUtilities module¶
Collection of utilities to be used with descriptors
"""
from _typeshed import Incomplete

def setDescriptorVersion(self, version="1.0.0"):
    """
    Set the version on the descriptor function.
    Use as a decorator"""
    ...

class VectorDescriptorNamespace(dict):
    ...

    def __init__(self, **kwargs) -> None: ...

class VectorDescriptorWrapper(object):
    """
    Wrap a function that returns a vector and make it seem like there
    is one function for each entry.  These functions are added to the global
    namespace with the names provided"""

    func: Incomplete
    names: Incomplete
    func_key: Incomplete
    namespace: Incomplete

    def __init__(self, func, names, version, namespace) -> None: ...
    def call_desc(self, mol, index): ...

def setDescriptorVersion(self, version="1.0.0"):
    """
    Set the version on the descriptor function.
    Use as a decorator"""
    ...
