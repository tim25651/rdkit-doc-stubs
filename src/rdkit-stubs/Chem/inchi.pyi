"""
rdkit.Chem.inchi module¶
"""
from _typeshed import Incomplete

INCHI_AVAILABLE: bool

class InchiReadWriteError(Exception):
    ...
    ...

def InchiToInchiKey(self, inchi):
    """
    Return the InChI key for the given InChI string. Return None on error"""
    ...

def MolBlockToInchi(
    self, molblock, options="", logLevel=None, treatWarningAsError=False
):
    """
    Returns the standard InChI string for a mol block
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    the standard InChI string returned by InChI API for the input molecule"""
    ...

def MolBlockToInchiAndAuxInfo(
    self, molblock, options="", logLevel=None, treatWarningAsError=False
):
    """
    Returns the standard InChI string and InChI auxInfo for a mol block
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    a tuple of the standard InChI string and the auxInfo string returned by
    InChI API, in that order, for the input molecule"""
    ...

def MolFromInchi(
    self, inchi, sanitize=True, removeHs=True, logLevel=None, treatWarningAsError=False
):
    """
    Construct a molecule from a InChI string
    Keyword arguments:
    sanitize – set to True to enable sanitization of the molecule. Default is
    True
    removeHs – set to True to remove Hydrogens from a molecule. This only
    makes sense when sanitization is enabled
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant
    molecule  and error message are part of the excpetion
    Returns:
    a rdkit.Chem.rdchem.Mol instance"""
    ...

def MolToInchi(self, mol, options="", logLevel=None, treatWarningAsError=False):
    """
    Returns the standard InChI string for a molecule
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    the standard InChI string returned by InChI API for the input molecule"""
    ...

def MolToInchiAndAuxInfo(
    self, mol, options="", logLevel=None, treatWarningAsError=False
):
    """
    Returns the standard InChI string and InChI auxInfo for a molecule
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    a tuple of the standard InChI string and the auxInfo string returned by
    InChI API, in that order, for the input molecule"""
    ...

def MolBlockToInchiAndAuxInfo(
    self, molblock, options="", logLevel=None, treatWarningAsError=False
):
    """
    Returns the standard InChI string and InChI auxInfo for a mol block
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    a tuple of the standard InChI string and the auxInfo string returned by
    InChI API, in that order, for the input molecule"""
    ...

def MolToInchi(self, mol, options="", logLevel=None, treatWarningAsError=False):
    """
    Returns the standard InChI string for a molecule
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    the standard InChI string returned by InChI API for the input molecule"""
    ...

def MolBlockToInchi(
    self, molblock, options="", logLevel=None, treatWarningAsError=False
):
    """
    Returns the standard InChI string for a mol block
    Keyword arguments:
    logLevel – the log level used for logging logs and messages from InChI
    API. set to None to diable the logging completely
    treatWarningAsError – set to True to raise an exception in case of a
    molecule that generates warning in calling InChI API. The resultant InChI
    string and AuxInfo string as well as the error message are encoded in the
    exception.
    Returns:
    the standard InChI string returned by InChI API for the input molecule"""
    ...

def InchiToInchiKey(self, inchi):
    """
    Return the InChI key for the given InChI string. Return None on error"""
    ...

def MolToInchiKey(self, mol, options=""):
    """
    Returns the standard InChI key for a molecule
    Returns:
    the standard InChI key returned by InChI API for the input molecule"""
    ...
