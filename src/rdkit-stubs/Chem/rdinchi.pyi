"""
rdkit.Chem.rdinchi module¶
"""
from typing import Any

from rdkit.Chem.rdchem import Mol

def InchiToInchiKey(self, inchi: str) -> str:
    """
    return the InChI key for an InChI string

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > InchiToInchiKey(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """
    ...

def InchiToMol(self, inchi: str, sanitize: bool = True, removeHs: bool = True) -> tuple:
    """
    return a ROMol for a InChI string
    Returns:
    a tuple with:
    the molecule
    the return code from the InChI conversion
    a string with any messages from the InChI conversion
    a string with any log messages from the InChI conversion

    C++ signature :boost::python::tuple InchiToMol(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,bool=True [,bool=True]])
    """
    ...

def MolBlockToInchi(self, molblock: str, options: str = "") -> tuple:
    """
    return the InChI for a ROMol molecule.

    Arguments:
    molblock: the mol block to use.
    options: the InChI generation options.
    Options should be prefixed with either a - or a /
    Available options are explained in the InChI technical FAQ:
    http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14
    and the User Guide:
    http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf

    Returns:
    a tuple with:
    the InChI
    the return code from the InChI conversion
    a string with any messages from the InChI conversion
    a string with any log messages from the InChI conversion
    a string with the InChI AuxInfo

    C++ signature :boost::python::tuple MolBlockToInchi(std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’])
    """
    ...

def MolToInchi(self, mol: Mol, options: str = "") -> tuple:
    """
    return the InChI for a ROMol molecule.

    Arguments:
    mol: the molecule to use.
    options: the InChI generation options.
    Options should be prefixed with either a - or a /
    Available options are explained in the InChI technical FAQ:
    http://www.inchi-trust.org/fileadmin/user_upload/html/inchifaq/inchi-faq.html#15.14
    and the User Guide:
    http://www.inchi-trust.org/fileadmin/user_upload/software/inchi-v1.04/InChI_UserGuide.pdf

    Returns:
    a tuple with:
    the InChI
    the return code from the InChI conversion
    a string with any messages from the InChI conversion
    a string with any log messages from the InChI conversion
    a string with the InChI AuxInfo

    C++ signature :boost::python::tuple MolToInchi(RDKit::ROMol [,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >=’’])
    """
    ...

def MolToInchiKey(self, mol: Mol, options: str = "") -> str:
    """
    return the InChI key for a ROMol molecule.

    Arguments:
    mol: the molecule to use.
    options: the InChI generation options.
    Options should be prefixed with either a - or a /
    Available options are explained in the InChI technical FAQ:
    https://www.inchi-trust.org/technical-faq-2/#15.14
    and the User Guide available from:
    https://www.inchi-trust.org/downloads/

    Returns: the InChI key

    C++ signature :std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> > MolToInchiKey(RDKit::ROMol [,char const*=’’])
    """
    ...
