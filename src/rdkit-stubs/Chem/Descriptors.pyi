"""
rdkit.Chem.Descriptors module¶
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem import rdPartialCharges as rdPartialCharges
from rdkit.Chem.EState.EState import MaxAbsEStateIndex as MaxAbsEStateIndex
from rdkit.Chem.EState.EState import MaxEStateIndex as MaxEStateIndex
from rdkit.Chem.EState.EState import MinAbsEStateIndex as MinAbsEStateIndex
from rdkit.Chem.EState.EState import MinEStateIndex as MinEStateIndex
from rdkit.Chem.QED import qed as qed

class PropertyFunctor(rdMolDescriptors.PythonPropertyFunctor):
    """
    Creates a python based property function that can be added to the
    global property list.  To use, subclass this class and override the
    __call__ method.  Then create an instance and add it to the
    registry.  The __call__ method should return a numeric value.
    Example:

    class NumAtoms(Descriptors.PropertyFunctor):
    def __init__(self):Descriptors.PropertyFunctor.__init__(self, “NumAtoms”, “1.0.0”)

    def __call__(self, mol):return mol.GetNumAtoms()

    numAtoms = NumAtoms()
    rdMolDescriptors.Properties.RegisterProperty(numAtoms)

    C++ signature :void __init__(_object*,_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    ...

    def __init__(self, name, version) -> None: ...
    def __call__(self, mol) -> None: ...

def CalcMolDescriptors(self, mol, missingVal=None, silent=True):
    """
    calculate the full set of descriptors for a molecule
    mol : RDKit molecule
    missingVal : float, optional

    This will be used if a particular descriptor cannot be calculated

    silentbool, optionalif True then exception messages from descriptors will be displayed

    dict A dictionary with decriptor names as keys and the descriptor values as values
    """
    ...

def MolWt(self, *x, **y):
    """
    The average molecular weight of the molecule
    >>> MolWt(Chem.MolFromSmiles('CC'))
    30.07
    >>> MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
    53.49..."""
    ...

def HeavyAtomMolWt(self, x):
    """
    The average molecular weight of the molecule ignoring hydrogens
    >>> HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
    24.02...
    >>> HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
    49.46"""
    ...

def ExactMolWt(self, *x, **y):
    """
    The exact molecular weight of the molecule
    >>> ExactMolWt(Chem.MolFromSmiles('CC'))
    30.04...
    >>> ExactMolWt(Chem.MolFromSmiles('[13CH3]C'))
    31.05..."""
    ...

def NumValenceElectrons(self, mol):
    """
    The number of valence electrons the molecule has
    >>> NumValenceElectrons(Chem.MolFromSmiles('CC'))
    14
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)O'))
    18
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)[O-]'))
    18
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)'))
    12"""
    ...

def NumRadicalElectrons(self, mol):
    """
    The number of radical electrons the molecule has(says nothing about spin state)

    >>> NumRadicalElectrons(Chem.MolFromSmiles('CC'))
    0
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH3]'))
    0
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH2]'))
    1
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH]'))
    2
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[C]'))
    3"""
    ...

def MaxPartialCharge(self, mol, force=False): ...
def MinPartialCharge(self, mol, force=False): ...
def MaxAbsPartialCharge(self, mol, force=False): ...
def MinAbsPartialCharge(self, mol, force=False): ...
def FpDensityMorgan1(self, x): ...
def FpDensityMorgan2(self, x): ...
def FpDensityMorgan3(self, x): ...
def HeavyAtomMolWt(self, x):
    """
    The average molecular weight of the molecule ignoring hydrogens
    >>> HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
    24.02...
    >>> HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
    49.46"""
    ...

def MaxAbsPartialCharge(self, mol, force=False): ...
def MaxPartialCharge(self, mol, force=False): ...
def MinAbsPartialCharge(self, mol, force=False): ...
def MinPartialCharge(self, mol, force=False): ...
def MolWt(self, *x, **y):
    """
    The average molecular weight of the molecule
    >>> MolWt(Chem.MolFromSmiles('CC'))
    30.07
    >>> MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
    53.49..."""
    ...

def NumRadicalElectrons(self, mol):
    """
    The number of radical electrons the molecule has(says nothing about spin state)

    >>> NumRadicalElectrons(Chem.MolFromSmiles('CC'))
    0
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH3]'))
    0
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH2]'))
    1
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH]'))
    2
    >>> NumRadicalElectrons(Chem.MolFromSmiles('C[C]'))
    3"""
    ...

def NumValenceElectrons(self, mol):
    """
    The number of valence electrons the molecule has
    >>> NumValenceElectrons(Chem.MolFromSmiles('CC'))
    14
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)O'))
    18
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)[O-]'))
    18
    >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)'))
    12"""
    ...

names: Incomplete
autocorr: Incomplete

def setupAUTOCorrDescriptors(self):
    """
    Adds AUTOCORR descriptors to the default descriptor lists"""
    ...

class PropertyFunctor(rdMolDescriptors.PythonPropertyFunctor):
    """
    Creates a python based property function that can be added to the
    global property list.  To use, subclass this class and override the
    __call__ method.  Then create an instance and add it to the
    registry.  The __call__ method should return a numeric value.
    Example:

    class NumAtoms(Descriptors.PropertyFunctor):
    def __init__(self):Descriptors.PropertyFunctor.__init__(self, “NumAtoms”, “1.0.0”)

    def __call__(self, mol):return mol.GetNumAtoms()

    numAtoms = NumAtoms()
    rdMolDescriptors.Properties.RegisterProperty(numAtoms)

    C++ signature :void __init__(_object*,_object*,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >,std::__cxx11::basic_string<char, std::char_traits<char>, std::allocator<char> >)
    """

    ...

    def __init__(self, name, version) -> None: ...
    def __call__(self, mol) -> None: ...

def CalcMolDescriptors(self, mol, missingVal=None, silent=True):
    """
    calculate the full set of descriptors for a molecule
    mol : RDKit molecule
    missingVal : float, optional

    This will be used if a particular descriptor cannot be calculated

    silentbool, optionalif True then exception messages from descriptors will be displayed

    dict A dictionary with decriptor names as keys and the descriptor values as values
    """
    ...
