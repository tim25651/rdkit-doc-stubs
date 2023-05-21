"""
rdkit.Chem.MolSurf module¶
Exposes functionality for MOE-like approximate molecular surface area
descriptors.

The MOE-like VSA descriptors are also calculated here
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit.Chem import Crippen as Crippen
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors
from rdkit.Chem import rdPartialCharges as rdPartialCharges

def LabuteASA(self, *x, **y): ...
def PEOE_VSA1(self, x, y=0):
    """
    MOE Charge VSA Descriptor 1 (-inf < x < -0.30)"""
    ...

def PEOE_VSA10(self, x, y=9):
    """
    MOE Charge VSA Descriptor 10 ( 0.10 <= x <  0.15)"""
    ...

def PEOE_VSA11(self, x, y=10):
    """
    MOE Charge VSA Descriptor 11 ( 0.15 <= x <  0.20)"""
    ...

def PEOE_VSA12(self, x, y=11):
    """
    MOE Charge VSA Descriptor 12 ( 0.20 <= x <  0.25)"""
    ...

def PEOE_VSA13(self, x, y=12):
    """
    MOE Charge VSA Descriptor 13 ( 0.25 <= x <  0.30)"""
    ...

def PEOE_VSA14(self, x, y=13):
    """
    MOE Charge VSA Descriptor 14 ( 0.30 <= x < inf)"""
    ...

def PEOE_VSA2(self, x, y=1):
    """
    MOE Charge VSA Descriptor 2 (-0.30 <= x < -0.25)"""
    ...

def PEOE_VSA3(self, x, y=2):
    """
    MOE Charge VSA Descriptor 3 (-0.25 <= x < -0.20)"""
    ...

def PEOE_VSA4(self, x, y=3):
    """
    MOE Charge VSA Descriptor 4 (-0.20 <= x < -0.15)"""
    ...

def PEOE_VSA5(self, x, y=4):
    """
    MOE Charge VSA Descriptor 5 (-0.15 <= x < -0.10)"""
    ...

def PEOE_VSA6(self, x, y=5):
    """
    MOE Charge VSA Descriptor 6 (-0.10 <= x < -0.05)"""
    ...

def PEOE_VSA7(self, x, y=6):
    """
    MOE Charge VSA Descriptor 7 (-0.05 <= x <  0.00)"""
    ...

def PEOE_VSA8(self, x, y=7):
    """
    MOE Charge VSA Descriptor 8 ( 0.00 <= x <  0.05)"""
    ...

def PEOE_VSA9(self, x, y=8):
    """
    MOE Charge VSA Descriptor 9 ( 0.05 <= x <  0.10)"""
    ...

def SMR_VSA1(self, x, y=0):
    """
    MOE MR VSA Descriptor 1 (-inf < x <  1.29)"""
    ...

def SMR_VSA10(self, x, y=9):
    """
    MOE MR VSA Descriptor 10 ( 4.00 <= x < inf)"""
    ...

def SMR_VSA2(self, x, y=1):
    """
    MOE MR VSA Descriptor 2 ( 1.29 <= x <  1.82)"""
    ...

def SMR_VSA3(self, x, y=2):
    """
    MOE MR VSA Descriptor 3 ( 1.82 <= x <  2.24)"""
    ...

def SMR_VSA4(self, x, y=3):
    """
    MOE MR VSA Descriptor 4 ( 2.24 <= x <  2.45)"""
    ...

def SMR_VSA5(self, x, y=4):
    """
    MOE MR VSA Descriptor 5 ( 2.45 <= x <  2.75)"""
    ...

def SMR_VSA6(self, x, y=5):
    """
    MOE MR VSA Descriptor 6 ( 2.75 <= x <  3.05)"""
    ...

def SMR_VSA7(self, x, y=6):
    """
    MOE MR VSA Descriptor 7 ( 3.05 <= x <  3.63)"""
    ...

def SMR_VSA8(self, x, y=7):
    """
    MOE MR VSA Descriptor 8 ( 3.63 <= x <  3.80)"""
    ...

def SMR_VSA9(self, x, y=8):
    """
    MOE MR VSA Descriptor 9 ( 3.80 <= x <  4.00)"""
    ...

def SlogP_VSA1(self, x, y=0):
    """
    MOE logP VSA Descriptor 1 (-inf < x < -0.40)"""
    ...

def SlogP_VSA10(self, x, y=9):
    """
    MOE logP VSA Descriptor 10 ( 0.40 <= x <  0.50)"""
    ...

def SlogP_VSA11(self, x, y=10):
    """
    MOE logP VSA Descriptor 11 ( 0.50 <= x <  0.60)"""
    ...

def SlogP_VSA12(self, x, y=11):
    """
    MOE logP VSA Descriptor 12 ( 0.60 <= x < inf)"""
    ...

def SlogP_VSA2(self, x, y=1):
    """
    MOE logP VSA Descriptor 2 (-0.40 <= x < -0.20)"""
    ...

def SlogP_VSA3(self, x, y=2):
    """
    MOE logP VSA Descriptor 3 (-0.20 <= x <  0.00)"""
    ...

def SlogP_VSA4(self, x, y=3):
    """
    MOE logP VSA Descriptor 4 ( 0.00 <= x <  0.10)"""
    ...

def SlogP_VSA5(self, x, y=4):
    """
    MOE logP VSA Descriptor 5 ( 0.10 <= x <  0.15)"""
    ...

def SlogP_VSA6(self, x, y=5):
    """
    MOE logP VSA Descriptor 6 ( 0.15 <= x <  0.20)"""
    ...

def SlogP_VSA7(self, x, y=6):
    """
    MOE logP VSA Descriptor 7 ( 0.20 <= x <  0.25)"""
    ...

def SlogP_VSA8(self, x, y=7):
    """
    MOE logP VSA Descriptor 8 ( 0.25 <= x <  0.30)"""
    ...

def SlogP_VSA9(self, x, y=8):
    """
    MOE logP VSA Descriptor 9 ( 0.30 <= x <  0.40)"""
    ...

def TPSA(self, *x, **y): ...
def pyLabuteASA(self, mol, includeHs=1):
    """
    calculates Labute’s Approximate Surface Area (ASA from MOE)
    Definition from P. Labute’s article in the Journal of the Chemical Computing Group
    and J. Mol. Graph. Mod.  _18_ 464-477 (2000)"""
    ...

def pyPEOE_VSA_(self, mol, bins=None, force=1):
    """
    Internal Use Only"""
    ...

ptable: Incomplete
bondScaleFacts: Incomplete
mrBins: Incomplete

def pySMR_VSA_(self, mol, bins=None, force=1):
    """
    Internal Use Only"""
    ...

SMR_VSA_: Incomplete
logpBins: Incomplete

def pySlogP_VSA_(self, mol, bins=None, force=1):
    """
    Internal Use Only"""
    ...

SlogP_VSA_: Incomplete
chgBins: Incomplete

def pyPEOE_VSA_(self, mol, bins=None, force=1):
    """
    Internal Use Only"""
    ...

PEOE_VSA_: Incomplete

def pyLabuteASA(self, mol, includeHs=1):
    """
    calculates Labute’s Approximate Surface Area (ASA from MOE)
    Definition from P. Labute’s article in the Journal of the Chemical Computing Group
    and J. Mol. Graph. Mod.  _18_ 464-477 (2000)"""
    ...

def LabuteASA(self, *x, **y): ...
def TPSA(self, *x, **y): ...
