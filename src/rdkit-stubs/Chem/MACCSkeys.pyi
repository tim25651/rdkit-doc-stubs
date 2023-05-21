"""
rdkit.Chem.MACCSkeys module¶
SMARTS definitions for the publicly available MACCS keys
and a MACCS fingerprinter
I compared the MACCS fingerprints generated here with those from two
other packages (not MDL, unfortunately). Of course there are
disagreements between the various fingerprints still, but I think
these definitions work pretty well. Some notes:

most of the differences have to do with aromaticity

2) there’s a discrepancy sometimes because the current RDKit
definitions do not require multiple matches to be distinct. e.g. the
SMILES C(=O)CC(=O) can match the (hypothetical) key O=CC twice in my
definition. It’s not clear to me what the correct behavior is.
3) Some keys are not fully defined in the MDL documentation
4) Two keys, 125 and 166, have to be done outside of SMARTS.
5) Key 1 (ISOTOPE) isn’t defined
Rev history:
2006 (gl): Original open-source release
May 2011 (gl): Update some definitions based on feedback from Andrew Dalke
"""
from _typeshed import Incomplete
from rdkit import Chem as Chem
from rdkit import DataStructs as DataStructs
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors

smartsPatts: Incomplete
maccsKeys: Incomplete
GenMACCSKeys: Incomplete
FingerprintMol: Incomplete
