"""
rdkit.Chem.PeriodicTable module¶
A class which stores information from the Periodic Table.
It is not possible to create a PeriodicTable object directly from Python,
use GetPeriodicTable() to get the global table.
The PeriodicTable object can be queried for a variety of properties:

GetAtomicWeight
GetAtomicNumber
GetElementSymbol
GetElementName
GetRvdw (van der Waals radius)
GetRCovalent (covalent radius)
GetDefaultValence
GetValenceList
GetNOuterElecs (number of valence electrons)
GetMostCommonIsotope
GetMostCommonIsotopeMass
GetRb0
GetAbundanceForIsotope
GetMassForIsotope

When it makes sense, these can be queried using either an atomic number (integer)
or an atomic symbol (string)
"""
