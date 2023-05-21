"""
rdkit.Chem.Descriptors3D module¶
Descriptors derived from a molecule’s 3D structure
"""
from _typeshed import Incomplete
from rdkit.Chem import rdMolDescriptors as rdMolDescriptors

def Asphericity(self, *x, **y):
    """
    molecular asphericity

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:0.5 * ((pm3-pm2)**2 + (pm3-pm1)**2 + (pm2-pm1)**2)/(pm1**2+pm2**2+pm3**2)

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def Eccentricity(self, *x, **y):
    """
    molecular eccentricity

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:sqrt(pm3**2 -pm1**2) / pm3

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def InertialShapeFactor(self, *x, **y):
    """
    Inertial shape factor

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:pm2 / (pm1*pm3)

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def NPR1(self, *x, **y):
    """
    Normalized principal moments ratio 1 (=I1/I3)

    from Sauer and Schwarz JCIM 43:987-1003 (2003)
    https://doi.org/10.1021/ci025599w

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def NPR2(self, *x, **y):
    """
    Normalized principal moments ratio 2 (=I2/I3)

    from Sauer and Schwarz JCIM 43:987-1003 (2003)
    https://doi.org/10.1021/ci025599w

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def PMI1(self, *x, **y):
    """
    First (smallest) principal moment of inertia
    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def PMI2(self, *x, **y):
    """
    Second principal moment of inertia
    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def PMI3(self, *x, **y):
    """
    Third (largest) principal moment of inertia
    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def NPR1(self, *x, **y):
    """
    Normalized principal moments ratio 1 (=I1/I3)

    from Sauer and Schwarz JCIM 43:987-1003 (2003)
    https://doi.org/10.1021/ci025599w

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def NPR2(self, *x, **y):
    """
    Normalized principal moments ratio 2 (=I2/I3)

    from Sauer and Schwarz JCIM 43:987-1003 (2003)
    https://doi.org/10.1021/ci025599w

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def RadiusOfGyration(self, *x, **y):
    """
    Radius of gyration

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:for planar molecules: sqrt( sqrt(pm3*pm2)/MW )
    for nonplanar molecules: sqrt( 2*pi*pow(pm3*pm2*pm1,1/3)/MW )

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def InertialShapeFactor(self, *x, **y):
    """
    Inertial shape factor

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:pm2 / (pm1*pm3)

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def Eccentricity(self, *x, **y):
    """
    molecular eccentricity

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:sqrt(pm3**2 -pm1**2) / pm3

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def Asphericity(self, *x, **y):
    """
    molecular asphericity

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:0.5 * ((pm3-pm2)**2 + (pm3-pm1)**2 + (pm2-pm1)**2)/(pm1**2+pm2**2+pm3**2)

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use
    useAtomicMasses: (optional) toggles use of atomic masses in the
    calculation. Defaults to True"""
    ...

def SpherocityIndex(self, *x, **y):
    """
    Molecular spherocityIndex

    from Todeschini and Consoni “Descriptors from Molecular Geometry”
    Handbook of Chemoinformatics
    https://doi.org/10.1002/9783527618279.ch37

    Definition:3 * pm1 / (pm1+pm2+pm3) where the moments are calculated without weights

    Arguments

    inMol: a molecule
    confId: (optional) the conformation ID to use"""
    ...
