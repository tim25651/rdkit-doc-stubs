"""
rdkit.ML.Descriptors.Parser module¶
The “parser” for compound descriptors.
I almost hesitate to document this, because it’s not the prettiest
thing the world has ever seen… but it does work (for at least some
definitions of the word).
Rather than getting into the whole mess of writing a parser for the
compound descriptor expressions, I’m just using string substitutions
and python’s wonderful ability to eval code.
It would probably be a good idea at some point to replace this with a
real parser, if only for the flexibility and intelligent error
messages that would become possible.
The general idea is that we’re going to deal with expressions where
atomic descriptors have some kind of method applied to them which
reduces them to a single number for the entire composition.  Compound
descriptors (those applicable to the compound as a whole) are not
operated on by anything in particular (except for standard math stuff).
Here’s the general flow of things:

Composition descriptor references ($a, $b, etc.) are replaced with the
corresponding descriptor names using string substitution.
(_SubForCompoundDescriptors)
Atomic descriptor references ($1, $2, etc) are replaced with lookups
into the atomic dict with “DEADBEEF” in place of the atom name.
(_SubForAtomicVars)
Calls to Calculator Functions are augmented with a reference to
the composition and atomic dictionary
(_SubMethodArgs)

NOTE:

anytime we don’t know the answer for a descriptor, rather than
throwing a (completely incomprehensible) exception, we just return
-666.  So bad descriptor values should stand out like sore thumbs.
"""
from math import *

from _typeshed import Incomplete
from rdkit import RDConfig as RDConfig

AVG = MEAN

def CalcMultipleCompoundsDescriptor(self, composVect, argVect, atomDict, propDictList):
    """
    calculates the value of the descriptor for a list of compounds
    ARGUMENTS:

    composVect: a vector of vector/tuple containing the compositioninformation.
    See _CalcSingleCompoundDescriptor()_ for an explanation of the elements.

    argVect: a vector/tuple with three elements:

    AtomicDescriptorNames:  a list/tuple of the names of the

    atomic descriptors being used. These determine the
    meaning of $1, $2, etc. in the expression

    CompoundDsscriptorNames:  a list/tuple of the names of the

    compound descriptors being used. These determine the
    meaning of $a, $b, etc. in the expression

    Expr: a string containing the expression to be used to

    evaluate the final result.

    atomDict:a dictionary of atomic descriptors.  Each atomic entry is
    another dictionary containing the individual descriptors
    and their values

    propVectList:a vector of vectors of descriptors for the composition.

    RETURNS:

    a vector containing the values of the descriptor for each
    compound.  Any given entry will be -666 if problems were
    encountered"""
    ...

def CalcSingleCompoundDescriptor(self, compos, argVect, atomDict, propDict):
    """
    calculates the value of the descriptor for a single compound
    ARGUMENTS:

    compos: a vector/tuple containing the compositioninformation… in the form:
    ‘[(“Fe”,1.),(“Pt”,2.),(“Rh”,0.02)]’

    argVect: a vector/tuple with three elements:

    AtomicDescriptorNames:  a list/tuple of the names of the

    atomic descriptors being used. These determine the
    meaning of $1, $2, etc. in the expression

    CompoundDescriptorNames:  a list/tuple of the names of the

    compound descriptors being used. These determine the
    meaning of $a, $b, etc. in the expression

    Expr: a string containing the expression to be used to

    evaluate the final result.

    atomDict:a dictionary of atomic descriptors.  Each atomic entry is
    another dictionary containing the individual descriptors
    and their values

    propVect:a list of descriptors for the composition.

    RETURNS:

    the value of the descriptor, -666 if a problem was encountered

    NOTE:

    because it takes rather a lot of work to get everything setup to calculate a descriptor, if you are calculating the
    same descriptor for multiple compounds, you probably want to
    be calling _CalcMultipleCompoundsDescriptor()_."""
    ...

def DEV(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the average deviation of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

knownMethods: Incomplete

def HAS(self, strArg, composList, atomDict):
    """
    Calculator Method
    does a string search
    Arguments

    strArg: the arguments in string form
    composList: the composition vector
    atomDict: the atomic dictionary

    Returns

    1 or 0"""
    ...

def MAX(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the maximum value of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

def SUM(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the sum of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

def MEAN(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the average of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

AVG = MEAN

def DEV(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the average deviation of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

def MIN(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the minimum value of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

def SUM(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the sum of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

def MAX(self, strArg, composList, atomDict):
    """
    Calculator Method
    calculates the maximum value of a descriptor across a composition
    Arguments

    strArg: the arguments in string form
    compos: the composition vector
    atomDict: the atomic dictionary

    Returns

    a float"""
    ...

def CalcSingleCompoundDescriptor(self, compos, argVect, atomDict, propDict):
    """
    calculates the value of the descriptor for a single compound
    ARGUMENTS:

    compos: a vector/tuple containing the compositioninformation… in the form:
    ‘[(“Fe”,1.),(“Pt”,2.),(“Rh”,0.02)]’

    argVect: a vector/tuple with three elements:

    AtomicDescriptorNames:  a list/tuple of the names of the

    atomic descriptors being used. These determine the
    meaning of $1, $2, etc. in the expression

    CompoundDescriptorNames:  a list/tuple of the names of the

    compound descriptors being used. These determine the
    meaning of $a, $b, etc. in the expression

    Expr: a string containing the expression to be used to

    evaluate the final result.

    atomDict:a dictionary of atomic descriptors.  Each atomic entry is
    another dictionary containing the individual descriptors
    and their values

    propVect:a list of descriptors for the composition.

    RETURNS:

    the value of the descriptor, -666 if a problem was encountered

    NOTE:

    because it takes rather a lot of work to get everything setup to calculate a descriptor, if you are calculating the
    same descriptor for multiple compounds, you probably want to
    be calling _CalcMultipleCompoundsDescriptor()_."""
    ...

def CalcMultipleCompoundsDescriptor(self, composVect, argVect, atomDict, propDictList):
    """
    calculates the value of the descriptor for a list of compounds
    ARGUMENTS:

    composVect: a vector of vector/tuple containing the compositioninformation.
    See _CalcSingleCompoundDescriptor()_ for an explanation of the elements.

    argVect: a vector/tuple with three elements:

    AtomicDescriptorNames:  a list/tuple of the names of the

    atomic descriptors being used. These determine the
    meaning of $1, $2, etc. in the expression

    CompoundDsscriptorNames:  a list/tuple of the names of the

    compound descriptors being used. These determine the
    meaning of $a, $b, etc. in the expression

    Expr: a string containing the expression to be used to

    evaluate the final result.

    atomDict:a dictionary of atomic descriptors.  Each atomic entry is
    another dictionary containing the individual descriptors
    and their values

    propVectList:a vector of vectors of descriptors for the composition.

    RETURNS:

    a vector containing the values of the descriptor for each
    compound.  Any given entry will be -666 if problems were
    encountered"""
    ...
