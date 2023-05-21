"""
rdkit.ML.Neural.ActFuncs module¶
Activation functions for neural network nodes
Activation functions should implement the following API:

_Eval(x)_: returns the value of the function at a given point
_Deriv(x)_: returns the derivative of the function at a given point

The current Backprop implementation also requires:

_DerivFromVal(val)_: returns the derivative of the function when itsvalue is val

In all cases _x_ is a float as is the value returned.
"""
from _typeshed import Incomplete

class ActFunc(object):
    """
    “virtual base class” for activation functions"""

    ...

    def __call__(self, x): ...

class Sigmoid(ActFunc):
    """
    the standard sigmoidal function"""

    def Eval(self, x): ...
    def Deriv(self, x): ...
    def DerivFromVal(self, val): ...
    def Eval(self, x): ...
    beta: Incomplete

    def __init__(self, beta: float = ...) -> None: ...

class TanH(ActFunc):
    """
    the standard hyperbolic tangent function"""

    def Eval(self, x): ...
    def Deriv(self, x): ...
    def DerivFromVal(self, val): ...
    def Eval(self, x): ...
    beta: Incomplete

    def __init__(self, beta: float = ...) -> None: ...
