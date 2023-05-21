"""
rdkit.ML.Neural.Trainers module¶
Training algorithms for feed-forward neural nets
Unless noted otherwise, algorithms and notation are taken from:
“Artificial Neural Networks: Theory and Applications”,

Dan W. Patterson, Prentice Hall, 1996
"""
from _typeshed import Incomplete

class BackProp(Trainer):
    """
    implement back propagation (algorithm on pp 153-154 of Patterson)

    I don’t think that I’ve made any assumptions about the connectivity ofthe net (i.e. full connectivity between layers is not required).

    NOTE: this code is currently making the assumption that the activationfunctions on the nodes in the network are capable of calculating their
    derivatives using only their values (i.e. a DerivFromVal method should
    exist).  This shouldn’t be too hard to change.

    Constructor
    Arguments

    speed: the speed parameter for back prop training
    momentum: the momentum term for back prop training
    Not currently used"""

    oldDeltaW: Incomplete

    def StepUpdate(self, example, net, resVect=None):
        """
        does a BackProp step based upon the example
        Arguments

        example: a 2-tuple:
        a list of variable values values
        a list of result values (targets)

        net: a _Network_ (or something supporting the same API)
        resVect: if this is nonzero, then the network is not required to
        classify the _example_

        Returns

        the backprop error from _network_ before the update

        Note

        In case it wasn’t blindingly obvious, the weights in _network_ are modified
        in the course of taking a backprop step."""
        ...
    def TrainOnLine(
        self, examples, net, maxIts=5000, errTol=0.1, useAvgErr=1, silent=0
    ):
        """
        carries out online training of a neural net

        The definition of online training is that the network is updated aftereach example is presented.

        Arguments

        examples: a list of 2-tuple:
        a list of variable values values
        a list of result values (targets)

        net: a _Network_ (or something supporting the same API)
        maxIts: the maximum number of training epochs (see below for definition) to be
        run
        errTol: the tolerance for convergence
        useAvgErr: if this toggle is nonzero, then the error at each step will be
        divided by the number of training examples for the purposes of checking
        convergence.
        silent: controls the amount of visual noise produced as this runs.

        Note

        a training epoch is one complete pass through all the training examples"""
        ...
    speed: Incomplete
    momentum: Incomplete

    def __init__(self, speed: float = ..., momentum: float = ...) -> None: ...

class Trainer(object):
    """
    “virtual base class” for network trainers"""

    ...
    ...

class BackProp(Trainer):
    """
    implement back propagation (algorithm on pp 153-154 of Patterson)

    I don’t think that I’ve made any assumptions about the connectivity ofthe net (i.e. full connectivity between layers is not required).

    NOTE: this code is currently making the assumption that the activationfunctions on the nodes in the network are capable of calculating their
    derivatives using only their values (i.e. a DerivFromVal method should
    exist).  This shouldn’t be too hard to change.

    Constructor
    Arguments

    speed: the speed parameter for back prop training
    momentum: the momentum term for back prop training
    Not currently used"""

    oldDeltaW: Incomplete

    def StepUpdate(self, example, net, resVect=None):
        """
        does a BackProp step based upon the example
        Arguments

        example: a 2-tuple:
        a list of variable values values
        a list of result values (targets)

        net: a _Network_ (or something supporting the same API)
        resVect: if this is nonzero, then the network is not required to
        classify the _example_

        Returns

        the backprop error from _network_ before the update

        Note

        In case it wasn’t blindingly obvious, the weights in _network_ are modified
        in the course of taking a backprop step."""
        ...
    def TrainOnLine(
        self, examples, net, maxIts=5000, errTol=0.1, useAvgErr=1, silent=0
    ):
        """
        carries out online training of a neural net

        The definition of online training is that the network is updated aftereach example is presented.

        Arguments

        examples: a list of 2-tuple:
        a list of variable values values
        a list of result values (targets)

        net: a _Network_ (or something supporting the same API)
        maxIts: the maximum number of training epochs (see below for definition) to be
        run
        errTol: the tolerance for convergence
        useAvgErr: if this toggle is nonzero, then the error at each step will be
        divided by the number of training examples for the purposes of checking
        convergence.
        silent: controls the amount of visual noise produced as this runs.

        Note

        a training epoch is one complete pass through all the training examples"""
        ...
    speed: Incomplete
    momentum: Incomplete

    def __init__(self, speed: float = ..., momentum: float = ...) -> None: ...
