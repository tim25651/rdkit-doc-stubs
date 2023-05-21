"""
rdkit.ML.Composite.AdjustComposite module¶
functionality to allow adjusting composite model contents
"""
from _typeshed import Incomplete

def BalanceComposite(
    self, model, set1, set2, weight, targetSize, names1=None, names2=None
):
    """
    adjusts the contents of the composite model so as to maximize
    the weighted classification accuracty across the two data sets.
    The resulting composite model, with _targetSize_ models, is returned.
    Notes:

    if _names1_ and _names2_ are not provided, _set1_ and _set2_ should
    have the same ordering of columns and _model_ should have already
    have had _SetInputOrder()_ called."""
    ...
