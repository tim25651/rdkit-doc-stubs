"""
rdkit.Chem.Subshape.SubshapeObjects module¶
"""
from _typeshed import Incomplete

class ShapeWithSkeleton(object):
    grid: None = ...
    skelPts: None = ...

    def __init__(self, *args, **kwargs) -> None: ...

class SkeletonPoint(object):
    featmapFeatures: None = ...
    fracVol: float = ...
    location: None = ...
    molFeatures: None = ...
    shapeDirs: None = ...
    shapeMoments: None = ...
    shapeDirs: None = ...
    molFeatures: None = ...
    featmapFeatures: None = ...
    fracVol: float = ...

    def __init__(self, *args, **kwargs) -> None: ...

class ShapeWithSkeleton(object):
    grid: None = ...
    skelPts: None = ...

    def __init__(self, *args, **kwargs) -> None: ...

class SubshapeShape(object):
    shapes: None = ...
    featMap: None = ...
    keyFeat: None = ...
    shapes: None = ...

    def __init__(self, *args, **kwargs) -> None: ...

def DisplaySubshape(self, viewer, shape, name, showSkelPts=True, color=(1, 0, 1)): ...
def DisplaySubshapeSkeleton(
    self, viewer, shape, name, color=(1, 0, 1), colorByOrder=False
): ...
def DisplaySubshape(self, viewer, shape, name, showSkelPts=True, color=(1, 0, 1)): ...