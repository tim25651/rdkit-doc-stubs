import ast


def __repr__(self):
    return self.__class__.__name__ + ": " + self.name


def __repr2__(self):
    return self.__class__.__name__ + ": " + f"{self.value}"


def __repr3__(self):
    return self.__class__.__name__ + ": " + self.module


def __repr4__(self):
    return self.__class__.__name__ + ": " + f"{self.names}"


def __repr5__(self):
    return self.__class__.__name__ + ": " + self.id


def __repr6__(self):
    return self.__class__.__name__ + ": " + f"{self.target}"


def __repr7__(self):
    return self.__class__.__name__ + ": " + f"{self.targets}"


ast.FunctionDef.__repr__ = __repr__
ast.ClassDef.__repr__ = __repr__
ast.Expr.__repr__ = __repr2__
ast.Constant.__repr__ = __repr2__
ast.ImportFrom.__repr__ = __repr3__
ast.Import.__repr__ = __repr4__
ast.alias.__repr__ = __repr__
ast.Name.__repr__ = __repr5__
ast.AnnAssign.__repr__ = __repr6__
ast.Assign.__repr__ = __repr7__
