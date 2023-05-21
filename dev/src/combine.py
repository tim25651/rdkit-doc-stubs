import ast
from dataclasses import dataclass
from enum import Enum
from typing import Any

from .names_align import NamesAlign

PRINT_WARNINGS = False


class TypeAST(Enum):
    FN = ast.FunctionDef
    CLS = ast.ClassDef
    ANN = ast.AnnAssign
    ASS = ast.Assign
    IMP = ast.Import
    EXP = ast.Expr

    def __repr__(self) -> str:
        return self.name


EXPnIMP = TypeAST.EXP, TypeAST.IMP


@dataclass(frozen=True)
class KeyAST:
    name: str
    typ: TypeAST
    private: bool
    overload: bool


class SortedAST:
    def __init__(self, ast: ast.Module | ast.ClassDef):
        self.ast = ast
        self.dict = self._get_dict()

    def names_for_type(
        self, _type: tuple[TypeAST] = None, _nottype: tuple[TypeAST] = None
    ):
        names = []
        for key in self.dict:
            if _type and not key.typ in _type:
                continue
            if _nottype and key.typ in _nottype:
                continue
            names.append(key.name)
        return names

    def get_by_name(self, name: str):
        items = [x for key, x in self.dict.items() if key.name == name]
        if len(items) > 2:
            raise ValueError((name, items))
        if len(items) == 2:
            types = {type(x) for x in items}
            if types == {ast.AnnAssign, ast.Assign}:
                return [x for x in items if isinstance(x, ast.Assign)][0]
            elif types == {ast.FunctionDef, ast.AnnAssign}:
                return [x for x in items if isinstance(x, ast.FunctionDef)][0]
            else:
                raise ValueError((types, name, items))
        elif len(items) == 1:
            return items[0]
        else:
            return None

    def __repr__(self) -> str:
        return self.dict.__repr__()

    def _get_dict(self):
        data: dict[KeyAST, Any] = {}
        for value in self.ast.body:
            name, typ = self._get_info(value)
            typ = TypeAST(typ)
            private = name.startswith("_")
            overload = (
                has_overload(value) if isinstance(value, ast.FunctionDef) else False
            )
            key = KeyAST(name, typ, private, overload)
            if key.overload:
                value = [value]
            if key not in data:
                data[key] = value
            else:
                if key.overload:
                    data[key] += value
                ## HARD CODED
                elif key.name in ("binaryHolder"):
                    data[key] = value
                else:
                    raise ValueError(f"Key already in data\n{key=}\n{data=}")
        return data

    @staticmethod
    def _get_info(x):
        if isinstance(x, (ast.FunctionDef, ast.ClassDef)):
            return str(x.name), type(x)
        elif isinstance(x, (ast.AnnAssign)):
            return str(x.target.id), ast.AnnAssign
        elif (
            isinstance(x, (ast.Assign))
            and len(x.targets) == 1
            and isinstance(x.targets[0], ast.Name)
        ):
            return str(x.targets[0].id), ast.Assign
        elif isinstance(x, (ast.Expr)) and isinstance(x.value, ast.Constant):
            return ast.unparse(x.value).strip(), ast.Expr
        elif isinstance(x, (ast.ImportFrom, ast.Import)):
            return ast.unparse(x).strip(), ast.Import

        raise NotImplementedError(x)


def init_combine(new: ast.Module | ast.ClassDef, base: ast.Module | ast.ClassDef):
    assert type(new) == type(base), f"Can't combine module with class: {new=} {base=}"
    module_mode = isinstance(new, ast.Module)
    cls = ast.Module if module_mode else ast.ClassDef
    kwargs = (
        {"type_ignores": new.type_ignores}
        if module_mode
        else {
            "name": new.name,
            "bases": new.bases,
            "keywords": new.keywords,
            "decorator_list": new.decorator_list,
        }
    )
    return cls, kwargs


def has_decorator(fn: ast.FunctionDef, decorator: str):
    return any(isinstance(d, ast.Name) and d.id == decorator for d in fn.decorator_list)


has_overload = lambda fn: has_decorator(fn, "overload")
has_property = lambda fn: has_decorator(fn, "property")


def combine(new: ast.ClassDef | ast.Module, base: ast.ClassDef | ast.Module):
    cls, kwargs = init_combine(new, base)
    body = update_body(new, base)
    return cls(**kwargs, body=body)


def update_function(new: list | ast.FunctionDef, base: list | ast.FunctionDef):
    if isinstance(new, ast.FunctionDef):
        if isinstance(base, list):
            if PRINT_WARNINGS:
                print(
                    f"""
    Warning: multiple definitions of function in base
    {new=}
    {base=}
    """
                )
        new = [new]

    return new


def update_body(new: ast.Module | ast.ClassDef, base: ast.Module | ast.ClassDef):
    new_, base_ = SortedAST(new), SortedAST(base)

    new_elems, base_elems = (
        [x.get_by_name(name) for name in x.names_for_type(_type=EXPnIMP)]
        for x in (new_, base_)
    )
    body = new_elems + base_elems

    new_names, base_names = (x.names_for_type(_nottype=EXPnIMP) for x in (new_, base_))
    combined_names = NamesAlign.align(new_names, base_names)

    other = []
    for name in combined_names:
        new_elem, base_elem = (x.get_by_name(name) for x in (new_, base_))

        if None in (new_elem, base_elem):
            other.append(new_elem or base_elem)

        elif isinstance(new_elem, ast.ClassDef) and isinstance(base_elem, ast.ClassDef):
            other.append(combine(new_elem, base_elem))

        elif isinstance(new_elem, (ast.FunctionDef, list)) and isinstance(
            base_elem, (list, ast.FunctionDef)
        ):
            other.extend(update_function(new_elem, base_elem))
        elif isinstance(new_elem, ast.AnnAssign) and isinstance(
            base_elem, ast.AnnAssign
        ):
            other.append(new_elem)
        elif isinstance(new_elem, ast.FunctionDef) and isinstance(
            base_elem, ast.Assign
        ):
            other.append(base_elem)  ## Could also be the other way around
        elif isinstance(new_elem, ast.FunctionDef) and isinstance(
            base_elem, ast.AnnAssign
        ):
            other.append(new_elem)
        else:
            if isinstance(new_elem, ast.AnnAssign) and isinstance(
                base_elem, ast.FunctionDef
            ):
                if has_property(base_elem):
                    x = new_elem.annotation.id
                    new_elem.annotation.id = "property[" + str(x) + "]"
                    other.append(new_elem)
                    continue

            # raise NotImplementedError((new_elem, base_elem))
            print("NOTIMPL", (new_elem, base_elem))

    return body + other
