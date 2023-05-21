import ast
import os
import re

import black
import isort

NONE_PATTERN = re.compile(r"(?<!>\s)None:")


def format_and_sort(code: str):
    code = isort.code(code)
    code = black.format_str(code, mode=black.Mode(is_pyi=True))
    return code


def parse_base(path: str) -> ast.Module:
    base = open(path).read()

    ## HARD CODED
    base = base.replace("(anonymousnamespace)", "anonymousnamespace")
    base = re.sub(NONE_PATTERN, "None_:", base)
    base = base.replace("used to filter multiple conformations", "Any")

    return ast.parse(base)


def unparse_and_format(x: ast.AST) -> str:
    code = ast.unparse(x).strip()
    return format_and_sort(code)


def find_path(dir: str, name: str) -> str:
    path = dir + "/" + name.replace(".", "/")

    if os.path.isdir(path):
        path += "/__init__.pyi"
    else:
        path += ".pyi"

    return path
