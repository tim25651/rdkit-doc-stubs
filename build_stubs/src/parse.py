import re

from .find import FoundAttribute, FoundClass, FoundFunction, FoundModule, FoundParam
from .infer_types import InferTypes
from .modules import BUILTIN_TYPES, MODULES

PARAM_PATTERN = re.compile(r"^(\((.*?)\))?(.*?)(=(.*))?$")
ATTRIBUTE_PATTERN = re.compile(
    r"^(.*?)(?:\s?(?:=|\r\n|\r|\n)\s?(.*)((?:(?:(?:\r\n|\r|\n)?)(?:.*))*))?$"
)

FN_DEF_PATTERN = re.compile(r"(?:.*?\()(.*)(?:\))(?:(?:\s*->\s*)(.*)(?:\s*:))?")


class Parse:
    """Collection of functions to parse found elements to their python equivalent"""

    def __init__(self, module: str) -> None:
        self.inferer = InferTypes(module)
        self.name = module

    @staticmethod
    def remove_keywords(name: str) -> str:
        if name in ("from", "in", "class", "return", "def", "as", "import"):
            name = name + "_"
        return name

    def remove_module(self, x: str) -> str | None:
        if x.startswith(self.name + "."):
            x = x.removeprefix(self.name + ".")
            return x

    def remove_rdkit(self, x: str) -> str | None:
        if x.startswith("rdkit."):
            self.inferer.imports.add(x)
            return ".".join(x.split(".")[-2:])

    def base(self, value: str) -> str:
        if not value:
            return ""

        if class_in_current_module := self.remove_module(value):
            value = class_in_current_module

        elif class_in_rdkit := self.remove_rdkit(value):
            value = class_in_rdkit
        else:
            for module in ("Boost.Python", "sqlalchemy", "enum"):
                if value.startswith(module + "."):
                    self.inferer.imports.add(module)
                    return value

        return value

    def return_(self, ret: str) -> str:
        if not ret:
            return ""

        x = ret.removeprefix("→").removesuffix(":").strip()

        ## HARD CODED
        # Errors in the documentation

        if "in the interval" in x:
            x = "float"
        if "returns True" in x:
            x = "bool"

        # Add to the import list
        self.inferer.imports.add(x)

        return x

    def classes(self, classes: list[FoundClass]) -> list[FoundClass]:
        return [self.class_(c) for c in classes]

    def functions(self, functions: list[FoundFunction]) -> list[FoundClass]:
        return [self.function(f) for f in functions]

    def class_(self, class_: FoundClass, level=0) -> FoundClass:
        class_.base = self.base(class_.base)
        class_.subclasses = [self.class_(c, level=level + 1) for c in class_.subclasses]
        class_.methods = [self.function(m, method=True) for m in class_.methods]
        class_.methods = [m for fn in class_.methods for m in fn]
        class_.attributes = [self.attribute(a) for a in class_.attributes]
        class_.level = level
        return class_

    def function(
        self,
        fn: FoundFunction,
        method=False,
    ) -> list[FoundFunction]:
        # params = [parse_param(p, module_name) for p in params]

        if not fn.overloads:
            fn.params = [self.param(p) for p in fn.params]
            fn.returns = self.return_(fn.returns)
            fn.overloads = []
            fn.method = method
            fn.overload = False
            return [fn]

        else:
            self.inferer.imports.add("overload")

            fns = []
            doc = fn.doc
            for ix, overload in enumerate([None] + fn.overloads):
                d, doc = (
                    (doc, None)
                    if ix == len(fn.overloads)
                    else doc.split(fn.overloads[ix], 1)
                )

                params, ret = (
                    self.overload(overload) if ix != 0 else (fn.params, fn.returns)
                )

                params = [self.param(p) for p in params]
                ret = self.return_(ret)

                fns.append(
                    FoundFunction(
                        fn.name, params, ret, doc, d, method=method, overload=True
                    )
                )

            return fns

    def overload(self, overload: str) -> tuple[list[str], str]:
        # Split function def in overload in params and return value
        matched_params, matched_ret = FN_DEF_PATTERN.match(overload).groups()
        # Optional params are seperated with brackets
        # Split them from the required params
        req_params, *opt_params = matched_params.split("[,")
        # Split the required params by comma
        req_params = req_params.split(",")
        combined_params = req_params + opt_params
        # Clean up the params, remove whitespaces and brackets and wrong apostrophes
        params = [p.strip().replace("’", "'") for p in combined_params]
        # Remove the brackets from the last optional params
        if "[" in params[-1]:
            print("WARNING: last optional param is a list, potential bug")
        while params[-1].endswith("]"):
            params[-1] = params[-1].removesuffix("]").strip()

        ret = matched_ret.strip()

        return params, ret

    def attribute(self, text: str) -> FoundAttribute:
        name, value, annotation = ATTRIBUTE_PATTERN.match(text).groups()
        value = value or ""
        type = self.inferer.infer_type(value.strip())
        name = name.strip()

        ## HARD CODED
        if name == "None":
            name = "None_"

        if annotation:
            annotation = annotation.strip()
            self.inferer.imports.add("Annotated")

        return FoundAttribute(name, type, annotation)

    def param(self, text: str) -> FoundParam:
        _, typ, name, _, default = PARAM_PATTERN.match(text).groups()

        # Incositent documentation: type annotated python-style
        if ":" in name:
            name, typ = name.split(":")

        name = name.strip()
        name = self.remove_keywords(name)  # Reserved keywords

        if typ:
            typ = typ.strip()

            if class_in_current_module := self.remove_module(typ):
                typ = class_in_current_module

            ## HARD CODED
            elif typ == "FilterCatalogs":
                typ = "FilterCatalogParams.FilterCatalogs"
            else:
                # Class in other module
                self.inferer.imports.add(typ)

        if default:
            default = default.strip()

            if class_in_current_module := self.remove_module(default):
                default = class_in_current_module

            elif class_in_other_module := self.remove_rdkit(default):
                default = class_in_other_module

            else:
                if default.startswith("<") and default.endswith(">"):
                    default = default.removeprefix("<").removesuffix(">").strip()

                    if "Boost.Python.function" in default:
                        default, typ = "...", "Callable"
                        self.inferer.imports.add("Callable")

                    elif "function" in default:
                        default = default.replace("function", "").strip()

                        if default == "EuclideanDist":
                            default = "DistFunctions.EuclideanDist"
                            self.inferer.imports.add(
                                "rdkit.ML.KNN.DistFunctions.EuclideanDist"
                            )
                        elif default == "ID3Boot":
                            default = "ID3.ID3Boot"
                            self.inferer.imports.add("rdkit.ML.DecTree.ID3.ID3Boot")

                    elif "class" in default:
                        default = default.replace("class", "").replace("'", "").strip()

                    elif "TextIOWrapper" in default:
                        default, typ = "...", "_io.TextIOWrapper"
                        self.inferer.imports.add("_io")
                    else:
                        default, typ = "...", self.inferer.infer_classvar(default)

        if "-->" in name:
            name, typ = "x", "list[int]"

        if " " and "MHFP" in name:
            if "integers" in name:
                name, typ = "mhfp", "list[int]"
            elif "strings" in name:
                name, typ = "mhfp", "list[str]"

        return FoundParam(name, typ, default)

    @staticmethod
    def module(
        name: str, doc: str, functions: list[FoundFunction], classes: list[FoundClass]
    ):
        parser = Parse(name)

        classes = parser.classes(classes)
        functions = parser.functions(functions)
        functions = [f for fn in functions for f in fn]
        imports = parser.imports([c.name for c in classes])

        return FoundModule(name, doc, imports, functions, classes)

    def imports(self, classes: list[str]):
        missing_imports = set()
        for import_ in self.inferer.imports:
            if not import_:
                continue

            cls = None

            if import_.startswith("rdkit."):
                cls = import_.split(".")[-2]

            if (
                import_ in BUILTIN_TYPES
                or import_ in classes
                or (cls and cls in classes)
            ):
                continue

            if import_.startswith("Boost.Python."):
                missing_imports.add("Boost.Python")
                continue

            missing_imports.add(import_)

        out = set()
        for module in missing_imports:
            if module in MODULES:
                parent = MODULES[module]
            else:
                try:
                    parent = ".".join(module.split(".")[:-2])
                    module = module.split(".")[-2]
                except IndexError as e:
                    print(module)
                    raise e

            if module == self.name:
                continue

            if parent:
                out.add(f"from {parent} import {module}")
            else:
                out.add(f"import {module}")

        return out
