from .find import FoundAttribute, FoundClass, FoundFunction, FoundModule, FoundParam

# Indentation: 4 spaces
PAD = "    "


def unparse_param(self: FoundParam):
    type = f": {self.type}" if self.type else ""
    default = f" = {self.default}" if self.default else ""
    return f"{self.name}{type}{default}"


def unparse_function(self: FoundFunction):
    if not self.params or self.params[0].name != "self":
        params = [FoundParam("self", None, None)] + self.params
    else:
        params = self.params

    pad = PAD if self.method else ""
    decorator = (pad + f"@overload\n") if self.overload else ""
    params = ", ".join(map(unparse_param, params))
    ret = f" -> {self.returns}" if self.returns else ""
    doc = f'{pad}{PAD}"""\n{self.doc}\n"""\n' if self.doc and self.doc.strip() else ""
    fn_def = f"{pad}def {self.name}({params}){ret}:\n"
    body = f"{pad}{PAD}...\n"

    return decorator + fn_def + doc + body


def unparse_class(self: FoundClass):
    class_def = f"class {self.name}({self.base}):\n"
    doc = f'{PAD}"""\n{self.doc}\n"""\n' if self.doc.strip() else ""
    # body = f"{PAD}...\n"
    attributes = "\n".join(map(unparse_attribute, self.attributes)) + "\n"
    subclasses = "\n".join(map(unparse_class, self.subclasses)) + "\n"
    methods = "\n".join(map(unparse_function, self.methods)) + "\n"

    body = attributes + subclasses + methods
    body = body if body.strip() else f"{PAD}...\n"
    return shift(class_def + doc + body, self.level)


def shift(string: str, level: int):
    pad = level * PAD
    return pad + f"\n{pad}".join(string.split("\n")) + "\n"


def unparse_attribute(self: FoundAttribute):
    annotation = (
        self.annotation.removeprefix("'").removesuffix("'") if self.annotation else None
    )
    type = (
        self.type if not annotation else f'Annotated[{self.type}, """{annotation}"""]'
    )
    return f"{PAD}{self.name}: {type} = ..."


def unparse_module(self: FoundModule):
    header = f"# Path: https://rdkit.org/docs/source/{self.name}.html\n"
    doc = f'"""\n{self.doc}\n"""\n' if self.doc and self.doc.strip() else ""
    imports = "\n".join(self.imports) + "\n"
    classes = "\n".join(map(unparse_class, self.classes)) + "\n"
    functions = "\n".join(map(unparse_function, self.functions)) + "\n"

    return header + doc + imports + classes + functions
