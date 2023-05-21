# %%
import re


class InferTypes:
    def __init__(self, module_name: str) -> None:
        self.imports: set[str] = set()
        self.module_name = module_name
        self.repl_counter = 0
        self.repl: dict[int, str] = {}

    def store_match(self, match: re.Match) -> str:
        self.repl_counter += 1
        self.repl[self.repl_counter] = match.group(0)
        return f"__§§{self.repl_counter}__"

    def replace(self, x: str):
        prev = None
        while prev != x:
            prev = x
            x = re.sub(r"('.*?')", self.store_match, x, 1)
            for p, s in zip("[(}", "])}"):
                x = re.sub(f"\\{p}.*?\\{s}", self.store_match, x, 1)

        return x

    def first_element(self, x: str):
        x = self.replace(x)
        return x.split(",")[0].strip()

    def first_element_dict(self, x: str):
        x = self.first_element(x)
        k, v = x.split(":")
        return k.strip(), v.strip()

    def infer_lts(self, x: str, ps="()"):
        p, s = ps
        x = x.removeprefix(p).removesuffix(s)
        if x == "":
            self.imports.add("Any")
            return "Any"
        x = self.first_element(x)
        return self.infer_type(x)

    def infer_dict(self, x: str):
        x = x.removeprefix("{").removesuffix("}")
        if x == "":
            self.imports.add("Any")
            return "Any, Any"
        key, value = self.first_element_dict(x)

        return f"{self.infer_type(key)}, {self.infer_type(value)}"

    def infer_classvar(self, x: str):
        x = x.removeprefix("<").removesuffix(">")
        if "object" in x:
            value = x.replace("object", "").strip()
            value = re.sub(r"(at 0x.*)", "", value).strip()
            if value.startswith(self.module_name + "."):
                value = value.removeprefix(self.module_name + ".")
                return value
            else:
                self.imports.add(value)
                return ".".join(value.split(".")[-2:])
        if "module" in x:
            self.imports.add("ModuleType")
            return "ModuleType"
        if x.startswith("HashLayer"):
            return "HashLayer"

        raise ValueError("CCC" + x)

    def infer_stacked_type(self, value: str):
        if value.startswith("[") and value.endswith("]"):
            return f"list[{self.infer_lts(value, ps='[]')}]"

        if value.startswith("(") and value.endswith(")"):
            return f"tuple[{self.infer_lts(value, ps='()')}]"

        if value.startswith("{") and value.endswith("}"):
            return f"dict[{self.infer_dict(value)}]"

        if value.startswith("<") and value.endswith(">"):
            return self.infer_classvar(value)

    def infer_type(self, value: str):
        if value is None or value in ("", "Any"):
            self.imports.add("Any")
            return "Any"

        if "§§" in value:
            for k, v in self.repl.items():
                value = value.replace(f"__§§{k}__", v)
            return self.infer_type(value)

        if value == "None":
            return "None"

        if value in ("True", "False"):
            return "bool"

        if value.startswith("'") and value.endswith("'"):
            return "str"

        if value.startswith("Alias for"):
            self.imports.add("_Alias")
            return "_Alias"
        if value.startswith("rdkit."):
            self.imports.add(value)
            return value.split(".")[-2]

        try:
            int(value)
            return "int"
        except:
            pass
        try:
            float(value)
            return "float"
        except:
            pass

        match = re.match(r"^([a-z, A-Z]*)\(.*\)$", value)
        if match and (match_value := match.group(1)):
            self.imports.add(match_value)
            return match_value

        stacked_type = self.infer_stacked_type(value)

        if stacked_type:
            return stacked_type

        else:
            raise ValueError(f"Cannot infer type from value: {value}")


# %%
