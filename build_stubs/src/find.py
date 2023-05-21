import re
from dataclasses import dataclass
from typing import Callable

from bs4 import BeautifulSoup, Tag

DOC_PATTERN: Callable[[str], str] = (
    lambda x: r"(" + x + r"\(.*\)(?:\s*->\s*(.*?))?\s*:(?!:))"
)


@dataclass
class FoundFunction:
    name: str
    params: list[str]
    returns: str
    doc: str
    overloads: list[str] | None = None
    method: bool | None = None
    overload: bool | None = None


@dataclass
class FoundClass:
    name: str
    base: str
    doc: str
    subclasses: list["FoundClass"]
    methods: list[FoundFunction]
    attributes: list[str]
    level: int | None = None


@dataclass
class FoundModule:
    name: str
    doc: str
    imports: set[str]
    functions: list[FoundFunction]
    classes: list[FoundClass]


@dataclass
class FoundParam:
    name: str
    type: str
    default: str


@dataclass
class FoundAttribute:
    name: str
    type: str
    annotation: str


class Find:
    """Collection of functions to find elements in the HTML tree of the RDKit documentation"""

    @staticmethod
    def name(node: Tag) -> str:
        return node.find("span", class_="sig-name descname").get_text().strip()

    @staticmethod
    def base(node: Tag) -> str:
        base = node.find("code", class_="py-class")
        if base:
            return base.get_text().strip()
        else:
            return ""

    @staticmethod
    def params(node: Tag) -> list[str]:
        params = node.find_all("em", class_="sig-param")
        return [p.get_text().strip() for p in params]

    @staticmethod
    def attributes(node: Tag) -> list[str]:
        attributes = node.findChildren("dl", class_="py attribute", recursive=False)
        return [a.get_text().strip().replace("¶", "") for a in attributes]

    @staticmethod
    def doc(node: Tag, classes: list[str]) -> str:
        first = node.find("dl", {"class": classes})

        if first:
            doc_nodes = first.findPreviousSiblings()[::-1]
        else:
            doc_nodes = node.findChildren(recursive=False)

        doc = "\n".join([x.get_text() for x in doc_nodes[1:]])

        return re.sub(r"\n\n\n+", "\n\n", doc).strip()

    @staticmethod
    def fn_doc(node: Tag) -> str:
        text = node.find("dd").get_text()
        return re.sub(r"\n\n\n+", "\n\n", text).strip()

    @staticmethod
    def return_(node: Tag) -> str:
        ret = node.find("span", class_="sig-return")
        if ret:
            return ret.get_text().strip()
        else:
            return ""

    @staticmethod
    def overloads(doc: str, name: str) -> list[str]:
        overloads: list[str] = []
        while doc_contains_func := re.search(DOC_PATTERN(name), doc):
            fn_def = doc_contains_func.group(1)
            # :: = C++ synthax
            if not "::" in fn_def:
                ## HARD CODED (Alternative: exlude lines starting with >>>)
                if not "EnumerateHeterocycles" in fn_def:
                    overloads.append(fn_def.strip())
            _, doc = doc.split(fn_def, 1)
        return overloads

    @staticmethod
    def functions(node: Tag, method=False) -> list[FoundFunction]:
        class_ = "py method" if method else "py function"
        return [
            Find.function(node)
            for node in node.findChildren("dl", class_=class_, recursive=False)
        ]

    @staticmethod
    def function(node: Tag) -> FoundFunction:
        name = Find.name(node)
        params = Find.params(node)
        ret = Find.return_(node)
        doc = Find.fn_doc(node)
        overloads = Find.overloads(doc, name)
        return FoundFunction(name, params, ret, doc, overloads)

    @staticmethod
    def classes(node: Tag) -> list[FoundClass]:
        classes = node.findChildren("dl", class_="py class", recursive=False)
        exceptions = node.findChildren("dl", class_="py exception", recursive=False)
        classes = classes + exceptions
        return [Find.class_(c) for c in classes]

    @staticmethod
    def class_(node: Tag) -> FoundClass:
        name = Find.name(node)

        node = node.find("dd")

        base = Find.base(node)
        doc = Find.doc(node, ["py method", "py attribute", "py class"])

        methods = Find.functions(node, method=True)
        attributes = Find.attributes(node)
        subclasses = Find.classes(node)

        return FoundClass(name, base, doc, subclasses, methods, attributes)

    def module(
        data: str,
    ) -> tuple[str, list[FoundClass], list[FoundFunction],]:
        soup = BeautifulSoup(data, "html.parser")
        node = soup.find("section")

        doc = Find.doc(node, ["py exception", "py class", "py function"])

        classes = Find.classes(node)

        functions = Find.functions(node)

        return doc, functions, classes
