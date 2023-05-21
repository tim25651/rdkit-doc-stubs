# %%
import ast
import os
import shutil

from tqdm import tqdm

from .src import *


def main():
    out = "out"

    if os.path.exists(out):
        print(f"Output directory {out} already exists!")
        exit(1)
    else:
        print("Output directory: " + out)

    print("Generating stubs for rdkit")
    os.system("stubgen -p rdkit -o " + out)
    shutil.rmtree(".mypy_cache")

    print("Download documentation from RDKit")
    try:
        file = __file__
    except NameError:
        file = os.getcwd()

    TMP_FILE = os.path.dirname(file) + "/data.pkl"

    data = get_data(TMP_FILE)

    for link in tqdm(data):
        doc, functions, classes = Find.module(data[link])
        name = link.removesuffix(".html")
        module = Parse.module(name, doc, functions, classes)

        base_path = find_path(out, name)

        new = ast.parse(format_and_sort(unparse_module(module)))

        if os.path.exists(base_path):
            base = parse_base(base_path)
            combined = combine(new, base)
        else:
            combined = new

        open(base_path, "w").write(unparse_and_format(combined))

    print("Running formatting and sorting")
    os.system("isort " + out)
    os.system("black " + out)

    print("Moving directory to src")
    shutil.move(out + "/rdkit", "src/rdkit-stubs")
    os.rmdir(out)


if __name__ == "__main__":
    main()
