# %%
import ast
import logging
import os
import shutil

from tqdm import tqdm

from .src import *

# If True, the script will exit if mypy's stubs already exist
# Keeping stubs can cause errors if files has been changed already
ERROR_ON_EXISTING_STUBS = False


def main():
    logging.basicConfig(level=logging.INFO)

    out = "tmp_stubs"

    if os.path.exists(out) and ERROR_ON_EXISTING_STUBS:
        logging.error(f"Temporary directory {out} already exists!")
        exit(1)
    else:
        logging.info("Temporary directory: " + out)

    # Generate stubs for rdkit
    if not os.path.exists(out):
        logging.info("Generating stubs for rdkit")
        os.system("stubgen -p rdkit -o " + out)
        shutil.rmtree(".mypy_cache")
    else:
        logging.info("Stubs for rdkit already exist")

    # Download documentation from RDKit
    try:
        file = __file__
    except NameError:
        file = os.getcwd()

    TMP_FILE = os.path.dirname(file) + "/data.pkl"
    logging.info("Using temporary file: " + TMP_FILE)

    logging.info("Downloading documentation from RDKit")
    data = get_data(TMP_FILE)

    logging.info("Parsing documentation and combining with mypy's stubs")
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

    logging.info("Format and source output")
    os.system("isort " + out)
    os.system("black " + out)

    if os.path.exists("src/rdkit-stubs"):
        logging.info("Removing old stubs")
        shutil.rmtree("src/rdkit-stubs")

    # Move directory to src, so that it can be used as a package
    # Remove the temporary directory
    logging.info("Moving directory to src")
    shutil.move(out + "/rdkit", "src/rdkit-stubs")
    os.rmdir(out)


if __name__ == "__main__":
    main()
