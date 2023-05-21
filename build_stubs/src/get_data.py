# %%
import os
import pickle

import requests
from bs4 import BeautifulSoup
from bs4.element import ResultSet, Tag
from tqdm import tqdm

BASE_URL = "https://www.rdkit.org/docs/api-docs.html"
BASE_SOURCE = "https://www.rdkit.org/docs/source/"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.12; rv:55.0) Gecko/20100101 Firefox/55.0",
}


def get_base() -> BeautifulSoup:
    data = download_url(BASE_URL)
    node = BeautifulSoup(data, "html.parser")
    return node


def get_links(node: BeautifulSoup) -> list[str]:
    links: ResultSet[Tag] = node.find_all("a")
    # Only source links
    href = [link.get("href") for link in links if link.get("href").startswith("source")]
    links = [l.removeprefix("source/").split("#")[0] for l in href]
    links = list(set(links))

    return links


def download_url(url: str) -> str:
    item = requests.get(url, headers=HEADERS)
    item.encoding = "utf-8"
    return item.text


def download_data(links: list[str]) -> dict[str, str]:
    data = {}
    for link in tqdm(links):
        url = BASE_SOURCE + link
        data[link] = download_url(url)

    return data


def dump_data(data: dict[str, str], dump: str) -> None:
    pickle.dump(data, open(dump, "wb"))


def load_data(dump: str) -> dict[str, str]:
    return pickle.load(open(dump, "rb"))


def get_data(file: str) -> dict[str, str]:
    node = get_base()
    links = get_links(node)

    if not os.path.exists(file):
        data = download_data(links)
        dump_data(data, file)
    else:
        data = load_data(file)

    return data


# %%
