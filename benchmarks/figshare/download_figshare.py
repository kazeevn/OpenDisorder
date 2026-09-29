"""Download figshare article 29145101 (generated crystals) into ``<dest>/data``.

All files are named ``data.csv.gz``; the directory tree is rebuilt from the
article's ``folder_structure`` field.  md5 checksums are verified.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

import requests

ARTICLE = "https://api.figshare.com/v2/articles/29145101"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("dest", type=Path, help="destination directory; files go to <dest>/data")
    args = ap.parse_args(argv)
    session = requests.Session()
    article = session.get(ARTICLE, timeout=60).json()
    folders = article.get("folder_structure", {})
    for f in article["files"]:
        target = args.dest / "data" / folders.get(str(f["id"]), "") / f["name"]
        md5 = f.get("supplied_md5") or f.get("computed_md5")
        if target.exists() and hashlib.md5(target.read_bytes()).hexdigest() == md5:
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        with session.get(f["download_url"], stream=True, timeout=600) as resp:
            resp.raise_for_status()
            with open(target, "wb") as fh:
                for chunk in resp.iter_content(1 << 20):
                    fh.write(chunk)
        if md5 and hashlib.md5(target.read_bytes()).hexdigest() != md5:
            raise RuntimeError(f"md5 mismatch for {target}")
        print(target)


if __name__ == "__main__":
    main()
