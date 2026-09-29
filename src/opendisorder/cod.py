"""Thin client for the Crystallography Open Database (https://www.crystallography.net/).

Only a handful of endpoints are used:

* ``/cod/result?...&format=csv`` -- search / metadata dump.  ``id=%`` returns
  metadata for every entry (~530k rows, ~400 MB, about a minute).
* ``/cod/result`` POST with ``id=<comma separated>`` and ``format=zip`` -- a zip
  archive with up to 300 CIFs.
* ``/cod/<id>.cif`` -- a single CIF.

The COD asks bulk users to be gentle: requests are serialised and spaced by
``min_interval`` seconds, and failures are retried with exponential backoff.
"""

from __future__ import annotations

import io
import tarfile
import logging
import time
import zipfile
from collections.abc import Iterable, Iterator
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

COD_URL = "https://www.crystallography.net/cod"
#: Snapshot of all COD CIFs (tar + xz, ~17 GB; also holds a MySQL dump, which is skipped).
COD_ARCHIVE_URL = "https://www.crystallography.net/cod/archives/cod-cifs-mysql.txz"
USER_AGENT = "opendisorder/0.1 (+https://github.com/kazeevn/OpenDisorder)"
MAX_ZIP_BATCH = 300


class CODClient:
    def __init__(
        self,
        base_url: str = COD_URL,
        min_interval: float = 1.0,
        retries: int = 5,
        timeout: float = 600.0,
        session: requests.Session | None = None,
    ):
        self.base_url = base_url.rstrip("/")
        self.min_interval = min_interval
        self.retries = retries
        self.timeout = timeout
        self.session = session or requests.Session()
        self.session.headers.setdefault("User-Agent", USER_AGENT)
        self._last = 0.0

    # -- low level ----------------------------------------------------------

    def _request(self, method: str, url: str, **kwargs) -> requests.Response:
        delay = 2.0
        for attempt in range(self.retries):
            wait = self._last + self.min_interval - time.monotonic()
            if wait > 0:
                time.sleep(wait)
            try:
                resp = self.session.request(method, url, timeout=self.timeout, **kwargs)
                self._last = time.monotonic()
                if resp.status_code == 404:
                    resp.raise_for_status()
                if resp.status_code >= 500 or resp.status_code == 429:
                    raise requests.HTTPError(f"HTTP {resp.status_code}", response=resp)
                resp.raise_for_status()
                return resp
            except requests.HTTPError as exc:
                if exc.response is not None and exc.response.status_code == 404:
                    raise
                err = exc
            except (requests.ConnectionError, requests.Timeout) as exc:
                err = exc
            self._last = time.monotonic()
            logger.warning("COD request failed (%s), attempt %d/%d", err, attempt + 1, self.retries)
            time.sleep(delay)
            delay *= 2
        raise err  # noqa: F821 - defined in the loop

    # -- search -------------------------------------------------------------

    def search(self, **params) -> list[dict]:
        """Run a COD search and return the JSON records (all values are strings).

        Useful parameters: ``formula`` (Hill order, space separated, exact
        match, e.g. ``"Fe Li O4 P"``), ``el1``..``el8``, ``nel1``..``nel4``,
        ``strictmin``/``strictmax`` (number of distinct elements), ``id``.
        """
        resp = self._request("GET", f"{self.base_url}/result", params={**params, "format": "json"})
        return resp.json()

    def download_metadata(self, path: str | Path, **params) -> Path:
        """Stream the CSV metadata of all entries (or of a search) to ``path``."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        query = {"id": "%", **params, "format": "csv"}
        tmp = path.with_suffix(path.suffix + ".part")
        resp = self._request("GET", f"{self.base_url}/result", params=query, stream=True)
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(1 << 20):
                fh.write(chunk)
        tmp.replace(path)
        return path

    # -- CIFs ---------------------------------------------------------------

    def get_cif(self, cod_id: str | int) -> str:
        return self._request("GET", f"{self.base_url}/{cod_id}.cif").text

    def iter_cifs(self, cod_ids: Iterable[str | int], batch_size: int = MAX_ZIP_BATCH) -> Iterator[tuple[str, str]]:
        """Yield ``(cod_id, cif_text)`` pairs, fetching in zip batches.

        A batch that keeps failing is split in halves; single entries are
        fetched as plain CIFs.  Entries that cannot be downloaded are skipped
        with a warning.
        """
        ids = [str(i) for i in cod_ids]
        batch_size = min(batch_size, MAX_ZIP_BATCH)
        for start in range(0, len(ids), batch_size):
            yield from self._fetch_batch(ids[start : start + batch_size])

    def _fetch_batch(self, batch: list[str]) -> Iterator[tuple[str, str]]:
        if len(batch) == 1:
            try:
                yield batch[0], self.get_cif(batch[0])
            except requests.RequestException as exc:
                logger.warning("could not download COD %s: %s", batch[0], exc)
            return
        retries, self.retries = self.retries, min(self.retries, 3)
        try:
            resp = self._request("POST", f"{self.base_url}/result", data={"id": ",".join(batch), "format": "zip"})
            zf = zipfile.ZipFile(io.BytesIO(resp.content))
        except (requests.RequestException, zipfile.BadZipFile) as exc:
            logger.warning("zip batch of %d failed (%s); splitting", len(batch), exc)
            half = len(batch) // 2
            yield from self._fetch_batch(batch[:half])
            yield from self._fetch_batch(batch[half:])
            return
        finally:
            self.retries = retries
        for name in zf.namelist():
            if name.endswith(".cif"):
                yield Path(name).stem, zf.read(name).decode("utf-8", errors="replace")


def iter_archive_cifs(source: str | Path = COD_ARCHIVE_URL, session: requests.Session | None = None) -> Iterator[tuple[str, str]]:
    """Stream ``(cod_id, cif_text)`` from a COD CIF archive without storing it.

    ``source`` is a URL or a local path to one of the archives in
    https://www.crystallography.net/cod/archives/ (``.txz``, ``.tgz``) or any
    tar archive whose members are named ``.../<cod_id>.cif``.
    """
    source = str(source)
    if source.startswith(("http://", "https://")):
        session = session or requests.Session()
        session.headers.setdefault("User-Agent", USER_AGENT)
        resp = session.get(source, stream=True, timeout=600)
        resp.raise_for_status()
        resp.raw.decode_content = True
        fileobj, close = resp.raw, resp.close
    else:
        fileobj = open(source, "rb")
        close = fileobj.close
    try:
        with tarfile.open(fileobj=fileobj, mode="r|*") as tf:
            for member in tf:
                if not member.isfile() or not member.name.endswith(".cif"):
                    continue
                fh = tf.extractfile(member)
                if fh is None:
                    continue
                yield Path(member.name).stem, fh.read().decode("utf-8", errors="replace")
    finally:
        close()
