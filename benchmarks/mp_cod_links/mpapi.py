"""Tiny MP API helper. Reads the key from $MP_API_KEY; never prints it."""
import json, os, sys, urllib.parse, urllib.request

_KEY = os.environ["MP_API_KEY"]
BASE = "https://api.materialsproject.org"


def get(path, **params):
    url = BASE + path + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"X-API-KEY": _KEY, "accept": "application/json", "User-Agent": "opendisorder-research/0.1 python-requests"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        return {"http_error": e.code, "body": e.read().decode()[:500]}


if __name__ == "__main__":
    path = sys.argv[1]
    params = dict(a.split("=", 1) for a in sys.argv[2:])
    print(json.dumps(get(path, **params), indent=1)[:8000])
