import http.server
import io
import os
import threading

import pytest

import opendisorder.cod as cod


class _FlakyHandler(http.server.BaseHTTPRequestHandler):
    """Serves DATA with Range support; the first two responses are cut short."""

    data = os.urandom(1_000_000)
    calls = 0

    def log_message(self, *args):
        pass

    def do_GET(self):
        type(self).calls += 1
        rng = self.headers.get("Range")
        start = int(rng.split("=")[1].split("-")[0]) if rng else 0
        self.send_response(206 if rng else 200)
        self.send_header("Content-Length", str(len(self.data) - start))
        self.end_headers()
        end = len(self.data) if self.calls > 2 else min(len(self.data), start + 300_000)
        self.wfile.write(self.data[start:end])
        self.wfile.flush()
        if end < len(self.data):
            self.connection.close()


@pytest.fixture()
def flaky_server(monkeypatch):
    monkeypatch.setenv("NO_PROXY", "127.0.0.1")
    monkeypatch.setattr(cod.time, "sleep", lambda s: None)
    _FlakyHandler.calls = 0
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _FlakyHandler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{srv.server_address[1]}/archive"
    srv.shutdown()


def test_resumable_reader(flaky_server):
    reader = io.BufferedReader(cod.ResumableHTTPReader(flaky_server))
    assert reader.read() == _FlakyHandler.data
    assert _FlakyHandler.calls == 3
