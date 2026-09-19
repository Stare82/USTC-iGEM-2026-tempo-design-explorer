"""Local HTTP server for the TEMPO Design Explorer and deterministic ODE API."""

from __future__ import annotations

import argparse
from functools import lru_cache
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import threading
import webbrowser

from backend.model_adapter import MODEL_VERSION, simulate_design


ROOT = Path(__file__).resolve().parent
SIMULATION_LOCK = threading.Lock()
REQUEST_SLOTS = threading.BoundedSemaphore(value=8)


@lru_cache(maxsize=48)
def _cached_simulation(normalized_json: str) -> dict:
    return simulate_design(json.loads(normalized_json))


class TempoHandler(SimpleHTTPRequestHandler):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(ROOT), **kwargs)

    def _send_json(self, payload: dict, status: HTTPStatus = HTTPStatus.OK) -> None:
        body = json.dumps(payload, ensure_ascii=False, allow_nan=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def end_headers(self) -> None:
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "strict-origin-when-cross-origin")
        super().end_headers()

    def do_GET(self) -> None:  # noqa: N802
        if self.path == "/api/health":
            self._send_json(
                {
                    "status": "ok",
                    "source": "deterministic_ode",
                    "model_version": MODEL_VERSION,
                    "states": {"oscillator": 8, "counter": 38, "shutdown": 3},
                }
            )
            return
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        if self.path != "/api/simulate":
            self._send_json({"error": "Not found"}, HTTPStatus.NOT_FOUND)
            return
        if not REQUEST_SLOTS.acquire(blocking=False):
            self._send_json(
                {"error": "The model queue is busy. Please retry shortly."},
                HTTPStatus.TOO_MANY_REQUESTS,
            )
            return
        try:
            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0 or length > 128_000:
                raise ValueError("Invalid request size")
            raw = json.loads(self.rfile.read(length).decode("utf-8"))
            normalized = json.dumps(raw, sort_keys=True, separators=(",", ":"))
            # SciPy's LSODA wrapper is not re-entrant. Browser slider requests
            # can overlap, so serialize model solves while retaining threaded
            # static-file serving.
            with SIMULATION_LOCK:
                result = _cached_simulation(normalized)
        except (ValueError, TypeError, json.JSONDecodeError) as error:
            self._send_json({"error": str(error)}, HTTPStatus.BAD_REQUEST)
            return
        except Exception as error:  # Keep the browser usable and surface a concise failure.
            self.log_error("Simulation failed: %s", error)
            self._send_json(
                {"error": "The deterministic ODE simulation failed for this design."},
                HTTPStatus.INTERNAL_SERVER_ERROR,
            )
            return
        finally:
            REQUEST_SLOTS.release()
        self._send_json(result)


class TempoServer(ThreadingHTTPServer):
    allow_reuse_address = True
    daemon_threads = True
    request_queue_size = 32


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default=os.environ.get("HOST", "127.0.0.1"))
    parser.add_argument("--port", type=int, default=int(os.environ.get("PORT", "4173")))
    parser.add_argument("--open", action="store_true", help="Open the default browser")
    args = parser.parse_args()

    server = TempoServer((args.host, args.port), TempoHandler)
    url = f"http://{args.host}:{args.port}/"
    print(f"TEMPO Design Explorer: {url}", flush=True)
    print(f"Model: {MODEL_VERSION} (8 + 38 + 3 states, SciPy LSODA)", flush=True)
    if args.open:
        threading.Timer(0.6, lambda: webbrowser.open(url)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
