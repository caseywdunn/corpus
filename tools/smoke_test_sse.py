#!/usr/bin/env python3
"""Smoke-test the MCP server's SSE transport + bearer-token auth.

Runs two layers against a locally-launched mcp_server.py:

  Layer 1 — raw HTTP: unauthenticated → 401, wrong-token → 401,
            correct-token → 200 with ``Content-Type: text/event-stream``.
            Catches middleware / uvicorn / binding issues.

  Layer 2 — real MCP client: initialize handshake, ``list_tools``,
            ``call_tool(bundle_info)``, ``call_tool(list_papers,
            limit=1)``.  Catches framing, keepalives, and tool-
            dispatch over SSE that stdio masks.

Starts the server as a subprocess, waits for it to bind, runs the
checks, shuts it down cleanly.  No external deps — stdlib for layer
1, the ``mcp`` Python SDK (already a dep) for layer 2.

Usage:
    python tools/smoke_test_sse.py output
    python tools/smoke_test_sse.py output --port 18080 --skip-layer2
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import secrets
import socket
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SERVER_PY = REPO_ROOT / "mcp_server.py"


# ── Tiny TAP-ish reporter ───────────────────────────────────────────

_USE_COLOR = sys.stdout.isatty()


def _color(s: str, code: str) -> str:
    return f"\033[{code}m{s}\033[0m" if _USE_COLOR else s


def _ok(msg: str) -> int:
    print(f"  {_color('PASS', '32')}  {msg}")
    return 0


def _fail(msg: str) -> int:
    print(f"  {_color('FAIL', '31')}  {msg}")
    return 1


def _info(msg: str) -> None:
    print(f"  {_color('•', '36')}     {msg}")


def _section(title: str) -> None:
    print(f"\n{_color('══', '34')} {title}")


# ── Server lifecycle ────────────────────────────────────────────────

def _wait_for_bind(host: str, port: int, timeout: float = 15.0) -> bool:
    deadline = time.time() + timeout
    while time.time() < deadline:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                return True
        except OSError:
            time.sleep(0.2)
    return False


# ── Layer 1: raw HTTP ───────────────────────────────────────────────

def layer1_http(host: str, port: int, token: str) -> int:
    _section("Layer 1: raw HTTP + bearer auth")
    url = f"http://{host}:{port}/sse"
    rc = 0

    # 1a: unauthenticated GET → 401
    try:
        urllib.request.urlopen(url, timeout=3)
        rc |= _fail("unauth GET /sse did not raise")
    except urllib.error.HTTPError as e:
        if e.code == 401:
            _ok(f"unauth GET /sse → 401")
        else:
            rc |= _fail(f"unauth GET /sse → {e.code} (expected 401)")
    except Exception as e:
        rc |= _fail(f"unauth GET /sse raised {type(e).__name__}: {e}")

    # 1b: wrong-token GET → 401
    req = urllib.request.Request(url, headers={"Authorization": "Bearer nope"})
    try:
        urllib.request.urlopen(req, timeout=3)
        rc |= _fail("wrong-token GET did not raise")
    except urllib.error.HTTPError as e:
        if e.code == 401:
            _ok("wrong-token GET /sse → 401")
        else:
            rc |= _fail(f"wrong-token GET /sse → {e.code} (expected 401)")

    # 1c: authenticated GET → 200 + text/event-stream Content-Type
    # SSE is an open-ended stream; we only need to confirm the response
    # line + headers arrive cleanly, not drain the body.
    req = urllib.request.Request(
        url, headers={"Authorization": f"Bearer {token}"}
    )
    try:
        resp = urllib.request.urlopen(req, timeout=3)
    except urllib.error.HTTPError as e:
        rc |= _fail(f"auth GET /sse → {e.code}")
        return rc
    try:
        ctype = resp.headers.get("Content-Type", "")
        if resp.status == 200 and "text/event-stream" in ctype:
            _ok(f"auth GET /sse → 200 {ctype}")
        else:
            rc |= _fail(
                f"auth GET /sse → status={resp.status} ctype={ctype!r}"
            )
    finally:
        resp.close()

    return rc


# ── Layer 2: MCP client over SSE ────────────────────────────────────

async def layer2_mcp_client(host: str, port: int, token: str) -> int:
    _section("Layer 2: MCP client over SSE")
    try:
        from mcp import ClientSession
        from mcp.client.sse import sse_client
    except ImportError as e:
        _info(f"mcp SDK client not importable: {e}")
        _info("skipping layer 2 (install `mcp` or run with --skip-layer2)")
        return 0

    url = f"http://{host}:{port}/sse"
    headers = {"Authorization": f"Bearer {token}"}
    rc = 0

    try:
        async with sse_client(url, headers=headers) as (read, write):
            async with ClientSession(read, write) as session:
                await session.initialize()
                _ok("initialize handshake")

                # list_tools
                try:
                    listing = await session.list_tools()
                    tools = listing.tools
                    if len(tools) < 1:
                        rc |= _fail(f"list_tools returned {len(tools)} tools")
                    else:
                        names = sorted(t.name for t in tools)
                        _ok(f"list_tools → {len(tools)} tools "
                            f"(e.g. {', '.join(names[:3])}, …)")
                except Exception as e:
                    rc |= _fail(f"list_tools raised: {e}")

                # call bundle_info — returns a dict regardless of deploy
                # state.  Parse the first text block as JSON.
                try:
                    result = await session.call_tool("bundle_info", {})
                    parsed = _parse_tool_result(result)
                    if parsed is not None and "bundle_version" in parsed:
                        _ok(
                            f"bundle_info → bundle_version="
                            f"{parsed.get('bundle_version')!r}"
                        )
                    else:
                        rc |= _fail(
                            f"bundle_info result missing bundle_version: "
                            f"{parsed!r}"
                        )
                except Exception as e:
                    rc |= _fail(f"call_tool(bundle_info) raised: {e}")

                # call list_papers(limit=1) — cheap and should succeed
                # on any processed corpus.
                try:
                    result = await session.call_tool(
                        "list_papers", {"limit": 1}
                    )
                    parsed = _parse_tool_result(result)
                    if isinstance(parsed, list):
                        _ok(f"list_papers(limit=1) → list of {len(parsed)}")
                    else:
                        _ok(
                            f"list_papers(limit=1) returned "
                            f"{type(parsed).__name__}"
                        )
                except Exception as e:
                    rc |= _fail(f"call_tool(list_papers) raised: {e}")

    except Exception as e:
        rc |= _fail(f"SSE session setup failed: {type(e).__name__}: {e}")

    return rc


def _parse_tool_result(result):
    """Best-effort extraction of a JSON payload from an MCP tool
    result's content blocks.  Returns the parsed value or None."""
    for block in getattr(result, "content", []) or []:
        text = getattr(block, "text", None)
        if not text:
            continue
        try:
            return json.loads(text)
        except Exception:
            # Some tools return plain text — hand that back as-is.
            return text
    return None


# ── Orchestration ───────────────────────────────────────────────────

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("output_dir", type=Path,
                        help="Processed corpus output dir to serve")
    parser.add_argument("--host", default="127.0.0.1",
                        help="Bind interface (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=18080,
                        help="Bind port (default: 18080, to avoid clashing "
                             "with a real 8080 service you may be running)")
    parser.add_argument("--server-python", default=sys.executable,
                        help="Python executable for the mcp_server subprocess")
    parser.add_argument("--skip-layer2", action="store_true",
                        help="Skip the MCP-client layer (useful when the "
                             "mcp SDK isn't installed in the current env)")
    parser.add_argument("--server-startup-timeout", type=float, default=15.0,
                        help="Seconds to wait for the server to bind (default: 15)")
    args = parser.parse_args()

    token = secrets.token_hex(16)
    tok_file = tempfile.NamedTemporaryFile(
        mode="w", delete=False, suffix=".token"
    )
    try:
        tok_file.write(token)
        tok_file.close()

        stderr_log = tempfile.NamedTemporaryFile(
            mode="w+", delete=False, suffix=".stderr"
        )
        stderr_log.close()

        cmd = [
            args.server_python, str(SERVER_PY), str(args.output_dir),
            "--transport", "sse", "--host", args.host, "--port", str(args.port),
            "--auth-token-file", tok_file.name,
        ]
        _section("Launching mcp_server.py --transport sse")
        _info(f"cmd: {' '.join(cmd)}")
        _info(f"stderr -> {stderr_log.name}")

        # stderr goes to a file so post-mortem diagnostics are readable
        # without blocking on a pipe.
        with open(stderr_log.name, "w") as errfh:
            server = subprocess.Popen(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=errfh,
            )

        rc = 0
        try:
            if not _wait_for_bind(args.host, args.port,
                                  args.server_startup_timeout):
                _fail(
                    f"server did not bind on {args.host}:{args.port} "
                    f"within {args.server_startup_timeout}s"
                )
                tail = Path(stderr_log.name).read_text()[-2000:]
                print(f"\nstderr tail:\n{tail}")
                return 1
            _ok(f"server bound on {args.host}:{args.port}")

            rc |= layer1_http(args.host, args.port, token)
            if not args.skip_layer2:
                rc |= asyncio.run(
                    layer2_mcp_client(args.host, args.port, token)
                )

        finally:
            server.terminate()
            try:
                server.wait(timeout=5)
            except subprocess.TimeoutExpired:
                server.kill()
                server.wait()

        _section("Result")
        if rc == 0:
            print(f"  {_color('All layers passed.', '32')}")
        else:
            print(f"  {_color('Some layers FAILED.', '31')}")
            print(f"  stderr log: {stderr_log.name}")
        return rc
    finally:
        os.unlink(tok_file.name)


if __name__ == "__main__":
    sys.exit(main())
