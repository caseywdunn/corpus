"""Figure-only authorization must not leak or replace the MCP bearer secret."""
import asyncio
import io
import json
import os
from urllib.parse import parse_qs, urlencode, urlsplit

import httpx
import pytest
from PIL import Image

from mcpsrv.figure_http import make_figure_app
from mcpsrv.figure_urls import FigureURLSigner, validate_public_base
from mcpsrv.transport import _BearerAuthASGI, _RouteMuxASGI
from mcpsrv.tools.figures import get_figure_url
from tests.test_figure_licensing_states import _make_index, _CLEARED, _UNCLEARED, HASH


@pytest.fixture
def server(tmp_path):
    idx = _make_index(tmp_path, work=dict(_CLEARED))
    idx.figure_url_base = "https://corpus.example.test"
    idx.figure_auth_token = "shared-mcp-secret-never-disclose"
    clock = [1000.0]
    signer = FigureURLSigner(clock=lambda: clock[0])
    idx._figure_url_signer = signer
    async def sse(scope, receive, send):
        await send({"type": "http.response.start", "status": 200, "headers": []})
        await send({"type": "http.response.body", "body": b"MCP"})
    app = _BearerAuthASGI(_RouteMuxASGI(sse, make_figure_app(idx)), idx.figure_auth_token, figure_signer=signer)
    return idx, clock, app


def request(app, url, method="GET", **kwargs):
    async def call():
        async with httpx.AsyncClient(transport=httpx.ASGITransport(app=app), base_url="https://corpus.example.test") as client:
            return await client.request(method, url, **kwargs)
    return asyncio.run(call())


def test_signed_download_without_bearer_and_frozen_result_fields(server):
    idx, _, app = server
    result = get_figure_url(HASH, "docling_1", "A", profile="manuscript")
    assert result["auth_header"] is None
    assert idx.figure_auth_token not in json.dumps(result)
    response = request(app, result["url"])
    assert response.status_code == 200
    assert Image.open(io.BytesIO(response.content)).size == (20, 20)
    assert response.headers["cache-control"] == "private, no-store"
    assert request(app, result["url"], "HEAD").content == b""


@pytest.mark.parametrize("change", ["profile", "label", "path", "expires", "duplicate", "extra", "method"])
def test_tampered_capability_cannot_change_its_scope(server, change):
    _, _, app = server
    url = get_figure_url(HASH, "docling_1", "A", profile="report")["url"]
    parts = urlsplit(url)
    params = {k: v[0] for k, v in parse_qs(parts.query).items()}
    path, method = parts.path, "GET"
    if change in {"profile", "label", "expires"}:
        params[change] = {"profile": "manuscript", "label": "B", "expires": "999999"}[change]
    elif change == "path":
        path = path.replace("docling_1", "docling_2")
    elif change == "extra":
        params["other"] = "ignored"
    elif change == "method":
        method = "POST"
    query = urlencode(params) + ("&profile=manuscript" if change == "duplicate" else "")
    assert request(app, path + "?" + query, method).status_code == 401


def test_expired_and_previous_process_links_fail(server):
    idx, clock, app = server
    url = get_figure_url(HASH, "docling_1")["url"]
    clock[0] += 301
    assert request(app, url).status_code == 401
    clock[0] = 1000
    replacement = _BearerAuthASGI(make_figure_app(idx), idx.figure_auth_token, figure_signer=FigureURLSigner(clock=lambda: 1000))
    assert request(replacement, url).status_code == 401


def test_figure_signature_never_authorizes_mcp_routes(server):
    idx, _, app = server
    url = get_figure_url(HASH, "docling_1")["url"]
    query = urlsplit(url).query
    for route in ("/sse", "/messages/session"):
        assert request(app, route + "?" + query).status_code == 401
        assert request(app, route, headers={"Authorization": "Bearer " + idx.figure_auth_token}).status_code == 200


def test_download_rechecks_licensing_at_boundary(server):
    idx, _, app = server
    url = get_figure_url(HASH, "docling_1", profile="manuscript")["url"]
    idx.biblio_db._work = dict(_UNCLEARED)
    assert request(app, url).status_code == 403
    assert get_figure_url(HASH, "docling_1", profile="manuscript")["code"] == "forbidden"


def test_public_proxy_preserves_signed_path_and_query(server):
    _, _, upstream = server
    async def proxy(scope, receive, send):
        # Same relevant contract as deploy/nginx.conf: forward /figures/
        # unchanged to a fixed loopback upstream, never a caller-chosen host.
        if not scope["path"].startswith("/figures/"):
            await send({"type": "http.response.start", "status": 404, "headers": []})
            await send({"type": "http.response.body", "body": b""})
            return
        forwarded = dict(scope, headers=[(b"host", b"localhost:8080")])
        await upstream(forwarded, receive, send)
    for label, size in ((None, (40, 40)), ("A", (20, 20))):
        response = request(proxy, get_figure_url(HASH, "docling_1", label)["url"])
        assert response.status_code == 200
        assert Image.open(io.BytesIO(response.content)).size == size


@pytest.mark.parametrize("url", ["https://user:password@example.org", "http://0.0.0.0:8080", "https://example.org?token=x", "file:///tmp/corpus", "https://example.org/#fragment", "https://example.org:bad", "https://example.org\n/evil"])
def test_public_base_validation_rejects_unsafe_or_unusable_values(url):
    with pytest.raises(ValueError):
        validate_public_base(url)


def test_deployed_proxy_exposes_downloads_without_logging_capabilities():
    from pathlib import Path
    text = (Path(__file__).parents[1] / "deploy" / "nginx.conf").read_text()
    route = text.split("location /figures/ {", 1)[1].split("}", 1)[0]
    assert "proxy_pass http://127.0.0.1:8080;" in route
    assert "access_log off;" in route
    assert "proxy_cache off;" in route


@pytest.mark.skipif(os.environ.get("CORPUS_TEST_NGINX") != "1", reason="opt-in live Docker/nginx acceptance")
def test_live_nginx_downloads(server, tmp_path):
    """Exercise the deployed route through real TCP/nginx, not an ASGI mock.

    Run with CORPUS_TEST_NGINX=1. Requires Docker and a locally available
    nginx:stable-alpine image; uses loopback high ports, no existing services.
    """
    import socket
    import subprocess
    import threading
    import time
    import uuid
    from pathlib import Path
    import uvicorn

    idx, clock, app = server
    upstream = socket.socket()
    upstream.bind(("127.0.0.1", 0))
    upstream.listen()
    upstream_port = upstream.getsockname()[1]
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        proxy_port = probe.getsockname()[1]
    config = uvicorn.Config(app, log_level="error", access_log=False)
    uvicorn_server = uvicorn.Server(config)
    worker = threading.Thread(target=uvicorn_server.run, kwargs={"sockets": [upstream]}, daemon=True)
    worker.start()
    deployed = (Path(__file__).parents[1] / "deploy/nginx.conf").read_text()
    route = "location /figures/ {" + deployed.split("location /figures/ {", 1)[1].split("}", 1)[0] + "}"
    route = route.replace("127.0.0.1:8080", f"127.0.0.1:{upstream_port}")
    nginx_config = tmp_path / "nginx.conf"
    nginx_config.write_text(
        "pid /tmp/nginx.pid; error_log /dev/stderr warn; events {} http { "
        "access_log off; client_body_temp_path /tmp/client; proxy_temp_path /tmp/proxy; "
        "fastcgi_temp_path /tmp/fastcgi; uwsgi_temp_path /tmp/uwsgi; scgi_temp_path /tmp/scgi; "
        f"server {{ listen 127.0.0.1:{proxy_port}; {route} }} }}"
    )
    name = "corpus-figure-test-" + uuid.uuid4().hex
    container_started = False
    try:
        subprocess.run([
            "docker", "run", "--pull=never", "--rm", "-d", "--name", name,
            "--network=host", "--read-only", "--cap-drop=ALL", "--tmpfs", "/tmp",
            "--user", f"{os.getuid()}:{os.getgid()}", "--entrypoint", "nginx",
            "-v", f"{nginx_config}:/etc/nginx/nginx.conf:ro", "nginx:stable-alpine",
            "-e", "stderr", "-g", "daemon off;",
        ], check=True, capture_output=True)
        container_started = True
        idx.figure_url_base = f"http://127.0.0.1:{proxy_port}"
        deadline = time.monotonic() + 10
        while True:
            try:
                with socket.create_connection(("127.0.0.1", proxy_port), timeout=0.2):
                    break
            except OSError:
                if time.monotonic() >= deadline:
                    pytest.fail("nginx did not bind its test port")
                time.sleep(0.05)
        with httpx.Client(trust_env=False) as client:
            for label, size in ((None, (40, 40)), ("A", (20, 20))):
                url = get_figure_url(HASH, "docling_1", label, profile="manuscript")["url"]
                response = client.get(url)
                assert response.status_code == 200, response.text
                assert Image.open(io.BytesIO(response.content)).size == size
                assert client.head(url).status_code == 200
                assert client.get(url.replace("profile=manuscript", "profile=report")).status_code == 401
            clock[0] += 301
            assert client.get(url).status_code == 401
    finally:
        if container_started:
            # --rm may already have removed a failed startup; don't mask
            # the original failure with a second cleanup exception.
            subprocess.run(["docker", "stop", "-t", "1", name], capture_output=True)
        uvicorn_server.should_exit = True
        worker.join(timeout=10)
        upstream.close()
