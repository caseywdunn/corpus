"""Expiring figure-only capabilities, independent of the MCP bearer secret."""
import hashlib
import hmac
import json
import re
import secrets
import threading
import time
from urllib.parse import parse_qs, urlencode, urlsplit

_PATH = re.compile(r"/figures/[a-f0-9]{12}/[A-Za-z0-9._-]{1,256}\Z")
_LOCK = threading.Lock()
URL_LIFETIME_SECONDS = 300


def validate_public_base(value):
    parts = urlsplit(value)
    if (parts.scheme not in {"http", "https"} or not parts.hostname
            or parts.username or parts.password or parts.query or parts.fragment
            or any(ord(c) <= 32 for c in value)
            or parts.hostname in {"0.0.0.0", "::"}):
        raise ValueError("public base URL must be an HTTP(S) client-reachable URL without credentials, query or fragment")
    # Validate a malformed port here, rather than at the first download.
    parts.port
    return value.rstrip("/")


class FigureURLSigner:
    def __init__(self, *, key=None, clock=time.time):
        self._key = key if key is not None else secrets.token_bytes(32)
        self._clock = clock

    def _signature(self, path, params):
        payload = json.dumps([path, params], sort_keys=True, separators=(",", ":")).encode()
        return hmac.new(self._key, payload, hashlib.sha256).hexdigest()

    def url(self, base, path, profile, label=None):
        if not _PATH.fullmatch(path):
            raise ValueError("invalid figure download path")
        if label is not None and not re.fullmatch(r"[A-Za-z0-9._-]{1,128}", label):
            raise ValueError("invalid figure label")
        params = {"profile": profile, "expires": str(int(self._clock()) + URL_LIFETIME_SECONDS)}
        if label is not None:
            params["label"] = label
        params["signature"] = self._signature(path, params)
        return validate_public_base(base) + path + "?" + urlencode(params)

    def accepts(self, scope):
        if scope.get("type") != "http" or scope.get("method") not in {"GET", "HEAD"}:
            return False
        path = scope.get("path", "")
        query = scope.get("query_string", b"")
        if not _PATH.fullmatch(path) or len(query) > 2048:
            return False
        try:
            parsed = parse_qs(query.decode("ascii"), keep_blank_values=True, max_num_fields=4)
            if any(len(values) != 1 for values in parsed.values()):
                return False
            params = {key: values[0] for key, values in parsed.items()}
            if set(params) not in ({"profile", "expires", "signature"}, {"profile", "expires", "signature", "label"}):
                return False
            signature = params.pop("signature")
            remaining = int(params["expires"]) - self._clock()
            if not 0 < remaining <= URL_LIFETIME_SECONDS:
                return False
            return hmac.compare_digest(signature.encode("ascii"), self._signature(path, params).encode("ascii"))
        except (ValueError, TypeError, UnicodeError):
            return False


def figure_signer(idx):
    with _LOCK:
        signer = getattr(idx, "_figure_url_signer", None)
        if signer is None:
            signer = FigureURLSigner()
            idx._figure_url_signer = signer
    return signer
