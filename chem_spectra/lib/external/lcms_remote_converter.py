import os
from typing import Optional
from urllib.parse import urljoin, urlparse

import requests
from flask import current_app, has_app_context
from requests.auth import HTTPBasicAuth


class RemoteConverterError(RuntimeError):
    def __init__(self, message: str, status_code: int = 502):
        super().__init__(message)
        self.status_code = status_code
        self.message = message


def _env(primary: str, fallback: Optional[str] = None, default: str = "") -> str:
    value = os.getenv(primary)
    if value:
        return value
    if fallback:
        value = os.getenv(fallback)
        if value:
            return value
    return default


def _config(primary: str, fallback: Optional[str] = None, default: str = "") -> str:
    if not has_app_context():
        return default
    value = current_app.config.get(primary)
    if value:
        return str(value)
    if fallback:
        value = current_app.config.get(fallback)
        if value:
            return str(value)
    return default


def _setting(primary: str, fallback: Optional[str] = None, default: str = "") -> str:
    config_value = _config(primary, fallback, "")
    if config_value:
        return config_value
    return _env(primary, fallback, default)


def _to_float(raw: str, default: float) -> float:
    try:
        return float(raw)
    except (TypeError, ValueError):
        return default


def _converter_url(base_url_override: Optional[str] = None) -> str:
    configured = _setting("CONVERTER_BASE_URL", default="http://193.196.39.140:4000").strip()
    raw = (base_url_override or "").strip() or configured
    parsed = urlparse(raw)
    if parsed.hostname == "host.docker.internal":
        raw = configured
        parsed = urlparse(raw)
    default_path = "/conversions"
    path = _setting("CONVERTER_PATH", default=default_path).strip() or default_path
    return urljoin(raw.rstrip("/") + "/", path.lstrip("/"))


def _multipart_content_type(source_path: str) -> str:
    lowered = source_path.lower()
    if lowered.endswith((".tar.gz", ".tgz", ".gz")):
        return "application/x-gzip"
    if lowered.endswith(".tar"):
        return "application/x-tar"
    if lowered.endswith(".zip"):
        return "application/zip"
    return "application/octet-stream"


def convert_file_to_jcampzip(
    source_path: str,
    converter_url: Optional[str] = None,
) -> bytes:
    if not source_path or not os.path.exists(source_path) or os.path.isdir(source_path):
        raise RemoteConverterError("Invalid source file for remote conversion.", status_code=400)

    username = _setting("CONVERTER_USER", default="chemotion").strip() or "chemotion"
    token = _setting("CONVERTER_TOKEN", default="chemotion").strip()
    timeout_raw = _setting("CONVERTER_TIMEOUT_SEC", default="180")
    timeout_sec = _to_float(timeout_raw, 180.0)
    endpoint = _converter_url(converter_url)
    auth = HTTPBasicAuth(username, token) if token else None

    try:
        with open(source_path, "rb") as handle:
            filename = os.path.basename(source_path) or "lcms_input.raw"
            content_type = _multipart_content_type(source_path)
            files = {"file": (filename, handle, content_type)}
            data = {"format": "jcampzip"}
            response = requests.post(
                endpoint,
                files=files,
                data=data,
                auth=auth,
                timeout=timeout_sec,
            )
    except requests.Timeout as exc:
        raise RemoteConverterError(
            "Timeout while calling remote converter.",
            status_code=504,
        ) from exc
    except requests.RequestException as exc:
        raise RemoteConverterError(
            "Network error while calling remote converter.",
            status_code=502,
        ) from exc

    if response.status_code != 200:
        detail = ""
        try:
            body = (response.text or "").strip()
            if body:
                detail = f" Body: {body[:800]}"
        except Exception:
            pass
        if 400 <= response.status_code < 500:
            raise RemoteConverterError(
                f"Converter returned a client error ({response.status_code}).{detail}",
                status_code=422,
            )
        raise RemoteConverterError(
            f"Converter returned a server error ({response.status_code}).{detail}",
            status_code=502,
        )

    if not response.content:
        raise RemoteConverterError("Remote converter returned an empty response.", status_code=502)

    return response.content


# Backward compatibility
LcmsRemoteConverterError = RemoteConverterError
