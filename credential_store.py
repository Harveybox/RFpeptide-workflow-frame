#!/usr/bin/env python3
import base64
import ctypes
import json
import sys
from ctypes import wintypes
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_STORE = SCRIPT_DIR / "credentials.local.json"


class CredentialStoreError(RuntimeError):
    pass


class DATA_BLOB(ctypes.Structure):
    _fields_ = [
        ("cbData", wintypes.DWORD),
        ("pbData", ctypes.POINTER(ctypes.c_char)),
    ]


def is_available() -> bool:
    return sys.platform.startswith("win")


def _require_windows():
    if not is_available():
        raise CredentialStoreError("Secure password storage currently requires Windows DPAPI.")


def _blob_from_bytes(data: bytes) -> DATA_BLOB:
    buffer = ctypes.create_string_buffer(data)
    blob = DATA_BLOB(len(data), ctypes.cast(buffer, ctypes.POINTER(ctypes.c_char)))
    blob._buffer = buffer
    return blob


def _bytes_from_blob(blob: DATA_BLOB) -> bytes:
    try:
        return ctypes.string_at(blob.pbData, blob.cbData)
    finally:
        ctypes.windll.kernel32.LocalFree(blob.pbData)


def protect_text(text: str) -> str:
    _require_windows()
    input_blob = _blob_from_bytes(text.encode("utf-8"))
    output_blob = DATA_BLOB()
    ok = ctypes.windll.crypt32.CryptProtectData(
        ctypes.byref(input_blob),
        None,
        None,
        None,
        None,
        0,
        ctypes.byref(output_blob),
    )
    if not ok:
        raise CredentialStoreError("CryptProtectData failed.")
    return base64.b64encode(_bytes_from_blob(output_blob)).decode("ascii")


def unprotect_text(encoded: str) -> str:
    _require_windows()
    encrypted = base64.b64decode(encoded.encode("ascii"))
    input_blob = _blob_from_bytes(encrypted)
    output_blob = DATA_BLOB()
    ok = ctypes.windll.crypt32.CryptUnprotectData(
        ctypes.byref(input_blob),
        None,
        None,
        None,
        None,
        0,
        ctypes.byref(output_blob),
    )
    if not ok:
        raise CredentialStoreError("CryptUnprotectData failed.")
    return _bytes_from_blob(output_blob).decode("utf-8")


def load_store(path: Path = DEFAULT_STORE) -> dict:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def save_store(store: dict, path: Path = DEFAULT_STORE) -> None:
    path.write_text(json.dumps(store, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def make_key(profile: dict) -> str:
    host = str(profile.get("host") or "").strip()
    port = str(profile.get("port") or 22).strip()
    user = str(profile.get("user") or "").strip()
    if not host or not user:
        raise CredentialStoreError("Credential key requires host and user.")
    return f"{user}@{host}:{port}"


def save_password(profile: dict, password: str, path: Path = DEFAULT_STORE) -> str:
    if not password:
        raise CredentialStoreError("Password is empty.")
    key = make_key(profile)
    store = load_store(path)
    store[key] = {
        "backend": "windows-dpapi",
        "password": protect_text(password),
    }
    save_store(store, path)
    return key


def load_password(profile: dict, path: Path = DEFAULT_STORE) -> str | None:
    key = make_key(profile)
    store = load_store(path)
    entry = store.get(key)
    if not entry:
        return None
    if entry.get("backend") != "windows-dpapi":
        raise CredentialStoreError(f"Unsupported credential backend for {key}.")
    return unprotect_text(entry["password"])


def delete_password(profile: dict, path: Path = DEFAULT_STORE) -> bool:
    key = make_key(profile)
    store = load_store(path)
    existed = key in store
    store.pop(key, None)
    save_store(store, path)
    return existed
