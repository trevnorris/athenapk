#!/usr/bin/env python3
"""Sanity check CUDA visibility from Python via CuPy."""

from __future__ import annotations

import json
import sys


def main() -> int:
    payload: dict[str, object] = {}
    try:
        import cupy as cp  # type: ignore
    except Exception as exc:  # pragma: no cover - probe utility
        payload["status"] = "FAIL"
        payload["error"] = f"cupy import failed: {exc}"
        print(json.dumps(payload))
        return 1

    payload["cupy_version"] = cp.__version__

    try:
        count = int(cp.cuda.runtime.getDeviceCount())
        payload["cuda_device_count"] = count
        if count <= 0:
            payload["status"] = "FAIL"
            payload["error"] = "No CUDA devices reported by CuPy runtime"
            print(json.dumps(payload))
            return 2

        props = cp.cuda.runtime.getDeviceProperties(0)
        name = props["name"]
        if isinstance(name, bytes):
            name = name.decode("utf-8", errors="replace")
        payload["device0_name"] = name
        payload["status"] = "PASS"
        print(json.dumps(payload))
        return 0
    except Exception as exc:  # pragma: no cover - probe utility
        payload["status"] = "FAIL"
        payload["error"] = str(exc)
        print(json.dumps(payload))
        return 3


if __name__ == "__main__":
    sys.exit(main())
