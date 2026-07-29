#!/usr/bin/env python3
"""Headless test of the FUSE container Jupyter kernel (no browser/tunnel).

Starts the kernelspec installed by install_kernel.sh through jupyter_client,
runs a FUSE smoke cell, and checks the output. Usage:

    python3 test_kernel_headless.py fuse-<version>

Requires jupyter_client (available in the omega FUSE conda environment).
"""

import sys

from jupyter_client.manager import KernelManager

CODE = (
    'using FUSE; '
    'ini, act = FUSE.case_parameters(:D3D, :L_mode); '
    'println("KERNEL OK ", pkgversion(FUSE))'
)


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    kernel_name = sys.argv[1]

    km = KernelManager(kernel_name=kernel_name)
    km.start_kernel()
    kc = km.client()
    kc.start_channels()
    try:
        print(f"waiting for kernel '{kernel_name}' to become ready ...")
        kc.wait_for_ready(timeout=300)
        print("kernel ready, executing smoke cell ...")
        msg_id = kc.execute(CODE)
        output = []
        while True:
            msg = kc.get_iopub_msg(timeout=600)
            if msg["parent_header"].get("msg_id") != msg_id:
                continue
            msg_type = msg["msg_type"]
            if msg_type == "stream":
                output.append(msg["content"]["text"])
            elif msg_type == "error":
                print("\n".join(msg["content"]["traceback"]))
                return 1
            elif msg_type == "status" and msg["content"]["execution_state"] == "idle":
                break
        text = "".join(output)
        print(text, end="")
        if "KERNEL OK" not in text:
            print("FAIL: expected 'KERNEL OK' in kernel output")
            return 1
        print("HEADLESS KERNEL TEST PASSED")
        return 0
    finally:
        kc.stop_channels()
        km.shutdown_kernel(now=False)


if __name__ == "__main__":
    sys.exit(main())
