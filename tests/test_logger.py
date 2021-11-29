import sys
import subprocess

CODE = r"""\
import logging
from assertionlib import assertion

assertion.len_eq(logging.root.handlers, 0)
import CAT
assertion.len_eq(logging.root.handlers, 0)
"""


def test_basic_config() -> None:
    p = subprocess.run([sys.executable, '-c', CODE], encoding="utf8", stderr=subprocess.PIPE)
    if p.returncode:
        raise AssertionError(f"Non-zero return code: {p.returncode!r}\n\n{p.stderr}")
