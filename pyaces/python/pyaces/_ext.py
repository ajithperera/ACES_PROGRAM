"""Locate and import the aces2py extension built by pyaces/build.sh."""

import importlib
import os
import sys

# python/pyaces/_ext.py -> pyaces/Main/build, where build.sh leaves the .so.
_DEFAULT_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Main', 'build'))

_module = None


def ext_dir():
    """Directory searched for aces2py*.so: $PYACES_EXT_DIR, else pyaces/Main/build."""
    return os.environ.get('PYACES_EXT_DIR', _DEFAULT_DIR)


def ext():
    """Return the raw aces2py module (Init, Runcc, Dgetrecpy, ...), importing it once."""
    global _module
    if _module is None:
        d = ext_dir()
        if d not in sys.path:
            sys.path.insert(0, d)
        try:
            _module = importlib.import_module('aces2py')
        except ImportError as e:
            raise ImportError(
                f'cannot import aces2py from {d}: {e}. Build it with '
                f'./build.sh pyaces and load the same toolchain modules '
                f'(Intel runtime and MKL) before importing.') from e
    return _module


def recname(name):
    """JOBARC record names are exactly 8 characters, blank padded."""
    if len(name) > 8:
        raise ValueError(f'JOBARC record name longer than 8 characters: {name!r}')
    return name.ljust(8)
