import builtins
import itertools

from fr3d.data import Component


# Compatibility shim for legacy tests that still call infer_hydrogens().
if not hasattr(Component, "infer_hydrogens"):
    Component.infer_hydrogens = Component.infer_NA_hydrogens


# Compatibility shim for legacy Python 2-era code paths used by tests.
if not hasattr(builtins, "xrange"):
	builtins.xrange = range

if not hasattr(itertools, "imap"):
	itertools.imap = map
