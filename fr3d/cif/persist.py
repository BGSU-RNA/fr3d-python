"""This package deals with writing and resotring cif file objects.
"""

import cPickle as pickle

from fr3d.cif.reader import Cif


def serialize(handle, cif):
    """Serailze a cif file. This will write out the data that a cif file object
    contains into a pickle formatted file. The object can then be rebuilt using
    deserialize.

    :handle: The filehandle to write to.
    :cif: The Cif object to persist.
    """
    data = cif.data
    catalog = data._ContainerBase__objCatalog
    for name, category in catalog.items():
        category._DataCategory__lfh = None
    pickle.dump(data, handle)


def deserialize(handle):
    """Load a serialized cif file.

    :handle: The filehandle to read from.
    :returns: A new Cif object from the persisted data.
    """
    data = pickle.load(handle)
    return Cif(data=data)
