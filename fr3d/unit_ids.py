"""This is a package for generating unit ids for units from data. This checks
if the id is valid, which means it contains all the required fields.
"""

SEPERATOR = '|'

FIELDS = ['pdb', 'model', 'chain', 'component_id', 'component_number',
          'atom_name', 'alt_id', 'insertion_code', 'symmetry']

DEFAULTS = {
    'pdb': None,
    'model': None,
    'chain': None,
    'component_id': None,
    'component_number': None,
    'atom_name': '',
    'alt_id': '',
    'insertion_code': '',
    'symmetry': '1_555'
}


class InvalidUnitId(Exception):
    """This is generated whenever we attempt to encode an invalid unit id. This
    means we raise if this exception if we are missing required fields such as
    PDB, model and chain.
    """
    pass


def encode(data, full=False):
    """Generate a unit ID for the given data. All possible fields in the
    given dictonary will be used to generate the id.

    :data: Dictonary of fields to use in unit id.
    :returns: A string of the new unit id.
    """

    ordered = []

    defaults = DEFAULTS

    for field in FIELDS:
        default = defaults[field]
        if default is None and field not in data:
            raise InvalidUnitId("Missing required field: " + field)

        value = str(data.get(field, default))
        ordered.append(value)

    if not full:
        possible = ['symmetry', 'insertion_code', 'alt_id', 'atom_name']
        while possible and ordered[-1] == DEFAULTS[possible[0]]:
            ordered.pop()
            possible.pop(0)

    return SEPERATOR.join(ordered)


def decode(unit_id):
    """Turn a unit id into a dictonary of it's components. This will infer any
    fields which are not shown, such as symmetry and atom name.

    :unit_id: A unit id string.
    :returns: A dictonary of field: value pairs for the given unit id.
    """

    parts = unit_id.split(SEPERATOR)
    total = dict(DEFAULTS)
    fields = FIELDS[slice(0, len(parts))]

    total.update(dict(zip(fields, parts)))

    if total['model'] is not None:
        total['model'] = int(total['model'])
    if total['component_number'] is not None:
        total['component_number'] = int(total['component_number'])

    return total
