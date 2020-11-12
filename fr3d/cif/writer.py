import itertools as it

from pdbx.writer.PdbxWriter import PdbxWriter as Writer

from pdbx.reader.PdbxContainers import DataCategory
from pdbx.reader.PdbxContainers import DataContainer


class CifAtom(object):
    """
    This is a class to write cifatom files. These are partial cif files that
    contain all atoms, including those created by symmetry operations. These
    files differ from normal cif files in that they only contain the atom_site
    block and this block differs from that in the standard cif files.

    1. This file conatins all atoms, even those produced by rotations.
    2. All label_* and auth_* fields will be the same and will contain the
        values used to compute the unit id of the atom.
    3. The entity id is '?'.
    4. All *_esd, charge, *_esi, occupancy and such entries are '?'.
    5. An additional final field contains the component unit id for each atom.
    """

    def __init__(self, handle, unit_ids=True, protect_lists_of_lists=False):
        self.writer = Writer(handle)
        self.unit_ids = unit_ids
        self.protect_lists_of_lists = protect_lists_of_lists


    def atom_container(self, structure, protect_lists_of_lists):
        atoms = DataCategory('atom_site')
        fields = ['group_PDB', 'id', 'type_symbol', 'label_atom_id',
                  'label_alt_id', 'label_comp_id', 'label_asym_id',
                  'label_entity_id', 'label_seq_id', 'pdbx_PDB_ins_code',
                  'Cartn_x', 'Cartn_y', 'Cartn_z', 'occupancy',
                  'B_iso_or_equiv', 'Cartn_x_esd', 'Cartn_y_esd',
                  'Cartn_z_esd', 'occupancy_ies', 'B_iso_or_equiv_esd',
                  'pdbx_formal_charge', 'auth_seq_id', 'auth_comp_id',
                  'auth_asym_id', 'auth_atom_id', 'pdbx_PDB_model_num']

        if self.unit_ids:
            fields.append('unit_id')

        for field in fields:
            atoms.appendAttribute(field)

        def key(atom):
            return (atom.symmetry, atom.model, atom.chain,
                    atom.component_number, atom.insertion_code)

        all_atoms = it.imap(lambda r: r.atoms(),
                            structure.residues(polymeric=None))
        all_atoms = it.chain.from_iterable(all_atoms)
        for index, atom in enumerate(sorted(all_atoms, key=key)):
            alt_id = getattr(atom, 'alt_id', '.')
            data = [atom.group, index, atom.type, atom.name,
                    alt_id, atom.component_id, atom.chain,
                    '?', atom.component_number, atom.insertion_code,
                    atom.x, atom.y, atom.z, '?',
                    '?', '?', '?',
                    '?', '?', '?',
                    '.', atom.component_number, atom.component_id,
                    atom.chain, atom.name, atom.model]

            # atoms added by infer_hydrogens don't have information to make
            # the unit_id and should not be written out anyway
            # So if the next line fails, don't write the atom out.
            try:
                unit_id = atom.component_unit_id()

                if self.unit_ids:
                    data.append(atom.component_unit_id())

                atoms.append(data)
            except:
                continue

        if protect_lists_of_lists is True:
            # Kludge fix for single atom residues.
            # Handles cases where single atom residues are not handled
            #     correctly as lists of lists, but as simple lists instead,
            #     which has downstream implications in units.coordinates.
            #
            # Here, we force a fix, using the line-skipping logic from
            #     units.coordinates to keep the kludge line out of the output
            #     data.
            dummy = [ 'loop_foo', 'foo', 'foo', 'foo', 'foo', 'foo', 'foo',
                        '?', 'foo', 'foo', 'foo', 'foo', 'foo', '?',
                        '?', '?', '?', '?', '?', '?',
                        '.', 'foo', 'foo', 'foo', 'foo', 'foo']

            atoms.append(dummy)

        return atoms

    def __call__(self, structure):
        atoms = self.atom_container(structure, self.protect_lists_of_lists)
        container = DataContainer(structure.pdb)
        container.append(atoms)
        self.writer.writeContainer(container)
