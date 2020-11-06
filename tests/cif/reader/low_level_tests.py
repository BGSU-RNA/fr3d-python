import pytest

import numpy as np
import itertools as it

from fr3d.definitions import RNAbaseheavyatoms as HEAVY

from fr3d.cif.reader import MissingColumn
from fr3d.cif.reader import MissingBlockException

from tests.cif import ReaderTest


class SimpleCIFTest(ReaderTest):
    name = '1FAT'

    def setUp(self):
        self.cif = self.__class__.cif

    def test_sets_the_pdb(self):
        val = self.cif.pdb
        ans = '1FAT'
        self.assertEqual(ans, val)

    def test_gets_table_with_leading_underscore(self):
        val = self.cif.table('_pdbx_poly_seq_scheme')
        self.assertTrue(val is not None)

    def test_gets_table_without_leading_underscore(self):
        val = self.cif.table('pdbx_poly_seq_scheme')
        self.assertTrue(val is not None)

    def test_attribute_gives_table(self):
        val = self.cif.pdbx_poly_seq_scheme
        self.assertTrue(val is not None)

    def test_fails_getting_unknown_table(self):
        self.assertRaises(MissingBlockException, self.cif.table, 'bob')

    def test_raises_key_when_getting_missing_attribute(self):
        self.assertRaises(AttributeError, lambda: self.cif.bob)

    def test_knows_if_something_is_water(self):
        self.assertTrue(self.cif.is_water('5'))

    def test_knows_if_something_is_not_water(self):
        self.assertFalse(self.cif.is_water('1'))

    def test_knows_if_an_atom_is_polymeric(self):
        self.assertTrue(self.cif.is_polymeric('1'))

    def test_knows_if_an_atom_is_not_polymeric(self):
        self.assertFalse(self.cif.is_polymeric('2'))

    def test_can_get_symmetry_operator_by_asym_id(self):
        val = [op['name'] for op in self.cif.operators('A')]
        ans = ['1_555']
        self.assertEqual(ans, val)

    def test_loads_all_symmetry_operators(self):
        self.assertEqual(2, len(self.cif._operators))

    @pytest.mark.skip()
    def test_loads_correct_symmetry_operatrs(self):
        pass

    def test_knows_if_has_block(self):
        assert self.cif.has_table('atom_site')
        assert self.cif.has_table('_atom_site')

    def test_knows_if_does_not_have_a_block(self):
        assert self.cif.has_table('bob') is False


class SimpleTableTest(ReaderTest):
    name = '1FAT'

    def setUp(self):
        self.data = self.__class__.cif.table('pdbx_poly_seq_scheme')

    def test_gets_all_columns(self):
        val = self.data.columns
        ans = ['asym_id', 'entity_id', 'seq_id', 'mon_id', 'ndb_seq_num',
               'pdb_seq_num', 'auth_seq_num', 'pdb_mon_id', 'auth_mon_id',
               'pdb_strand_id', 'pdb_ins_code', 'hetero']
        self.assertEquals(val, ans)

    def test_len_is_row_count(self):
        val = len(self.data)
        ans = 1008
        self.assertEqual(val, ans)

    def test_gets_size(self):
        ans = (1008, 12)
        val = self.data.size()
        self.assertEqual(val, ans)

    def test_gets_a_row(self):
        ans = {
            'asym_id': 'A',
            'entity_id': '1',
            'seq_id': '1',
            'mon_id': 'SER',
            'ndb_seq_num': '1',
            'pdb_seq_num': '1',
            'auth_seq_num': '1',
            'pdb_mon_id': 'SER',
            'auth_mon_id': 'SER',
            'pdb_strand_id': 'A',
            'pdb_ins_code': '.',
            'hetero': 'n'
        }
        val = self.data.rows[0]
        self.assertEqual(val, ans)

    def test_fails_getting_too_large_row(self):
        self.assertRaises(IndexError, lambda: self.data.rows[9000])

    def test_iterates_over_all_rows(self):
        ans = 1008
        val = len(self.data.rows)
        self.assertEqual(val, ans)

    def test_gets_a_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data.column('asym_id'))))
        self.assertEqual(val, ans)

    def test_fails_getting_missing_column(self):
        self.assertRaises(MissingColumn, self.data.column, 'bob')

    def test_get_item_can_give_row(self):
        ans = {
            'asym_id': 'A',
            'entity_id': '1',
            'seq_id': '1',
            'mon_id': 'SER',
            'ndb_seq_num': '1',
            'pdb_seq_num': '1',
            'auth_seq_num': '1',
            'pdb_mon_id': 'SER',
            'auth_mon_id': 'SER',
            'pdb_strand_id': 'A',
            'pdb_ins_code': '.',
            'hetero': 'n'
        }
        val = self.data[0]
        self.assertEqual(val, ans)

    def test_get_item_can_give_subtable(self):
        ans = ['SER', 'ASN']
        val = self.data[0:2].mon_id
        self.assertEqual(val, ans)

    def test_get_item_can_give_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data['asym_id'])))
        self.assertEqual(val, ans)

    def test_dot_gives_column(self):
        ans = ['A', 'B', 'C', 'D']
        val = sorted(list(set(self.data.asym_id)))
        self.assertEqual(val, ans)

    def test_get_item_on_missing_string_gives_key(self):
        self.assertRaises(KeyError, lambda: self.data['bob'])

    def test_dot_on_missing_column_gives_attribute(self):
        self.assertRaises(AttributeError, lambda: self.data.bob)

    def test_get_item_on_too_big_int_gives_index(self):
        self.assertRaises(IndexError, lambda: self.data[90000])


class StructureWithSymmetry(ReaderTest):
    name = '1WMQ'

    def test_loads_the_vector(self):
        val = [op['vector'] for op in self.cif.operators('A')]
        ans = [np.array([0.0, 0.0, 0.0])] * 3
        np.testing.assert_array_almost_equal(ans, val)

    def test_loads_the_matrix(self):
        val = [op['matrix'] for op in self.cif.operators('A')]
        ans = [np.array([[1.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0],
                         [0.0, 0.0, 1.0]]),
               np.array([[-0.5, -0.8660254,  0],
                         [0.8660254, -0.5, 0],
                         [0, 0, 1]]),
               np.array([[-0.5,  0.8660254,  0],
                         [-0.8660254, -0.5, 0],
                         [0, 0, 1]])]
        np.testing.assert_array_almost_equal(ans, val)


class StructureWithTransformationVector(ReaderTest):
    name = '4NGG'

    def test_it_builds_the_correct_rotation_matrix(self):
        val = self.cif._operators['1']['matrix']
        ans = np.array([[1.0, 0.0, 0.0],
                        [0.0, -1.0, 0.0],
                        [0.0, 0.0, -1.0]])
        np.testing.assert_array_almost_equal(ans, val)

    def test_it_builds_the_corect_symmetry_vector(self):
        val = self.cif._operators['1']['vector']
        ans = np.array([0.0, 97.240, 0.0])
        np.testing.assert_array_almost_equal(ans, val)

    def test_it_builds_the_correct_transformation_matrix(self):
        val = self.cif._operators['1']['transform']
        ans = np.array([[1.0, 0.0, 0.0, 0.0],
                        [0.0, -1.0, 0.0, 97.240],
                        [0.0, 0.0, -1.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0]])
        np.testing.assert_array_almost_equal(ans, val)


class ProblematicReadingTest(ReaderTest):
    name = '1AQ3'

    def test_can_set_default_name(self):
        atoms = it.imap(lambda r: r.atoms(), self.structure.residues())
        atoms = it.chain.from_iterable(atoms)
        atom = next(atom for atom in atoms if atom.symmetry != 'I')
        self.assertEquals('P_1', atom.symmetry)


class MissingAssemblyReadingTest(ReaderTest):
    name = '2UUA'

    def test_will_give_all_operators_if_unknown_chain(self):
        assert [op['name'] for op in self.cif.operators('X')] == ['1_555']


class AltIdReadingTest(ReaderTest):
    name = '3R1E'

    def residue(self, unit_id):
        residues = self.structure.residues()
        return next(r for r in residues if r.unit_id() == unit_id)

    def test_it_loads_all_atoms_with_alt_ids(self):
        unit_ids = [
            '3R1E|1|A|C|5||B',
            '3R1E|1|A|C|5||A',
        ]
        for unit_id in unit_ids:
            residues = self.structure.residues()
            val = next(r for r in residues if r.unit_id() == unit_id)
            msg = "Incomplete: %s" % unit_id
            assert val.is_complete(HEAVY[val.sequence]), msg

    def test_does_not_generate_unneeded_unit_without_alt_ids(self):
        unit_ids = set(r.unit_id() for r in self.structure.residues())
        assert '3R1E|1|A|C|5' not in unit_ids

    def test_loads_all_residues(self):
        ans = [
            '3R1E|1|A|G|1',
            '3R1E|1|A|C|2',
            '3R1E|1|A|GRB|3',
            '3R1E|1|A|G|4',
            '3R1E|1|A|C|5||A',
            '3R1E|1|A|C|5||B',
            '3R1E|1|A|G|6||A',
            '3R1E|1|A|G|6||B',
            '3R1E|1|A|G|7',
            '3R1E|1|A|C|8',
        ]
        residues = self.structure.residues(chain='A')
        assert sorted([r.unit_id() for r in residues]) == sorted(ans)

    def test_it_does_not_copy_over_unneeded_atoms(self):
        residue = self.residue('3R1E|1|A|G|6||A')
        atom = list(residue.atoms(name='N1'))[0]
        np.testing.assert_almost_equal(atom.x, -11.480)
