import numpy as np

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


class SequenceMappingTest(ReaderTest):
    name = '1GID'

    def setUp(self):
        super(SequenceMappingTest, self).setUp()
        self.data = self.cif.experimental_sequence_mapping('A')

    def test_can_compute_mapping(self):
        val = self.data[-1]
        ans = ('C', '1GID|Sequence|A|C|158', '1GID|1|A|C|260')
        self.assertEqual(ans, val)
