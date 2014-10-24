import os
from unittest import TestCase

from fr3d.cif.reader import CIF
from fr3d.cif.reader import MissingColumn
from fr3d.cif.reader import MissingBlockException


class ReaderTest(TestCase):
    name = None

    @classmethod
    def setUpClass(cls):
        with open(os.path.join('files', cls.name + '.cif'), 'rb') as raw:
            cls.cif = CIF(raw)
            cls.structure = cls.cif.structure()

    def setUp(self):
        self.structure = self.__class__.structure


class SimpleCIFTest(ReaderTest):
    name = '1FAT'

    def setUp(self):
        self.cif = self.__class__.cif

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



class ReaderStructureTest(ReaderTest):
    name = '1GID'

    def test_assigns_pdb_id(self):
        val = self.structure.pdb
        ans = '1GID'
        self.assertEqual(ans, val)

    def test_loads_all_models(self):
        val = len(self.structure.models)
        ans = 1
        self.assertEqual(ans, val)

    def test_loads_all_chains(self):
        print(list(self.structure.chains()))
        val = len(list(self.structure.chains()))
        ans = 2
        self.assertEqual(ans, val)

    def test_loads_all_components(self):
        val = len(list(self.structure.residues()))
        ans = 350
        self.assertEqual(ans, val)

    # def test_loads_all_atoms(self):
    #     self.fail()

    def test_can_get_a_model(self):
        val = self.structure.model(1).model
        ans = 1
        self.assertEqual(ans, val)

    def test_fails_getting_invalid_model(self):
        self.assertRaises(IndexError, self.structure.model, 2)

    def test_can_get_a_chain(self):
        val = self.structure.chain(0, 'A').chain
        ans = 'A'
        self.assertEqual(ans, val)

    def test_fails_getting_invalid_chain(self):
        self.assertTrue(self.structure.chain(1, 'C') is None)


class ReaderResidueTest(ReaderTest):
    name = '1GID'

    def setUp(self):
        super(ReaderResidueTest, self).setUp()
        self.residues = list(self.structure.residues())

    def test_assigns_atoms_to_residues(self):
        self.residues.sort(key=lambda r: r.number)
        val = len(self.residues[49].atoms())
        ans = 20
        self.assertEqual(ans, val)

    def test_assigns_numbers_correctly(self):
        self.residues.sort(key=lambda r: r.number)
        val = self.residues[0].number
        ans = 103
        self.assertEqual(ans, val)

    def test_assigns_pdb(self):
        val = self.residues[0].pdb
        ans = '1GID'
        self.assertEqual(ans, val)

    def test_assigns_model(self):
        val = self.residues[0].model
        ans = 1
        self.assertEqual(ans, val)

    def test_assigns_chain(self):
        self.residues.sort(key=lambda r: '%s%s' % (r.chain, r.number))
        val = self.residues[0].chain
        ans = 'A'
        self.assertEqual(ans, val)

    def test_assigns_symmetry(self):
        val = self.residues[0].symmetry
        ans = '1_555'
        self.assertEqual(ans, val)

    def test_assigns_ins_code(self):
        val = self.residues[0].ins_code
        ans = None
        self.assertEqual(ans, val)

    def test_assigns_sequence(self):
        self.residues.sort(key=lambda r: r.number)
        val = self.residues[0].sequence
        ans = 'G'
        self.assertEqual(ans, val)

    def test_can_generate_unit_id(self):
        self.residues.sort(key=lambda r: '%s%s' % (r.chain, r.number))
        val = self.residues[0].unit_id()
        ans = '1GID|1|A|G|103'
        self.assertEqual(ans, val)

    #def test_orders_using_poly_seq(self):
        #val = [(res.number, res.chain) for res in self.residues]
        #ans = list(chain(izip(repeat(103, 260), repeat('A')),
                         #izip(repeat(103, 260), repeat('B'))))
        #self.assertEquals(val, ans)


# class LoadingPolymersTest(ReaderTest):
#     name = '2UUA'

#     def test_detects_all_polymers(self):
#         val = len(self.structure.polymers())
#         ans = 2
#         self.assertEqual(ans, val)

#     def test_detects_breaks_correctly(self):
#         val = [poly.endpoints() for poly in self.structure.polymers()]
#         ans = None
#         self.assertEqual(ans, val)


# class ReaderAtomTest(ReaderTest):

#     def setUp(self):
#         self.atoms = []
#         for residue in READER.structures()[0].residues():
#             self.atoms.extend(residue.atoms())

#     def test_loads_all_atoms(self):
#         val = len(self.atoms)
#         ans = 6824
#         self.assertEqual(ans, val)


# class ReaderSymmetryTest(ReaderTest):

#     def setUp(self):
#         self.atoms = []
#         self.residues = self.__class__.structure.residues()
#         for residue in self.residues:
#             self.atoms.extend(residue.atoms())

#     def test_loads_all_atoms(self):
#         self.fail()

#     def test_loads_all_residues(self):
#         val = len(self.residues)
#         ans = 77287
#         self.assertEqual(val, ans)
