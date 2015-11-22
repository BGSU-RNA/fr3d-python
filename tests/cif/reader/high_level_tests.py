import numpy as np
from nose import SkipTest

from tests.cif import ReaderTest


class StructureTest(ReaderTest):
    name = '1GID'

    def test_assigns_pdb_id(self):
        val = self.structure.pdb
        ans = '1GID'
        self.assertEqual(ans, val)

    def test_loads_all_components(self):
        val = len(list(self.structure.residues(polymeric=None)))
        ans = 350
        self.assertEqual(ans, val)

    def test_loads_all_rna(self):
        val = len(list(self.structure.residues(type='RNA linking')))
        ans = 316
        self.assertEqual(ans, val)

    def test_can_get_a_model(self):
        raise SkipTest()
        val = self.structure.select(model=1).unit_id()
        self.assertEqual('1GID|1', val)

    def test_can_get_a_chain(self):
        raise SkipTest()
        val = self.structure.select(model=0, chain='A').unit_id()
        self.assertEqual('1GID|0|A', val)

    def test_knows_if_structure_is_true(self):
        self.assertTrue(self.structure)


class ResidueTest(ReaderTest):
    name = '1GID'

    def setUp(self):
        super(ResidueTest, self).setUp()
        self.residues = list(self.structure.residues())

    def test_assigns_atoms_to_residues(self):
        self.residues.sort(key=lambda r: r.number)
        val = len(list(self.residues[49].atoms()))
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
        val = self.residues[0].insertion_code
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

    def test_assigns_type(self):
        val = self.residues[0].type
        ans = 'RNA linking'
        self.assertEqual(ans, val)


class StructureWithSymmetry(ReaderTest):
    name = '1WMQ'

    def setUp(self):
        super(StructureWithSymmetry, self).setUp()
        self.data = list(self.structure.residues())

    def test_can_compute_correct_standard_coordinates(self):
        val = self.structure.residue('1WMQ|1|C|U|3').centers['N1']
        ans = np.array([12.788, 17.797, 104.281])
        np.testing.assert_array_almost_equal(ans, val, decimal=3)

    def test_can_compute_correct_transformed_coordinates(self):
        val = self.structure.residue('1WMQ|1|C|U|1||||2_555').centers['N1']
        ans = np.array([-20.04338, -10.465848, 108.06])
        np.testing.assert_array_almost_equal(ans, val)


class StructureWithTransformationMatrix(ReaderTest):
    name = '4NGG'

    def test_can_compute_correct_transformed_coordinates(self):
        val = self.structure.residue('4NGG|1|B|A|1||||4_565').centers["C4'"]
        ans = np.array([-19.794, 41.314, 11.949])
        np.testing.assert_array_almost_equal(ans, val, decimal=3)

    def test_can_compute_standard_coordinates(self):
        val = self.structure.residue('4NGG|1|B|A|1').centers["C4'"]
        ans = np.array([-19.794, 55.926, -11.949])
        np.testing.assert_array_almost_equal(ans, val, decimal=3)


class ReadsAtomsTest(ReaderTest):
    name = '4NMG'

    def setUp(self):
        super(ReadsAtomsTest, self).setUp()
        self.data = list(self.structure.residue('4NMG|1|A|G|2648||A').atoms())

    def test_reads_in_atom_group(self):
        val = self.data[0].group
        ans = 'ATOM'
        self.assertEquals(ans, val)

    def test_reads_in_alt_id(self):
        val = self.data[0].alt_id
        ans = 'A'
        self.assertEquals(ans, val)

    def test_reads_in_type(self):
        val = self.data[0].type
        ans = 'P'
        self.assertEquals(ans, val)

    def test_reads_in_polymeric(self):
        val = self.data[0].polymeric
        ans = True
        self.assertEquals(ans, val)

    def test_reads_in_name(self):
        val = self.data[0].name
        ans = "P"
        self.assertEquals(ans, val)

    def test_reads_in_the_symmetry(self):
        val = self.data[0].symmetry
        ans = '1_555'
        self.assertEquals(ans, val)

    def test_generates_the_correct_unit_id(self):
        val = self.data[0].unit_id()
        ans = '4NMG|1|A|G|2648|P|A'
        self.assertEquals(ans, val)

    def test_generates_correct_componet_unit_id(self):
        val = self.data[0].component_unit_id()
        ans = '4NMG|1|A|G|2648||A'
        self.assertEquals(ans, val)


class ReadsComponentsTest(ReaderTest):
    name = '4NMG'

    def setUp(self):
        super(ReadsComponentsTest, self).setUp()
        self.data = self.structure.residue('4NMG|1|A|G|2648||A')

    def test_reads_in_alt_id(self):
        self.assertEquals('A', self.data.alt_id)


class UnsortedAtomsTest(ReaderTest):
    name = '4CS1'

    def test_loads_all_residues(self):
        residues = self.structure.residues(polymeric=None)
        val = len([res.unit_id() for res in residues])
        self.assertEquals(88, val)


class DuplicateOperatorsTest(ReaderTest):
    name = '4MCE'

    def test_loads_all_residues(self):
        residues = self.structure.residues(polymeric=None)
        self.cif.structure()
        val = [res.unit_id() for res in residues]
        self.assertEquals(291, len(val))


class DuplicateWithPointSymmetryTest(ReaderTest):
    name = '4OQ8'

    def test_loads_all_residues(self):
        residues = self.structure.residues(polymeric=None)
        val = [res.unit_id() for res in residues]
        self.assertEquals(954, len(val))
        # self.assertEquals(894, len(val))


class UniqueUnitIdsTest(ReaderTest):
    name = '5AJ3'

    def test_it_has_no_duplicate_unit_ids(self):
        val = len(set(r.unit_id() for r in self.structure.residues()))
        ans = len(list(r.unit_id() for r in self.structure.residues()))
        self.assertEquals(ans, val)


class AltIdTest(ReaderTest):
    name = '1CGM'

    def test_it_has_no_duplicate_unit_ids(self):
        val = len(set(r.unit_id() for r in self.structure.residues()))
        ans = len(list(r.unit_id() for r in self.structure.residues()))
        self.assertEquals(ans, val)

    def test_it_builds_all_units(self):
        val = list(self.structure.residues(symmetry='P_25'))
        self.assertEquals(164, len(val))

    def test_it_build_atoms_with_empty_alt_id(self):
        residues = list(self.structure.residues(symmetry='P_25'))
        val = list(residues[0].atoms())[0]
        self.assertEquals(None, val.alt_id)

    def test_it_builds_atoms_with_correct_component_unit_id(self):
        residues = list(self.structure.residues(symmetry='P_25'))
        val = list(residues[0].atoms())[0]
        self.assertEquals(residues[0].unit_id(), val.component_unit_id())
