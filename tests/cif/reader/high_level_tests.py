import pytest
import numpy as np

from tests.cif import ReaderTest


class SolitaryTest(ReaderTest):
    name = '4FEO'

    def test_loads_all_rna(self):
        val = len(list(self.structure.residues(type=['RNA'])))
        ans = 68
        self.assertEqual(ans, val)


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
        val = len(list(self.structure.residues(type=['RNA'])))
        ans = 316
        self.assertEqual(ans, val)

    def test_loads_all_rna_linking(self):
        val = len(list(self.structure.residues(type=['RNA linking'])))
        ans = 316
        self.assertEqual(ans, val)

    def test_can_get_a_model(self):
        val = self.structure.select(model=1).unit_id()
        self.assertEqual('1GID|1', val)

    def test_can_get_a_chain(self):
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
        ans = 24
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
        ans = '1'
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
        self.residues.sort(key=lambda r: '%s%s' % (r.chain, r.number))

        print(self.residues[0])

        val = self.residues[0].type
        ans = 'RNA linking'
        self.assertEqual(ans, val)


class ChainIndexTest(ReaderTest):
    name = '4A3J'

    def setUp(self):
        super(ChainIndexTest, self).setUp()

        # for r in list(self.structure.residues(polymeric=False)):
        #     print(r.unit_id(), r.index)

        # self.data = sorted(self.structure.residues(polymeric=False),
        #                    key=lambda c: c.index)

        self.data = sorted(self.structure.residues(chains='A'),
                           key=lambda c: c.index)
        # self.data = [d for d in self.data if d.chain == 'A']

    def test_it_assigns_numbers_to_all_indexes(self):
        # actually, solitary nucleotides will have no numeric index
        # 4A3J has 4A3J|1|T|G2P|1024 which is not polymeric

        for d in self.data:
            if d.index is not None:
                assert isinstance(d.index, int)
            else:
                assert d.polymeric is False

    def test_it_assigns_the_correct_index(self):
        # Non polymeric entries have no index so should set to none.
        for d in self.structure.residues(polymeric=False):
            assert d.index is None


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
        self.assertEqual(ans, val)

    def test_reads_in_alt_id(self):
        val = self.data[0].alt_id
        ans = 'A'
        self.assertEqual(ans, val)

    def test_reads_in_type(self):
        val = self.data[0].type
        ans = 'P'
        self.assertEqual(ans, val)

    def test_reads_in_polymeric(self):
        val = self.data[0].polymeric
        ans = True
        self.assertEqual(ans, val)

    def test_reads_in_name(self):
        val = self.data[0].name
        ans = "P"
        self.assertEqual(ans, val)

    def test_reads_in_the_symmetry(self):
        val = self.data[0].symmetry
        ans = '1_555'
        self.assertEqual(ans, val)

    def test_generates_the_correct_unit_id(self):
        val = self.data[0].unit_id()
        ans = '4NMG|1|A|G|2648|P|A'
        self.assertEqual(ans, val)

    def test_generates_correct_componet_unit_id(self):
        val = self.data[0].component_unit_id()
        ans = '4NMG|1|A|G|2648||A'
        self.assertEqual(ans, val)


class ReadsComponentsTest(ReaderTest):
    name = '4NMG'

    def setUp(self):
        super(ReadsComponentsTest, self).setUp()
        self.data = self.structure.residue('4NMG|1|A|G|2648||A')

    def test_reads_in_alt_id(self):
        self.assertEqual('A', self.data.alt_id)


class UnsortedAtomsTest(ReaderTest):
    name = '4CS1'

    def test_loads_all_residues(self):
        residues = self.structure.residues(polymeric=None)
        val = len([res.unit_id() for res in residues])
        self.assertEqual(88, val)


class DuplicateOperatorsTest(ReaderTest):
    name = '4MCE'
    # original version from PDB has two non-trivial symmetry operators

    def test_loads_all_residues(self):
        residues = self.structure.residues(polymeric=None)
        # self.cif.structure()
        val = [res.unit_id() for res in residues]
        print(val)
        self.assertEqual(291, len(val))


class DuplicateWithPointSymmetryTest(ReaderTest):
    name = '4OQ8'

    def test_loads_all_chain_A(self):
        residues = self.structure.residues(chain='A', symmetry='1_555')
        val = [res.unit_id() for res in residues]
        assert len(val) == 428

    def test_loads_all_chain_B(self):
        residues = self.structure.residues(chain='B', symmetry='1_555')
        val = [res.unit_id() for res in residues]
        assert len(val) == 16

    def test_loads_all_chain_C(self):
        residues = self.structure.residues(chain='C', symmetry='1_555')
        val = [res.unit_id() for res in residues]
        assert len(val) == 16

    def test_loads_all_chain_D(self):
        residues = self.structure.residues(chain='D', symmetry='1_555')
        val = [res.unit_id() for res in residues]
        assert len(val) == 3

    def test_loads_all_residues(self):
        residues = self.structure.residues(polymeric=None)
        val = [res.unit_id() for res in residues]
        self.assertEqual(27780, len(val))


class UniqueUnitIdsTest(ReaderTest):
    name = '5AJ3'

    def test_it_has_no_duplicate_unit_ids(self):
        val = len(set(r.unit_id() for r in self.structure.residues()))
        ans = len(list(r.unit_id() for r in self.structure.residues()))
        self.assertEqual(ans, val)


class AltIdTest(ReaderTest):
    name = '1CGM'

    def test_it_has_no_duplicate_unit_ids(self):
        val = len(set(r.unit_id() for r in self.structure.residues()))
        ans = len(list(r.unit_id() for r in self.structure.residues()))

        # seen = set()
        # for r in self.structure.residues():
        #     u = r.unit_id()
        #     print('test_it_has_no_duplicate_unit_ids',u,r.index)
        #     if u in seen:
        #         print('Already saw %s' % u)
        #     else:
        #         seen.add(u)

        self.assertEqual(ans, val)

    def test_it_builds_all_units(self):
        val = list(self.structure.residues(symmetry='ASM_24'))
        self.assertEqual(169, len(val))

    def test_it_build_atoms_with_empty_alt_id(self):
        residues = list(self.structure.residues(symmetry='ASM_24'))
        val = list(residues[0].atoms())[0]
        print(val)
        self.assertEqual(None, val.alt_id)

    def test_it_builds_atoms_with_correct_component_unit_id(self):
        residues = list(self.structure.residues(symmetry='ASM_24'))
        val = list(residues[0].atoms())[0]
        self.assertEqual(residues[0].unit_id(), val.component_unit_id())


class DistanceTreeTest(ReaderTest):
    name = '1GID'

    def setUp(self):
        super(DistanceTreeTest, self).setUp()
        self.tree = self.structure.distances()

    def test_it_can_create_a_residue_tree(self):
        counts = self.tree.count_neighbors(self.tree, 3)
        assert counts == 350


class AtomDistanceTreeTest(ReaderTest):
    name = '1GID'

    def setUp(self):
        super(AtomDistanceTreeTest, self).setUp()
        self.tree = self.structure.atom_distances()

    def test_it_can_create_an_atom_tree(self):
        counts = self.tree.count_neighbors(self.tree, 3)
        assert counts >= 316


class NonpolymerChains(ReaderTest):
    name = '2UUA'

    def test_it_loads_paryomcin_residues(self):
        chains = set(r.sequence for r in self.structure._residues)
        assert 'PAR' in chains

    def test_it_loads_chain_z(self):
        residues = self.structure.residues(polymeric=None)
        chains = set(r.chain for r in residues)
        assert 'Z' in chains
