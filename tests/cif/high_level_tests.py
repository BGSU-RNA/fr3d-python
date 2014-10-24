from tests.cif import ReaderTest


class StructureTest(ReaderTest):
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


class ResidueTest(ReaderTest):
    name = '1GID'

    def setUp(self):
        super(ResidueTest, self).setUp()
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


class ChainTest(ReaderTest):
    name = '1FAT'

    def setUp(self):
        self.data = self.__class__.structure.chain(1, 'D')

    def test_it_knows_the_chain_id(self):
        val = self.data.chain
        ans = 'D'
        self.assertEqual(ans, val)

    def test_it_only_has_rows_from_correct_chain(self):
        val = list(set([a.chain for a in self.data.atoms()]))
        ans = ['D']
        self.assertEqual(ans, val)

    def test_it_can_find_unit_id_for_index(self):
        val = self.data[1].unit_id()
        ans = '1FAT|1|D|ASN|2'
        self.assertEqual(val, ans)

    def test_it_can_get_the_first_residue(self):
        val = self.data.first().unit_id()
        ans = '1FAT|1|D|SER|1'
        self.assertEqual(val, ans)

    def test_it_can_get_the_last_residue(self):
        val = self.data.last().unit_id()
        ans = '1FAT|1|D|SER|233'
        self.assertEqual(val, ans)


# class CifBuildingStructuresWithBreaks(ReaderTest):
#     def test_detects_all_polymers(self):
#         val = len(list(self.structure.chain(1, 'A').polymers()))
#         ans = 2
#         self.assertEqual(ans, val)

#     def test_detects_breaks_correctly(self):
#         val = [poly.endpoints() for poly in self.structure.polymers()]
#         ans = None
#         self.assertEqual(ans, val)
