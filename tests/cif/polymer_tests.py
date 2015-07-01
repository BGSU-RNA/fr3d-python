# class BasicChainPolymersTest(ReaderTest):
#     name = '1FAT'

#     def setUp(self):
#         self.chain = self.__class__.structure.chain(1, 'A')

#     def test_can_get_polymers(self):
#         val = [len(p) for p in self.chain.polymers()]
#         ans = [36, 232 - 36]
#         self.assertEqual(ans, val)

#     def test_it_has_no_breaks(self):
#         self.assertFalse(self.chain.polymer(0).has_breaks())
#         self.assertFalse(self.chain.polymer(1).has_breaks())


# class PolymersTest(ReaderTest):
#     name = '1FAT'

#     def setUp(self):
#         self.cif = self.__class__.cif
#         self.data = self.cif.chain('1_555', '1', 'D').polymers()

#     def test_gets_all_polymers_in_all_chains(self):
#         val = sorted(set([poly['chain'] for poly in self.cif.polymers()]))
#         ans = ['A', 'B', 'C', 'D']
#         self.assertEqual(val, ans)

#     def test_gets_requested_chain_polymer(self):
#         val = [poly['chain'] for poly in self.data]
#         ans = ['D', 'D']
#         self.assertEqual(val, ans)

#     def test_finds_polymer_sequence(self):
#         ans = [
#             ['SER', 'ASN', 'ASP', 'ILE', 'TYR', 'PHE', 'ASN', 'PHE', 'GLN',
#              'ARG', 'PHE', 'ASN', 'GLU', 'THR', 'ASN', 'LEU', 'ILE', 'LEU',
#              'GLN', 'ARG', 'ASP', 'ALA', 'SER', 'VAL', 'SER', 'SER', 'SER',
#              'GLY', 'GLN', 'LEU', 'ARG', 'LEU', 'THR', 'ASN', 'LEU'],
#             ['ASN', 'GLY', 'GLU', 'PRO', 'ARG', 'VAL', 'GLY', 'SER', 'LEU',
#              'GLY', 'ARG', 'ALA', 'PHE', 'TYR', 'SER', 'ALA', 'PRO', 'ILE',
#              'GLN', 'ILE', 'TRP', 'ASP', 'ASN', 'THR', 'THR', 'GLY', 'THR',
#              'VAL', 'ALA', 'SER', 'PHE', 'ALA', 'THR', 'SER', 'PHE', 'THR',
#              'PHE', 'ASN', 'ILE', 'GLN', 'VAL', 'PRO', 'ASN', 'ASN', 'ALA',
#              'GLY', 'PRO', 'ALA', 'ASP', 'GLY', 'LEU', 'ALA', 'PHE', 'ALA',
#              'LEU', 'VAL', 'PRO', 'VAL', 'GLY', 'SER', 'GLN', 'PRO', 'LYS',
#              'ASP', 'LYS', 'GLY', 'GLY', 'PHE', 'LEU', 'GLY', 'LEU', 'PHE',
#              'ASP', 'GLY', 'SER', 'ASN', 'SER', 'ASN', 'PHE', 'HIS', 'THR',
#              'VAL', 'ALA', 'VAL', 'GLU', 'PHE', 'ASP', 'THR', 'LEU', 'TYR',
#              'ASN', 'LYS', 'ASP', 'TRP', 'ASP', 'PRO', 'THR', 'GLU', 'ARG',
#              'HIS', 'ILE', 'GLY', 'ILE', 'ASP', 'VAL', 'ASN', 'SER', 'ILE',
#              'ARG', 'SER', 'ILE', 'LYS', 'THR', 'THR', 'ARG', 'TRP', 'ASP',
#              'PHE', 'VAL', 'ASN', 'GLY', 'GLU', 'ASN', 'ALA', 'GLU', 'VAL',
#              'LEU', 'ILE', 'THR', 'TYR', 'ASP', 'SER', 'SER', 'THR', 'ASN',
#              'LEU', 'LEU', 'VAL', 'ALA', 'SER', 'LEU', 'VAL', 'TYR', 'PRO',
#              'SER', 'GLN', 'LYS', 'THR', 'SER', 'PHE', 'ILE', 'VAL', 'SER',
#              'ASP', 'THR', 'VAL', 'ASP', 'LEU', 'LYS', 'SER', 'VAL', 'LEU',
#              'PRO', 'GLU', 'TRP', 'VAL', 'SER', 'VAL', 'GLY', 'PHE', 'SER',
#              'ALA', 'THR', 'THR', 'GLY', 'ILE', 'ASN', 'LYS', 'GLY', 'ASN',
#              'VAL', 'GLU', 'THR', 'ASN', 'ASP', 'VAL', 'LEU', 'SER', 'TRP',
#              'SER', 'PHE', 'ALA', 'SER', 'LYS', 'LEU', 'SER']
#         ]
#         val = [poly.sequence for poly in self.data]
#         self.assertEqual(val, ans)


# class BuildingStructuresWithBreaks(ReaderTest):
#     name = '2UUA'

#     def test_detects_all_polymers(self):
#         val = len(list(self.structure.chain(1, 'A').polymers()))
#         ans = 2
#         self.assertEqual(ans, val)

#     def test_detects_breaks_correctly(self):
#         polymers = self.structure.polymers()
#         val = [(p.first().unit_id(), p.last.unit_id()) for p in polymers]
#         ans = [('2UUA|1|A|U|5', '2UUA|1|A|A|1531'),
#                ('2UUA|1|A|A|1532', '2UUA|1|A|U|1544')]
#         self.assertEqual(ans, val)
