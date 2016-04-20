from fr3d.unit_ids import decode

from tests.cif import ReaderTest


class ExperimentalSequenceMappingTest(ReaderTest):
    chains = None

    def setUp(self):
        super(ExperimentalSequenceMappingTest, self).setUp()
        if self.chains:
            chains = self.chains
            self.data = list(self.cif.experimental_sequence_mapping(chains))
            self.mapping = self.__mapping(self.data)

    def __mapping(self, data):
        mapping = {}
        for entry in data:
            mapping[entry['unit_id']] = entry
            mapping[entry['seq_id']] = entry
        return mapping


class BasicTest(ExperimentalSequenceMappingTest):
    name = '1GID'
    chains = 'A'

    def test_it_creates_all_mappings(self):
        assert len(self.data) == 158

    def test_can_compute_mapping(self):
        assert self.mapping['1GID|1|A|C|260'] == {
            'index': 157,
            'unit_id': '1GID|1|A|C|260',
            'seq_id': '1GID|Sequence|A|C|260',
            'seq_unit': 'C',
            'number': 260,
        }

    def test_finds_no_duplicate_unit_ids(self):
        ind = [d['unit_id'] for d in self.data]
        assert len(ind) == len(set(ind))
        assert len(ind) == 158

    def test_finds_no_duplicate_indexes(self):
        ind = [d['index'] for d in self.data]
        assert len(ind) == len(set(ind))
        assert len(ind) == 158

    def test_finds_no_duplicate_seq_ids(self):
        sids = [d['seq_id'] for d in self.data]
        assert len(sids) == len(set(sids))
        assert len(sids) == 158


class MultipleChainTest(ExperimentalSequenceMappingTest):
    name = '1GID'
    chains = ('A', 'B')

    def test_can_get_for_both_chains(self):
        ans = set('AB')
        assert set(decode(d['unit_id'])['chain'] for d in self.data) == ans
        assert set(d['seq_id'].split('|')[2] for d in self.data) == ans

    def test_finds_no_duplicate_unit_ids(self):
        uids = [d['unit_id'] for d in self.data]
        assert len(uids) == len(set(uids))

    def test_finds_no_duplicate_seq_ids(self):
        sids = [d['seq_id'] for d in self.data]
        assert len(sids) == len(set(sids))

    def test_can_compute_all_mappings(self):
        # Two chains of length 158
        assert len(self.data) == 158 * 2

    def test_resets_the_index(self):
        chain_of = lambda d: decode(d['unit_id'])['chain']
        chain_a = [d['index'] for d in self.data if chain_of(d) == 'A']
        chain_b = [d['index'] for d in self.data if chain_of(d) == 'B']
        assert chain_a == range(0, 158)
        assert chain_b == range(0, 158)


class LargeSequenceMappingTest(ExperimentalSequenceMappingTest):
    name = '1S72'
    chains = '0'

    def test_can_compute_full_mapping(self):
        self.assertEqual(2922, len(self.data))

    def test_it_maps_unobserved_residus(self):
        uids = [d['unit_id'] for d in self.data]
        assert None in uids


class MutipleEntriesInExpSeqTest(ExperimentalSequenceMappingTest):
    name = '1I9K'
    chains = 'A'

    def test_can_map_exp_seq(self):
        assert self.data == [
            {'number': 1, 'index': 0, 'seq_unit': 'U', 'unit_id': '1I9K|1|A|U|1', 'seq_id': '1I9K|Sequence|A|U|1'},
            {'number': 2, 'index': 1, 'seq_unit': 'C', 'unit_id': '1I9K|1|A|C|2', 'seq_id': '1I9K|Sequence|A|C|2'},
            {'number': 3, 'index': 2, 'seq_unit': 'C', 'unit_id': '1I9K|1|A|C|3', 'seq_id': '1I9K|Sequence|A|C|3'},
            {'number': 4, 'index': 3, 'seq_unit': 'C', 'unit_id': '1I9K|1|A|C|4', 'seq_id': '1I9K|Sequence|A|C|4'},
            {'number': 5, 'index': 4, 'seq_unit': 'C', 'unit_id': '1I9K|1|A|C|5', 'seq_id': '1I9K|Sequence|A|C|5'},
            {'number': 6, 'index': 5, 'seq_unit': 'C', 'unit_id': '1I9K|1|A|C|6', 'seq_id': '1I9K|Sequence|A|C|6'},
        ]


class WithMissingTest(ExperimentalSequenceMappingTest):
    name = '1IBK'
    chains = 'A'

    def test_can_create_correct_mappings(self):
        uids = [d['unit_id'] for d in self.data]
        known = [d['unit_id'] for d in self.data if d['unit_id'] is not None]
        missing = [d['unit_id'] for d in self.data if d['unit_id'] is None]
        assert None in uids
        assert len(uids) == 1522
        assert len(known) == 1506
        assert len(missing) == 1522 - 1506


class WithNoIdentityOperator(ExperimentalSequenceMappingTest):
    name = '4OQ8'
    chains = 'B'

    def test_it_maps_all_residues(self):
        assert len(self.data) == 20

    def test_maps_both_symmetry_operators(self):
        ops = set(decode(d['unit_id'])['symmetry'] for d in self.data)
        assert ops == set(['P_P', 'P_1'])

    def test_it_creates_correct_mappings(self):
        assert self.data == [
            {'number': 161, 'index': 0, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|161||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|161'},
            {'number': 161, 'index': 0, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|161||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|161'},
            {'number': 162, 'index': 1, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|162||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|162'},
            {'number': 162, 'index': 1, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|162||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|162'},
            {'number': 163, 'index': 2, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|163||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|163'},
            {'number': 163, 'index': 2, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|163||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|163'},
            {'number': 164, 'index': 3, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|164||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|164'},
            {'number': 164, 'index': 3, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|164||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|164'},
            {'number': 165, 'index': 4, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|165||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|165'},
            {'number': 165, 'index': 4, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|165||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|165'},
            {'number': 166, 'index': 5, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|166||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|166'},
            {'number': 166, 'index': 5, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|166||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|166'},
            {'number': 167, 'index': 6, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|167||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|167'},
            {'number': 167, 'index': 6, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|167||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|167'},
            {'number': 168, 'index': 7, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|168||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|168'},
            {'number': 168, 'index': 7, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|168||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|168'},
            {'number': 169, 'index': 8, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|169||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|169'},
            {'number': 169, 'index': 8, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|169||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|169'},
            {'number': 170, 'index': 9, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|170||A||P_1', 'seq_id': '4OQ8|Sequence|B|A|170'},
            {'number': 170, 'index': 9, 'seq_unit': 'A', 'unit_id': '4OQ8|1|B|A|170||A||P_P', 'seq_id': '4OQ8|Sequence|B|A|170'},
        ]


class WithNoIdentityOperators2(ExperimentalSequenceMappingTest):
    name = '4OQ9'
    chains = '1'

    def test_maps_both_symmetry_operators(self):
        ops = set(decode(d['unit_id'])['symmetry'] for d in self.data)
        assert ops == set(['P_P', 'P_1'])


class WithNonStandardModel(ExperimentalSequenceMappingTest):
    name = '4R3I'
    chains = 'B'

    def test_can_generate_mapping_to_model_0(self):
        assert self.data == [
            {'index': 0, 'seq_id': '4R3I|Sequence|B|G|1', 'seq_unit': 'G', 'number': 1, 'unit_id': '4R3I|0|B|G|1'},
            {'index': 1, 'seq_id': '4R3I|Sequence|B|G|2', 'seq_unit': 'G', 'number': 2, 'unit_id': '4R3I|0|B|G|2'},
            {'index': 2, 'seq_id': '4R3I|Sequence|B|6MZ|3', 'seq_unit': '6MZ', 'number': 3, 'unit_id': '4R3I|0|B|6MZ|3'},
            {'index': 3, 'seq_id': '4R3I|Sequence|B|C|4', 'seq_unit': 'C', 'number': 4, 'unit_id': '4R3I|0|B|C|4'},
            {'index': 4, 'seq_id': '4R3I|Sequence|B|U|5', 'seq_unit': 'U', 'number': 5, 'unit_id': '4R3I|0|B|U|5'},
        ]


class WithDuplicateEntries(ExperimentalSequenceMappingTest):
    name = '4X4N'
    chains = 'G'

    def test_it_creates_correct_number_of_mappings(self):
        # 14 (unosbered) + 2 * 18 (observed with 2 alt ids each)
        self.assertEquals(50, len(self.data))

    def test_it_can_map_to_none(self):
        assert None in self.mapping
        assert self.mapping['4X4N|Sequence|G|C|20'] == {
            'seq_id': '4X4N|Sequence|G|C|20',
            'unit_id': None,
            'index': 19,
            'number': 20,
            'seq_unit': 'C'
        }

    def test_it_takes_the_first_entry(self):
        assert self.mapping['4X4N|1|G|A|29||A'] == {
            'index': 28,
            'seq_id': '4X4N|Sequence|G|A|29',
            'unit_id': '4X4N|1|G|A|29||A',
            'number': 29,
            'seq_unit': 'A',
        }
        assert self.mapping['4X4N|1|G|G|29||B'] == {
            'seq_id': '4X4N|Sequence|G|A|29',
            'unit_id': '4X4N|1|G|G|29||B',
            'seq_unit': 'A',
            'number': 29,
            'index': 28,
        }
        assert self.mapping['4X4N|1|G|U|30||A'] == {
            'seq_id': '4X4N|Sequence|G|U|30',
            'unit_id': '4X4N|1|G|U|30||A',
            'seq_unit': 'U',
            'number': 30,
            'index': 29,
        }
        assert self.mapping['4X4N|1|G|C|30||B'] == {
            'unit_id': '4X4N|1|G|C|30||B',
            'seq_id': '4X4N|Sequence|G|U|30',
            'seq_unit': 'U',
            'number': 30,
            'index': 29
        }


class WithAltidsTest(ExperimentalSequenceMappingTest):
    name = '2G32'
    chains = 'L'

    def test_it_builds_all_mappings(self):
        self.assertEquals(16, len(self.data))

    def test_it_build_mappings_using_alt_ids(self):
        assert self.mapping['2G32|1|L|0C|90||A'] == {
            'seq_id': '2G32|Sequence|L|0C|90',
            'seq_unit': '0C',
            'unit_id': '2G32|1|L|0C|90||A',
            'number': 90,
            'index': 0
        }


class WithInsertionCodes(ExperimentalSequenceMappingTest):
    name = '4V6F'
    chains = 'AB'

    def test_maps_strange_residue(self):
        assert '4V6F|1|AB|A|89|||A' in self.mapping

    def test_it_maps_all_units(self):
        assert len(self.data) == 122

    def test_it_maps_the_correct_order(self):
        assert self.mapping['4V6F|1|AB|A|1|||M'] == {
            'unit_id': '4V6F|1|AB|A|1|||M',
            'index': 0,
            'seq_id': '4V6F|Sequence|AB|A|1|||M',
            'seq_unit': 'A',
            'number': 1,
        }
