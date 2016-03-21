from unittest import TestCase

import pytest

from fr3d.unit_ids import decode
from fr3d.unit_ids import encode
from fr3d.unit_ids import InvalidUnitId


class ResidueUnitIDTest(TestCase):

    def test_decodes_a_full_residue_id(self):
        val = decode('2AVY|1|A|C|50||||1_555')
        ans = {
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': None,
            'alt_id': None,
            'insertion_code': None,
            'symmetry': '1_555'
        }
        self.assertEquals(val, ans)

    def test_decodes_short_with_default_data(self):
        val = decode('2AVY|1|A|C|50')
        ans = {
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': None,
            'alt_id': None,
            'insertion_code': None,
            'symmetry': '1_555'
        }
        self.assertEqual(val, ans)

    def test_encodes_full_residue_id(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': None,
            'alt_id': None,
            'insertion_code': None,
            'symmetry': '1_555'
        }, full=True)
        ans = '2AVY|1|A|C|50||||1_555'
        self.assertEqual(val, ans)

    def test_encodes_residue_id_with_non_default_symmetry(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'symmetry': '6_555'
        })
        ans = '2AVY|1|A|C|50||||6_555'
        self.assertEqual(val, ans)

    def test_encodes_residue_id_without_sequence(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_number': 50,
            'symmetry': '6_555'
        })
        ans = '2AVY|1|A||50||||6_555'
        self.assertEqual(val, ans)

    def test_encodes_short_residue_id(self):
        val = encode({
            'pdb': '2AVY',
            'model': '1',
            'chain': 'A',
            'component_id': 'C',
            'component_number': '50',
        })
        ans = '2AVY|1|A|C|50'
        self.assertEqual(val, ans)

    def test_encodes_short_residue_id_without_sequence(self):
        val = encode({
            'pdb': '2AVY',
            'model': '1',
            'chain': 'A',
            'component_number': '50',
        })
        ans = '2AVY|1|A||50'
        self.assertEqual(val, ans)

    def test_encodes_residue_id_with_numbers(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
        })
        ans = '2AVY|1|A|C|50'
        self.assertEqual(val, ans)

    def test_encodes_short_residue_id_with_default_symmetry(self):
        val = encode({
            'pdb': '2AVY',
            'model': '1',
            'chain': 'A',
            'component_id': 'C',
            'component_number': '50',
            'symmetry': '1_555'
        })
        ans = '2AVY|1|A|C|50'
        self.assertEqual(val, ans)

    def test_encodes_full_residue_id_using_defaults(self):
        val = encode({
            'pdb': '2AVY',
            'model': '1',
            'chain': 'A',
            'component_id': 'C',
            'component_number': '50',
        }, full=True)
        ans = '2AVY|1|A|C|50||||1_555'
        self.assertEqual(val, ans)

    def test_round_trips_an_id(self):
        ans = '2AVY|1|A|C|50'
        val = encode(decode(ans))
        self.assertEqual(val, ans)

    @pytest.mark.skip()
    def test_fails_encoding_residue_id_missing_pdb(self):
        self.assertRaises(InvalidUnitId, encode, {
            'model': '1',
            'chain': 'A',
            'component_id': 'C',
            'component_number': '50',
        })

    @pytest.mark.skip()
    def test_fails_encoding_residue_id_missing_model(self):
        self.assertRaises(InvalidUnitId, encode, {
            'pdb': '2AVY',
            'chain': 'A',
            'component_id': 'C',
            'component_number': '50',
        })

    @pytest.mark.skip()
    def test_fails_encoding_residue_id_missing_chain(self):
        self.assertRaises(InvalidUnitId, encode, {
            'pdb': '2AVY',
            'model': '1',
            'component_id': 'C',
            'component_number': '50',
        })

    def test_fails_encoding_residue_id_missing_component_number(self):
        self.assertRaises(InvalidUnitId, encode, {
            'pdb': '2AVY',
            'model': '1',
            'chain': 'A',
            'component_id': 'C',
        })


class AtomUnitIDTest(TestCase):
    def test_can_generate_full_atom_id(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': "C1'",
            'alt_id': '',
            'insertion_code': '',
            'symmetry': '1_555'
        }, full=True)
        ans = "2AVY|1|A|C|50|C1'|||1_555"
        self.assertEqual(val, ans)

    def test_can_generate_short_atom_id(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': "C1'",
            'alt_id': '',
            'insertion_code': '',
            'symmetry': '1_555'
        })
        ans = "2AVY|1|A|C|50|C1'"
        self.assertEqual(val, ans)

    def test_can_generate_full_atom_id_missing_stuff(self):
        val = encode({
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': "C1'",
            'symmetry': '6_555'
        }, full=True)
        ans = "2AVY|1|A|C|50|C1'|||6_555"
        self.assertEqual(val, ans)

    def test_generate_short_with_missing_and_extra(self):
        val = encode({
            'pdb': '1GID',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': "C1'",
            'symmetry': '6_555',
            'x': -1,
            'y': 0,
            'z': 0,
        })
        ans = "1GID|1|A|C|50|C1'|||6_555"
        self.assertEqual(val, ans)

    def test_decodes_short_atom_id(self):
        val = decode("2AVY|1|A|C|50|C1'")
        ans = {
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': "C1'",
            'alt_id': None,
            'insertion_code': None,
            'symmetry': '1_555'
        }
        self.assertEqual(val, ans)

    def test_decodes_full_atom_id(self):
        val = decode("2AVY|1|A|C|50|C1'|A|c|6_555")
        ans = {
            'pdb': '2AVY',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'atom_name': "C1'",
            'alt_id': 'A',
            'insertion_code': 'c',
            'symmetry': '6_555'
        }
        self.assertEqual(val, ans)
