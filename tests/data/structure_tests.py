import unittest as ut
from random import shuffle

from fr3d.data import Structure

from tests import define_1S72_nucleotides as nts


class BasicTest(ut.TestCase):
    def test_empty_structure_is_false(self):
        val = Structure([])
        self.assertFalse(bool(val))

    def test_non_empty_structure_is_true(self):
        val = Structure([1])
        self.assertTrue(bool(val))

    def test_length_is_number_of_residues(self):
        val = Structure([1])
        self.assertEquals(1, len(val))

    def test_no_residues_is_0_length(self):
        val = Structure([])
        self.assertEquals(0, len(val))


class UnitIdTest(ut.TestCase):
    def test_has_a_unit_id(self):
        val = Structure([], pdb='1S72').unit_id()
        self.assertEquals('1S72', val)

    def test_selecting_a_subset_updates_unit_id(self):
        val = Structure([], pdb='1S72').select(model=1).unit_id()
        self.assertEquals('1S72|1', val)


class ResiduesTest(ut.TestCase):
    def setUp(self):
        self.structure = Structure([nts.nt77_9, nts.nt78_9, nts.nt79_9,
                                    nts.nt80_9], pdb="1S72")

    def test_can_get_sequence_from_all_residues(self):
        self.assertEquals(['A', 'G', 'U', 'A'], self.structure.sequence)

    def test_can_iterate_over_all_residues(self):
        val = list(self.structure.residues())
        self.assertEquals([nts.nt77_9, nts.nt78_9, nts.nt79_9, nts.nt80_9],
                          val)

    def test_can_select_residues(self):
        val = list(self.structure.residues(sequence='A'))
        self.assertEquals([nts.nt77_9, nts.nt80_9], val)

    def test_skip_iterating_over_nonpolymeric_by_default(self):
        self.structure._residues.append(nts.nt212_0)
        val = list(self.structure.residues())
        self.assertEquals([nts.nt77_9, nts.nt78_9, nts.nt79_9, nts.nt80_9],
                          val)

    def test_giving_polymeric_is_none_iterates_overall(self):
        self.structure._residues.append(nts.nt212_0)
        val = list(self.structure.residues(polymeric=None))
        ans = [nts.nt77_9, nts.nt78_9, nts.nt79_9, nts.nt80_9, nts.nt212_0]
        self.assertEquals(ans, val)


class FindingAResidueTest(ut.TestCase):
    def setUp(self):
        self.structure = Structure([nts.nt77_9, nts.nt78_9, nts.nt79_9,
                                    nts.nt80_9], pdb="1S72")

    def test_can_find_residue_by_unit_id(self):
        val = self.structure.residue('1S72||9|U|79')
        self.assertEquals(nts.nt79_9, val)

    def test_can_find_residue_by_index(self):
        val = self.structure.residue(1)
        self.assertEquals(nts.nt78_9, val)

    def test_can_find_residues_by_negative_index(self):
        val = self.structure.residue(-1)
        self.assertEquals(nts.nt80_9, val)

    def test_can_find_by_unit_id_residue_even_when_unordered(self):
        shuffle(self.structure._residues)
        val = self.structure.residue('1S72||9|U|79')
        self.assertEquals(nts.nt79_9, val)

    def test_raises_exception_if_missing_unit_id(self):
        self.assertRaises(IndexError, lambda: self.structure.residue('bob'))

    def test_raises_exception_if_missing_index(self):
        self.assertRaises(IndexError, lambda: self.structure.residue(1000000))
