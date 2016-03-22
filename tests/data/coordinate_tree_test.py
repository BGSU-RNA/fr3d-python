import functools as ft
from unittest import TestCase

from fr3d.data.atoms import Atom
from fr3d.data.base import CoordinateTree
from fr3d.data.components import Component

import pytest


class CoordinateTreeTest(TestCase):

    def first(self):
        yield (
            Component([Atom(x=0.0, y=2.0, z=0.0),
                       Atom(x=0.0, y=2.0, z=1.0)],
                      pdb='0FJG', model=1, chain='A', number=3),
            [0.0, 2.0, 0.0]
        )

        yield (
            Component([Atom(x=0.0, y=2.0, z=-1.0),
                       Atom(x=0.0, y=2.0, z=-2.0)],
                      pdb='0FJG', model=1, chain='A', number=4),
            [0.0, 2.0, -1.0]
        )

    def second(self):
        yield (
            Component([Atom(x=0.0, y=0.0, z=0.0),
                       Atom(x=0.0, y=0.0, z=1.0)],
                      pdb='0FJG', model=1, chain='B', number=3),
            [0.0, 0.0, 0.0]
        )

        yield (
            Component([Atom(x=0.0, y=0.0, z=-1.0),
                       Atom(x=0.0, y=0.0, z=-2.0)],
                      pdb='0FJG', model=1, chain='B', number=4),
            [0.0, 0.0, -1.0]
        )

    def pairs(self, tree, *args):
        return [(x.unit_id(), y.unit_id()) for x, y in tree.pairs(*args)]

    def neighbors(self, tree1, tree2, *args, **kwargs):
        pairs = tree1.neighbors(tree2, *args, **kwargs)
        return [(x.unit_id(), y.unit_id()) for x, y in pairs]

    def setUp(self):
        self.tree1 = CoordinateTree(self.first())
        self.tree2 = CoordinateTree(self.second())

    def test_it_can_count_number_of_neighbors(self):
        assert self.tree1.count_neighbors(self.tree2, 0.1) == 0
        assert self.tree1.count_neighbors(self.tree2, 2.0) == 2

    def test_it_can_get_pairs(self):
        assert self.pairs(self.tree1, 0.0) == []
        assert self.pairs(self.tree1, 1.0) == [('0FJG|1|A||3', '0FJG|1|A||4')]

    def test_it_can_get_neighbors_between_trees(self):
        neighbors = ft.partial(self.neighbors, self.tree1, self.tree2)
        assert neighbors(0.0) == []
        assert self.neighbors(1.0) == [('0FJG|1|A|3', '0FJG|1|B|4')]

    @pytest.mark.skip()
    def test_it_can_get_only_unique_pairs(self):
        pass


class RealDataCoordinateTreeTest(TestCase):

    def setUp(self):
        pass

    @pytest.mark.skip()
    def test_it_can_build_a_tree_for_residues(self):
        pass

    @pytest.mark.skip()
    def test_it_can_build_a_tree_for_atoms(self):
        pass
