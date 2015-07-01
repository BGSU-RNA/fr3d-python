from unittest import TestCase

from fr3d.data import Atom
from fr3d.data.base import EntitySelector


class BasicSelectionTest(TestCase):
    def setUp(self):
        self.atoms = [
            Atom(type='C', name='a1', component_number=3,
                 x=1.0, y=0.0, z=0.0, pdb='1S72', model=1, chain='A',
                 component_id='C'),
            Atom(type='C', name='a2', component_number=2,
                 x=2.0, y=0.0, z=0.0, pdb='1S72', model=1, chain='A',
                 component_id='C'),
            Atom(type='N', name='b1', component_number=1,
                 x=3.0, y=0.0, z=0.0, pdb='1S72', model=1, chain='A',
                 component_id='C'),
            Atom(type='N', name='c2', component_number=0,
                 x=0.0, y=1.0, z=0.0, pdb='1S72', model=1, chain='A',
                 component_id='C')
        ]

    def test_can_iterate_over_all(self):
        val = len(list(EntitySelector(self.atoms)))
        self.assertEquals(4, val)

    def test_gives_nothing_if_no_match(self):
        val = list(EntitySelector(self.atoms, name='bob'))
        self.assertEquals([], val)

    def test_can_select_by_one_name(self):
        val = list(EntitySelector(self.atoms, name='a1'))
        self.assertEquals([self.atoms[0]], val)

    def test_can_select_with_list_of_names(self):
        val = sorted(EntitySelector(self.atoms, name=['a1', 'c2']),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[0], self.atoms[3]], val)

    def test_can_select_with_tuple_of_names(self):
        val = sorted(EntitySelector(self.atoms, name=('a1', 'c2')),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[0], self.atoms[3]], val)

    def test_can_select_with_set_of_names(self):
        val = sorted(EntitySelector(self.atoms, name=set(['a1', 'c2'])),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[0], self.atoms[3]], val)

    def test_can_filter_by_function(self):
        val = sorted(EntitySelector(self.atoms, x=lambda x: x >= 2),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[1], self.atoms[2]], val)

    def test_can_filter_by_function_using_object(self):
        val = sorted(EntitySelector(self.atoms, _=lambda a: a.x >= 2),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[1], self.atoms[2]], val)

    def test_can_filter_by_several(self):
        val = sorted(EntitySelector(self.atoms, name=('a1', 'c2'), type='N'),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[3]], val)

    def test_it_does_not_fail_if_missing_attribute_used(self):
        val = sorted(EntitySelector(self.atoms, name3=('a1', 'c2'), type='N'),
                     key=lambda a: a.name)
        self.assertEquals([], val)

    def test_can_filter_by_method_to_find_one(self):
        comp_id = '1S72|1|A|C|1'
        val = list(EntitySelector(self.atoms, component_unit_id=comp_id))
        self.assertEquals([self.atoms[2]], val)

    def test_can_filter_by_method_to_find_several(self):
        ids = ['1S72|1|A|C|1', '1S72|1|A|C|3']
        val = sorted(EntitySelector(self.atoms, component_unit_id=ids),
                     key=lambda a: a.name)
        self.assertEquals([self.atoms[0], self.atoms[2]], val)
