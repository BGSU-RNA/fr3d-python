import itertools as it

from fr3d.data.base import EntitySelector
from fr3d.data.pairs import Pairs
from fr3d.unit_ids import encode


class Container(object):

    def __init__(self, residues, pdb=None, model=None, chain=None,
                 symmetry=None, breaks=None):

        """Create a new Container.

        :residues: A list of Components that represent the residues in this
        Structure.
        :pdb: The pdb this container is a part of.
        :model: The model number of this container.
        :chain: The chain name of this container.
        :symmetry: The symmetry operator for this container.
        """

        self.pdb = pdb
        self.model = model
        self.chain = chain
        self.symmetry = symmetry
        self._residues = residues
        self._breaks = breaks

    def models(self):
        """Create an iterator over all models in this Container.
        """
        pass

    def symmetry(self, operator):
        """Create a new Container with the given symmetry operator.
        """

        return Container(self.residues(symmetry=operator, pdb=self.pdb),
                         pdb=self.pdb, symmetry=self.symmetry,
                         breaks=self.breaks)

    def model(self, model):
        """Get a model by number. The number is the same as in the cif file.
        Will raise an IndexException if asking for an unknown model.

        :model: Integer for the model to get.
        :returns: A model.
        """

        return Container(self.residues(model=model, pdb=self.pdb),
                         pdb=self.pdb, model=model, breaks=self.breaks)

    def chains(self):
        """Get an iterator over all chains in in this Structure.

        :returns: An iterator for all chains.
        """
        for model in self.models():
            for chain in model.chains:
                yield Container(self.residues(model=model, chain=chain),
                                pdb=self.pdb, model=model, chain=chain,
                                breaks=self.breaks)

    def chain(self, model_number, chain_id):
        """Get a specific chain.

        :model_number: The model number to use.
        :chain_id: The chain id to get.
        :returns: A Chain or None if no chain is present.
        """
        return self.select(model=model_number, chain=chain_id)

    def residues(self, **kwargs):
        """Get residues from this structure. The keyword arguments work as
        described by EntitySelector.

        :kwargs: Keywords for filtering and ordering
        :returns: The requested residues.
        """
        if 'polymeric' not in kwargs:
            kwargs['polymeric'] = True
        if kwargs.get('polymeric', False) is None:
            kwargs.pop('polymeric')

        return EntitySelector(self._residues, **kwargs)

    def infer_hydrogens(self):
        """ Infers hydrogen atoms for all bases.
        """
        for residue in self._residues:
            residue.infer_hydrogens()

    def atoms(self, **kwargs):
        for residue in self._residues:
            for atom in residue.atoms(**kwargs):
                yield atom

    def residue(self, unit_id):
        pass

    def select(self, **kwargs):
        """
        """

        kwargs['pdb'] = self.pdb
        data = {}
        for key in ['pdb', 'model', 'chain']:
            data[key] = kwargs.get(key)
        data['breaks'] = self._breaks
        return Container(self.residues(**kwargs), **data)

    def pairs(self, first={}, second={}, distance={}):
        """Create an iterator for all pairs in the structure. The first and
        second arguments are selectors, that is a dictionary that can contain
        the same arguments given to the residues method. The third argument
        distance is a dictionary

        :first: A selector for the first item in the pair.
        :second: A selector for the second item in the pair.
        :distance: A
        :returns: A iterator that will go over all matching pairs in the
        structure.
        """

        pairs = Pairs(self)
        pairs.first(**first)
        pairs.second(**second)
        pairs.distance(**distance)
        return pairs

    def unit_id(self):
        return encode({
            'pdb': self.pdb,
            'model': self.model,
            'chain': self.chain,
            'symmetry': self.symmetry
        })

    @property
    def sequence(self):
        if self._sequence is None:
            self._sequence = [r['residue'] for r in self.residue_iterator()]
        return self._sequence

    def __len__(self):
        """Compute the length of this Container. That is the number of residues
        in this structure.

        :returns: The number of atoms.
        """
        return sum(it.imap(lambda _: 1, self._residues))

    def __bool__(self):
        return bool(self._residues)

    def __repr__(self):
        return '<Container: %s>' % self.pdb
