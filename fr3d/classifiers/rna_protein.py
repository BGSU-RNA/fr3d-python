from fr3d import definitions as defs
from fr3d.classifiers.generic import Classifier as BaseClassifier


class Classifier(BaseClassifier):
    """A classifier for RNA protein interactions.
    """

    def __init__(self):
        first = {'sequence': defs.RNAbaseheavyatoms.keys()}
        second = {'sequence': defs.aa_backconnect.keys()}
        distance = {'use': 'center', 'cutoff': 4.0}
        super(Classifier, self).__init__(first=first, second=second,
                                         distance=distance)

    def classification(self, first, second):
        pass
