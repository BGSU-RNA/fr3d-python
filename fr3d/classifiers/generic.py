import abc


class Classifier(object):
    """This is a base classifier for all classifiers.
    This implements the generic behavior needed for dealing with a single
    structure and leaves the specifics for dealing with a pair up to each
    specific classifier.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, first=None, second=None, distance=None):
        """Create a new classifier. If no selectors are given then this will
        default to selecting all pairs regardless of distance or other
        properties.

        :first: The selector for the first set of components.
        :second: The selector for the second set of components.
        :distance: The distance criteria for the pairs.
        """

        self.first = first or {}
        self.second = second or {}
        self.distance = distance or {}

    @abc.abstractmethod
    def classification(self, first, second):
        """Compute the classification for a specific pair of components.

        :first: The first component in the pair.
        :second: The second component in the pair
        :returns: The classficiation or a false value if there is none.
        """
        pass

    def classify(self, structure):
        """Classify all pairs in a structure. This will select all pairs which
        match the given first, second and distance filters and then attempt to
        classify each one. If we cannot classify them or if there is no
        classification for the pair it is ignored. Otherwise it is added to a
        list of tuples. Each entry in the list is (component1, component2,
        classification).

        :structure: A structure object to classify all pairs of.
        :returns: A list of tuples of the classification.
        """

        pairs = structure.pairs(first=self.first, second=self.second,
                                distance=self.distance)
        classified = []
        for (first, second) in pairs:
            try:
                classification = self.classification(first, second)
                if classification:
                    classified.append((first, second, classification))
            except:
                pass

        return classified
