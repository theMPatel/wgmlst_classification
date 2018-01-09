import os
import sys
import warnings
from collections import defaultdict

class DistMatrix(object):

    def __init__(self):
        self._entry_to_index = {}
        self._distances = []

    def add(self, entry, dists):

        if entry in self._entry_to_index:
            warnings.warn( 'Tried to add distances for an'
                'entry already present in the matrix', RuntimeWarning )
            return

        self._distances.append(dists)
        self._entry_to_index[entry] = len(self._distances)-1


    def __getitem__(self, tup):

        if len(tup) != 2:
            raise RuntimeError('Tried to access too many'
                ' distances: {}'.format(', '.join(tup)))

        for entry in tup:
            if entry not in self._entry_to_index:
                raise RuntimeError('Missing key: {}'.format(
                    entry))

        index_zero = self._entry_to_index[tup[0]]
        index_one = self._entry_to_index[tup[1]]

        if index_zero == index_one:
            if tup[0] == tup[1]:
                return 0.
            else:
                raise RuntimeError('Unrecoverable error'
                    ' two different isolates have the same'
                    ' index in distance_matrix')

        elif index_zero < index_one:
            return self._distances[index_one][index_zero]

        else:
            return self._distances[index_zero][index_one]

