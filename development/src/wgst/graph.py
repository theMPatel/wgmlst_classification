import sys
import os
import json
import threading
import queue

from collections import defaultdict

# 
# Graph entrance
# 

##############################################################################
class Tree( object ):

    CURRENT = {}

    def __init__( self, thresholds, organism, loci ):

        if Tree.CURRENT[organism] is None:
            Tree.CURRENT[organism] = self
        else:
            raise RuntimeError('This organism has already been'
                'initialized')
        
        self._organism = organism
        self._thresholds = thresholds
        self._depth = len(self._thresholds)
        self._tree = Node( -1, 0, self )
        self._names = defaultdict(list)
        self._loci_count = loci

    def build_tree( self ):
        # *** Need to implement *** #
        pass

    def finalize_name(self, key, name):
        node = self.traverse( name )
        assert( node.is_named( key ) )
        self._names[key] = node.address()

    def get_height(self):
        return self._depth

    def get_loci(self):
        return self._loci_count

    def get_name(self, key):
        return self._names.get(key, None)

    def get_names(self):
        return self._names

    def get_part(self, key, depth):
        toreturn = self.get_name(key)

        if toreturn:
            return toreturn[:depth]
        else:
            return toreturn
    
    def get_str_name(self, key, depth=None):
        if depth is not None:
            if isinstance(depth, int):
                return Tree.name_to_str(self.get_part(key, depth))
            else:
                raise ValueError("Cannot give partial using provided"
                    " depth of type {}".format(type(depth)))
        else:
            return Tree.name_to_str( self.get_name(key) )

    def get_thresholds(self):
        return self._thresholds

    def init_name(self, key):
        return self._names.setdefault(key, [])

    def is_final(self, key):
        return key in self._names and \
            len(self.get_name(key)) == self._depth

    def load(self, flobj):
        # *** Need to implement *** #
        pass

    def remove_name(self, key):
        # *** Need to implement *** #
        pass

    def traverse(self, address):
        # *** Need to implement *** #
        pass

    @staticmethod
    def name_to_str( parts ):
        assert( isinstance( parts, list) or \
            isinstance( parts, tuple ) )

        return '.'.join( map(str, parts) )

    @staticmethod
    def name_from_str( string ):
        assert( isinstance( name, str ) )
        return list( map(int, name.split('.') ) )

    def tree(self):
        return self._tree

    def __len__(self):
        return len(self._names)

##############################################################################

class Node( object ):
    
    def __init__( self, ID, level, parent ):
        self._ID = ID
        self._level = level
        self._parent = parent
        self._children = {}
        self._reference = [ set() for _ in range(self.get_loci()) ]

    def get_height(self):
        return self._parent.get_height()

    def get_loci(self):
        return self._parent.get_loci()

    def _process(self, isolate, indices, distance):
        assert distance <= len(self._reference)
        assert (100. - 100.*( distance / len(self._reference )))

        if len(indices) == 0:
            if len( self._children ) > 1 and level != self.get_height():
                raise RuntimeError('Major algorithm error'
                    ' invalid logic in relation to indices'
                    ' and height' )

            else:
                (ID, child), = self._children.items()
                child._process(isolate, indices, distance)




class NamedNode( Node ):

    def __init__(self, ID, level, parent):
        super(self.__class__, self).__init__( ID, level, parent )