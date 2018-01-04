# 
# Module providing the graph functionality. Serves as the basis for the nomenclature
# 
# nomenclature/graph.py
# 
# 
# 


# 
# Imports
# 

import sys
import os
from collections import defaultdict

__all__ = [ 'Tree', 'Node', 'NamedNode' ]

# 
# Class that holds the graph entrance and helper functions
# 

class Tree( object ):
    
    CURRENT = None
    DEPTH = -1

    def __init__(self, depth):
        
        Tree.CURRENT = self
        Tree.DEPTH = depth
        self._tree = Node( 1, 0, None )
        self._names = defaultdict(list)

    def BuildTree( self ):
        lengthMismatch = 0
        for key, address in self._names.items():
            if len( address ) == Tree.DEPTH:
                
                currentNode = self._tree
                for i in address[:-1]:
                    nextNode = currentNode.GetChild( i )
                    
                    if nextNode is None:
                        nextNode = Node(i, currentNode.Level()+1, currentNode)
                        currentNode.AddChild( i, nextNode )
                    
                    currentNode = nextNode

                nextNode = currentNode.GetChild( address[-1] )
                if nextNode is None:
                    nextNode = NamedNode(address[-1], currentNode.Level()+1,
                        currentNode)

                    nextNode.AddChild( key )
                    currentNode.AddChild( address[-1], nextNode )

                else:
                    nextNode.AddChild( key )
            
            else:
                lengthMismatch += len( address )

        if lengthMismatch > 0:
            pass

    def FinalizeName( self, key, name ):
        node = self.Traverse(name)
        assert( node.IsNamed( key ) )
        self._names[ key ] = name

    def GetName( self, key ):
        return self._names.setdefault(key, [])

    def GetNames(self):
        return self._names

    def GetPart( self, key, level ):
        return self.GetName( key )[ :level ]

    def GetStrName( self, key ):
        return Tree.NameToStr( self.GetName(key) )

    def HasName( self, key ):
        return len(self.GetName(key)) == Tree.DEPTH

    def Load(self, flobj):
        database = json.load( flobj )
        assert isinstance( database, dict )
        self._names.update( database )

    @staticmethod
    def NameToStr( parts ):
        assert( isinstance( parts, list) or \
            isinstance( parts, tuple ) )

        return '.'.join( map(str, parts) )

    @staticmethod
    def NameFromStr( name ):
        try:
            assert( isinstance( name, str ) )
        except AssertionError:
            print( name )

        return list( map(int, name.split('.') ) )

    def RemoveName(self, key ):
        location = self._names.get( key, None )
        if location is not None:
            self.Traverse( location ).RemoveNamed( key )
            del self._names[key]

    def Save(self, flobj):
        json.dump(self._names, flobj)

    def Traverse( self, ids ):
        assert( isinstance( ids, tuple ) or \
            isinstance( ids, list ) )
        return self._tree.Traverse( ids )

    def Tree(self):
        return self._tree

    def __len__(self):
        return len( self._names )

class Node( object ):

    def __init__( self, ID, level, parent ):

        self._ID = ID
        self._level = level
        self._parent = parent
        self._children = {}

    def AddChild( self, ID, node ):
        assert self._children.get(ID, None) is None
        node.SetParent( self )
        self._children[ID] = node

    def DeleteChild( self, ID ):
        del self._children[ID]

    def DFSNamed(self):

        if isinstance( self, NamedNode ):
            yield self

        else:
            toCheck = list( self.GetChildrenNodes() )

            while( len(toCheck) ):

                ob = toCheck.pop()
                if isinstance( ob, NamedNode ):
                    yield ob
                else:
                    toCheck.extend( ob.GetChildrenNodes() )

    def GetChild(self, ID):
        return self._children.get( ID, None )

    def GetChildren(self):
        return self._children

    def GetChildrenKeys(self):
        return self._children.keys()

    def GetChildrenNodes(self):
        return self._children.values()

    def ID(self):
        return self._ID

    def Level(self):
        return self._level

    def NewChildNode( self ):

        if not len(self._children):
            nextCluster = 1
        
        else:
            nextCluster = max( self._children.keys() ) + 1
        
        if self._level == Tree.DEPTH - 1:
            self._children[ nextCluster ] = NamedNode( nextCluster,
                self._level+1, self )
        else:
            self._children[ nextCluster ] = Node( nextCluster, 
                self._level+1, self )
        
        return self.GetChild( nextCluster )

    def MergeNodes(self, nodes):

        allNodes = [ self.GetChild(i) for i in nodes ]
        minObjs = {}
        
        for node in allNodes:
            size = sum( n.Size() for n in node.DFSNamed() )
            minObjs[ node ] = size

        maxObj = max( minObjs, key=lambda x: minObjs[x] )
        del minObjs[ maxObj ]

        if self._level == Tree.DEPTH - 1:
            
            for node in minObjs:
                for key in node.GetChildrenKeys():
                    maxObj.AddChild( key )

                self.DeleteChild( node.ID() )
        else:

            toStart = max( maxObj.GetChildrenKeys() ) + 1
            for node in minObjs:

                for i, child in enumerate( node.GetChildrenNodes(), toStart ):
                    child.SetID( i )
                    maxObj.AddChild( i, child )

                toStart = i + 1

                self.DeleteChild( node.ID() )

        for node in maxObj.DFSNamed():
            wgst = node.ValidateNames()
            for key in node.GetChildrenKeys():
                Tree.CURRENT.FinalizeName( key, wgst )

        return maxObj

    def SetID(self, ID):
        assert isinstance( ID, int )
        self._ID = ID

    def SetParent( self, parent ):
        assert isinstance( parent, Node )
        self._parent = parent

    def Size(self):
        return len( self._children )

    def Traverse( self, ids ):
        if len( ids ) == 1:
            return self.GetChild( ids[0] )
        
        else:
            child = self.GetChild( ids[0] )
            
            if child:
                return child.Traverse( ids[1:] )
            else:
                return None

    def NTraverse( self ):
        if self._parent == None:
            if self.Size() == 1:
                return [-1]
            else:
                return [self.ID()]
        else:
            if self.Size() == 1:
                return [-1] + self._parent.NTraverse()
            else:
                return [ self.ID() ] + self._parent.NTraverse()

    def RTraverse( self ):
        if self._parent == None:
            return []
        else:
            return [self._ID] + self._parent.RTraverse()

    def __repr__(self):
        return ( 
            'Class: {} | ID: {} | Level: {} ' 
            '| Parent: {}'.format(
                                    type(self),
                                    self._ID,
                                    self._level,
                                    self._parent.ID() 
                                )
            )

class NamedNode( Node ):

    def __init__(self, ID, level, parent ):
        super(self.__class__, self).__init__( ID, level, parent )
        self._wgst = self.RTraverse()[-1::-1]

    def AddChild(self, key):
        self._children[ key ] = self._wgst
        return self.GetChild( key )
    
    def IsNamed(self, key):
        return key in self._children

    def RemoveNamed( self, key ):
        del self._children[ key ]

    def ValidateNames( self ):
        for i, val in enumerate( self.RTraverse()[-1::-1] ):
            self._wgst[i] = val

        return self._wgst

    def Address( self ):
        return self._wgst