##############################################################
# Nomenclature naming script
# Author: Hannes Pouseele, Milan Patel
# Last Modified: 08/16/2017
# Version: 2.0
# Contact: hannes_pouseele@applied-maths.com, mpatel5@cdc.gov
##############################################################
import os
import sys
from datetime import datetime
from collections import defaultdict
from queue import PriorityQueue
import pickle
import time

import json
import csv
import random
import numpy as np

from wgst.database import Database
from .tree import *
from tqdm import *

# ======================== DISTANCE FUNCTION ================================#
def GetDistance(p1, p2):
    
    common = np.multiply(p1>0, p2>0)
    nCommon = np.sum(common)
    nSame = np.sum(p1[common]==p2[common])
    if nCommon:
        return 100.0 * (float(nCommon) - float(nSame)) / float(nCommon) #nCommon-nSame
    else:
        return 100.0

# =========================== DATABASE OBJECT ===============================#
class Database(Database):

    CHILDREN = []

    def __init__( self, organism, scheme, metadata, outdir, **kwargs ):
        super().__init__( organism, scheme, metadata, outdir, **kwargs )
        Database.CHILDREN.append( self )

    def CreateSubset( self ):

        def Bisect():

            print('Bisecting database by date')

            allSet = list( self.GetEntries() )
            allSet.sort( key=lambda x: x._uploaddate )
            startingSet = allSet[ :(len(allSet)//2) ]
            addingSet = allSet[ (len(allSet)//2): ]
            addingSet = Chunkify( addingSet )

            return (( startingSet, addingSet ))

        def Chunkify(li):
            return [ li[i:i+30] for i in range(0, len(li), 30) ]

        if not len( self._test ):
            return Bisect()

        else:
            return self._train, Chunkify( self._test )

#=========================== WGST HISTORY CLASS =============================#
class WgstHistory(object):

    def __init__( self, thresholds, outdir, seedSize, addingSize ):
        
        self._thresholds = thresholds
        self._outdir = outdir
        self._seedSize = seedSize
        self._addingSize = addingSize
        self._attrs = [ '_totalChanges', '_totalObChanges', '_seedChanges', '_seedObChanges', '_postChanges', '_postObChanges' ]

        for attribute in self._attrs:
            setattr(self, attribute, list() )

    def CalculateDifference( self, oldVal, newVal ):

        old = list(map( int, oldVal.split('.') ))
        new = list(map( int, newVal.split('.') ))

        assert( len(old) == len(new) )

        theSame = 0

        for a, b in zip( old, new ):
            if a == b:
                theSame += 1
            else:
                break

        return theSame

    def Flatten(self, arr):
        tempList = [0]*len(self._thresholds)

        for i in range(len(arr[0])):
            for j in arr:
                tempList[i] += j[i]

        return tempList

    def GetTotal(self):
        return sum( self._totalChanges )

    def TrackChanges( self, oldVal, newVal, iteration, ob = False, seed = False ):

        theSame = self.CalculateDifference( oldVal, newVal )
        curr = { attr:[0]*len(self._thresholds) for attr in self._attrs }
        
        curr['_totalChanges'][ theSame ] += 1

        if ob:
            curr['_totalObChanges'][ theSame ] += 1
        
        if seed:

            if ob:

                curr['_seedObChanges'][ theSame ] += 1

            curr['_seedChanges'][ theSame ] += 1

        else:

            if ob:
                curr['_postObChanges'][ theSame ] += 1

            curr['_postChanges'][ theSame ] += 1

        for attribute, changes in curr.items():
            self.__getattribute__(attribute).append( changes )

    def Save(self):

        OUT_DIR = os.path.join( self._outdir, 'wgstHistory' )
        if not os.path.isdir( OUT_DIR ):
            os.mkdir( OUT_DIR )

        OUT_NAME = 'wgst_history' + '-{}'.format( datetime.now().strftime("%m-%d-%y@%H-%M-%S") )
        PATH = os.path.join( OUT_DIR, '{}.csv'.format( OUT_NAME ) )

        with open( PATH , 'w' ) as f:
            writer = csv.writer(f)

            oArray = [['Seed size:', self._seedSize, 'Adding size:', self._addingSize], [], self._thresholds]

            for attribute in self._attrs:
                tmp = [[attribute],[]]
                tmp[0].extend( self.Flatten( self.__getattribute__(attribute) ))
                tmp.extend( self.__getattribute__(attribute) )

                oArray.extend( tmp )

            writer.writerows( oArray )

#=========================== NAMING FUNCTION ================================#
def CalcName(named, tree, unNamedEntry, distances, thresholds ):

    # Make sure thresholds are sorted biggest first
    thresholds.sort(key=lambda x: -x)

    # Let's get the head of the tree
    headNode = tree.Tree()

    # UnNamed should have no name
    patternName = tree.GetName( unNamedEntry )
    assert len( patternName ) == 0

    # Will serve as the node place holder as we drill down
    # into the proper home for a our unNamedEntry
    currentNode = headNode

    # Create a local copy of the named list
    # namedEntries = list( named )    
    for level, threshold in enumerate( thresholds, 1 ):

        # Create the clusters at this level:
        closestClusters = set()
        toRemove = []
        for i, (dist, namedIndex) in enumerate( distances ):
            
            if dist <= threshold:
                keyToGet = named[namedIndex]
                partial = tuple( tree.GetPart( keyToGet, level ) )
                closestClusters.add( partial )

            else:
                toRemove.append( i )

        # If something isn't similar at this level
        # There's no way it will be similar at the lower
        # levels, let's remove those
        for index in reversed(toRemove):
            del distances[ index ]

        # Assigned to only one cluster
        if len( closestClusters ) == 1:
            pattern = closestClusters.pop()
            patternName.append( pattern[-1] )
            currentNode = headNode.Traverse( pattern )

        # New cluster!
        elif len( closestClusters ) == 0:
            currentNode = currentNode.NewChildNode()
            patternName.append( currentNode.ID() )

        # Merge time:
        else:
            toMerge =[ c[-1] for c in closestClusters ]
            currentNode = currentNode.MergeNodes( toMerge )
            patternName.append( currentNode.ID() )

    assert isinstance( currentNode, NamedNode )

    currentNode.AddChild( unNamedEntry )

    assert currentNode.GetChild( unNamedEntry ) == patternName

    return tree

#================== CALCULATOR OBJECT =======================================#
class Calculator(object):
    
    def __init__( self, args ):

        self._organism = args['organism']
        self._thresholds = args['thresholds']
        self._thresholds.sort( key=lambda x: -x )
        self._scheme = args['scheme']
        self._minPres = args['minPres']
        self._outdir = args['outdir']

        self._entryBase = Database( 
            self._organism,
            self._scheme,
            args['fields_path'],
            self._outdir,
            locisize=args['locisize'],
            schemespath=args['views_path']
        )

        self._entryBase.LoadCalls()
        self._entryBase.QC( self._minPres )

    def DoValidation(self):

        self._tree = Tree(len(self._thresholds))
        
        # Create starting and adding sets return the updated Names object to track history
        self._startingSet, addingSet = self._entryBase.CreateSubset()

        namedEntries = []

        self.DoCalc( self._startingSet, self._tree, namedEntries, 'Initializing' )

        self._wgstHistory = WgstHistory( self._thresholds, self._outdir, \
            len(self._startingSet), self._entryBase.TotalIsolates() - \
            len(self._startingSet) )


        for i in trange( len( addingSet), desc='Clustering' ):
            self.DoCalc( addingSet[i], self._tree, namedEntries, 'Incrementing', i+1 )

        self._wgstHistory.Save()
        self._entryBase.Save()

    def DoCalc(self, selection, tree, namedEntries, status, iteration = 0):

        nameFrequencies = {}

        for entry in namedEntries:
            nameFrequencies[ tree.GetStrName(entry) ] = nameFrequencies.get( tree.GetStrName(entry), 0) + 1

        #for each entry, calculate the name 
        for entry in tqdm( selection, desc=status ):
            
            qcStatus = entry._qc

            if not qcStatus:
                continue
                
            if entry._key in namedEntries:
                continue
            
            #calculate the distance between the unnamed sample and all the named samples
            dists = []
            for i in range(len(namedEntries)):
                dist = GetDistance( entry._allelecalls,
                    self._entryBase.GetEntry(namedEntries[i])._allelecalls )

                if dist <= self._thresholds[0]:
                    dists.append((dist, i))

            # dists = [ GetDistance( 
            #     entry._allelecalls, 
            #     self._entryBase.GetEntry( e )._allelecalls ) for e in namedEntries ]
            
            #calculate the name of the entry
            CalcName( namedEntries, tree, entry._key, dists, self._thresholds )
            
            #keep track of the data
            if tree.HasName( entry._key ) and entry._key not in namedEntries:
                namedEntries.append( entry._key )
                nameFrequencies[ tree.GetStrName( entry._key ) ] = nameFrequencies.get( tree.GetStrName( entry._key ), 0) + 1


        if iteration == 0:
            for entry in namedEntries:
                self._entryBase.GetEntry( entry )._wgst = tree.GetStrName( entry  )

        for entry in namedEntries:
            dbEntry = self._entryBase.GetEntry( entry  )

            oldName = dbEntry._wgst

            newName = tree.GetStrName( dbEntry._key )

            if oldName == newName:
                continue
                
            if iteration > 0 and oldName and nameFrequencies[oldName] > 1:
                ob = bool( dbEntry._outbreak )
                self._wgstHistory.TrackChanges( oldName, newName, iteration, ob = ob, seed = entry in self._startingSet )

            dbEntry._wgst = newName
            dbEntry._wgstHistory.append( newName )

def Main(args):

    try:
        calc = Calculator( args )
        calc.DoValidation()
    
    except:
        outname = 'state-{}.pickle'.format( datetime.now().strftime("%m-%d-%y@%H-%M-%S") )
        outfile = os.path.join( args['outdir'], 'state_data', outname )

        with open( outfile, 'wb' ) as f:
            pickle.dump( calc._wgstHistory, f )