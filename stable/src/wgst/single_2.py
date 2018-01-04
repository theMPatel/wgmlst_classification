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

import json
import csv
import random
import numpy as np

from wgst import database
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
class Dbase(database.DB):

    CHILDREN = []

    def __init__( self, scheme, fieldsPath, viewsPath, outdir, adding_file = None ):
        super().__init__( scheme, fieldsPath, viewsPath, outdir, adding_file )
        Dbase.CHILDREN.append( self )

    def CreateSubset( self ):

        def Bisect():

            print('Bisecting database by date')

            allSet = list( self.GetEntries().values() )
            allSet.sort( key=lambda x: x.UploadDate() )
            startingSet = allSet[ :(len(allSet)//2) ]
            addingSet = allSet[ (len(allSet)//2): ]
            addingSet = [ addingSet[i:i+30] for i in range(0, len(addingSet), 30) ]

            return (( startingSet, addingSet ))

        if self._addExternal:
            if len( self._addingSet ) > 0:
                print( 'Creating subset based on external adding set' )
                startingSet = [entry for entry in self._entries.values() if entry \
                    not in self._addingSet ]

                addingSet = [ self._addingSet[i:i+30] for i in range(0, len(self._addingSet), 30) ]

                return(( startingSet, addingSet ))
            else:
                print('External data is non-unique, bisecting database instead')
                return Bisect()

        else:
            return Bisect()

#=========================== WGST HISTORY CLASS =============================#
class WgstHistory(object):

    def __init__( self, thresholds, outdir, seedSize, addingSize ):
        
        self._thresholds = thresholds
        self._outdir = outdir
        self._seedSize = seedSize
        self._addingSize = addingSize

        self._seedChanges = [0]*len( self._thresholds )
        self._seedChangesWithTime = defaultdict(list)
        
        self._seedObChanges = [0]*len( self._thresholds )
        self._seedObChangesWithTime = defaultdict(list)

        self._postChanges = [0]*len( self._thresholds )
        self._postChangesWithTime = defaultdict(list)

        self._postObChanges = [0]*len( self._thresholds )
        self._postObChangesWithTime = defaultdict(list)

        self._totalChanges = [0]*len( self._thresholds )
        self._totalChangesWithTime = defaultdict(list)
        
        self._totalObChanges = [0]*len( self._thresholds )
        self._totalObChangesWithTime = defaultdict(list)

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

    def GetTotal(self):
        return sum( self._totalChanges )

    def TrackChanges( self, oldVal, newVal, iteration, ob = False, seed = False ):

        theSame = self.CalculateDifference( oldVal, newVal )
        
        self._totalChanges[ theSame ] += 1
        self._totalChangesWithTime[ iteration ].append( theSame )

        if ob:
            self._totalObChanges[ theSame ] += 1
            self._totalObChangesWithTime[ iteration ].append( theSame )
        
        if seed:

            if ob:

                self._seedObChanges[ theSame ] += 1
                self._seedObChangesWithTime[ iteration ].append( theSame )

            else:
                self._seedChanges[ theSame ] += 1
                self._seedChangesWithTime[ iteration ].append( theSame )

        else:

            if ob:
                self._postObChanges[ theSame ] += 1
                self._postObChangesWithTime[ iteration ].append( theSame )

            else:
                self._postChanges[ theSame ] += 1
                self._postChangesWithTime[ iteration ].append( theSame )


    def Save(self):

        OUT_DIR = os.path.join( self._outdir, 'wgstHistory' )
        if not os.path.isdir( OUT_DIR ):
            os.mkdir( OUT_DIR )

        OUT_NAME = '-'.join( map(str, self._thresholds ) ) + '-{}'.format( datetime.now().strftime("%m-%d-%y@%H-%M-%S") )
        PATH = os.path.join( OUT_DIR, '{}.csv'.format( OUT_NAME ) )

        with open( PATH , 'w' ) as f:
            writer = csv.writer(f)

            oArray = [['Seed size:', self._seedSize, 'Adding size:', self._addingSize]]
            oArray.append([])
            oArray.append( ['Total'] )
            oArray.append( self._thresholds )
            oArray.append( self._totalChanges )
            oArray.append( ['Total OB'] )
            oArray.append( self._totalObChanges )
            oArray.append( ['Seed'] )
            oArray.append( self._seedChanges )
            oArray.append( ['Seed OB'] )
            oArray.append( self._seedObChanges )
            oArray.append( ['Post'] )
            oArray.append( self._postChanges )
            oArray.append( ['Post OB'] )
            oArray.append( self._postObChanges )
            
            oArray.append( ['Total Over Time'])
            for iteration in self._totalChangesWithTime:
                temp = [0]*len(self._thresholds)

                for changes in self._totalChangesWithTime[iteration]:
                    temp[ changes ] += 1

                temp.insert(0, iteration)
                oArray.append( temp )

            oArray.append( ['Total OB Over Time'])
            for iteration in self._totalChangesWithTime:
                temp = [0]*len(self._thresholds)

                for changes in self._totalObChangesWithTime[iteration]:
                    temp[ changes ] += 1

                temp.insert(0, iteration)
                oArray.append( temp )

            oArray.append( ['Seed Over Time'])
            for iteration in self._seedChangesWithTime:
                temp = [0]*len(self._thresholds)

                for changes in self._seedChangesWithTime[iteration]:
                    temp[ changes ] += 1

                temp.insert(0, iteration)
                oArray.append( temp )

            oArray.append( ['Seed OB Over Time'])
            for iteration in self._seedObChangesWithTime:
                temp = [0]*len(self._thresholds)

                for changes in self._seedObChangesWithTime[iteration]:
                    temp[ changes ] += 1

                temp.insert(0, iteration)
                oArray.append( temp )
            
            oArray.append( ['Post Over Time'])
            for iteration in self._postChangesWithTime:
                temp = [0]*len(self._thresholds)

                for changes in self._totalChangesWithTime[iteration]:
                    temp[ changes ] += 1

                temp.insert(0, iteration)
                oArray.append( temp )
            
            oArray.append( ['Post OB Over Time'])
            for iteration in self._postObChangesWithTime:
                temp = [0]*len(self._thresholds)

                for changes in self._seedObChangesWithTime[iteration]:
                    temp[ changes ] += 1

                temp.insert(0, iteration)
                oArray.append( temp )

            writer.writerows( oArray )

#================== CALCULATOR OBJECT =======================================#
class Calculator(object):
    
    def __init__( self, args ):

        self.charViewId = args['scheme']
        self.thresholds = args['thresholds']
        self.thresholds.sort( key=lambda x: -x )
        self.minPresenceThreshold = args['minPres']
        self._args = args
        self._entryBase = Dbase( 
            self.charViewId,
            args['fields_path'],
            args['views_path'],
            args['outdir'],
            args['adding_file'] )

        self._entryBase.QC( self.minPresenceThreshold )

    def DoValidation(self):

        existingNames = Names()
        
        # Create starting and adding sets return the updated Names object to track history
        self._startingSet, addingSet = self._entryBase.CreateSubset()
        existingNames = self.DoCalc(self._startingSet, existingNames, 'Initializing' )

        self._wgstHistory = WgstHistory( self.thresholds, self._args['outdir'], \
            len(self._startingSet), self._entryBase.TotalIsolates() - \
            len(self._startingSet) )

        for i in trange( len( addingSet), desc='Clustering' ):
            existingNames = self.DoCalc(addingSet[i], existingNames, 'Incrementing', i+1 )

        self._entryBase.Save( clusterSearch = self._args['clustersearch'] )
        self._wgstHistory.Save()

    def DoCalc(self, selection, existingNames, status, iteration = 0):

        namedEntries = []
        nameFrequencies = {}
        
        for entry in self._entryBase.GetEntries().values():
            if existingNames.HasResolvedName( entry.Key() ):
                namedEntries.append( entry.Key() )
                nameFrequencies[ entry.Wgst() ] = nameFrequencies.get(entry.Wgst(), 0) + 1

        #for each entry, calculate the name 
        for entry in tqdm( selection, desc=status ):
            
            qcStatus = entry.QC()
                
            if qcStatus and entry.Key() in namedEntries:
                continue
                
            if existingNames.HasName( entry.Key() ):
                existingNames.DropName( entry.Key() )
            
            #calculate the distance between the unnamed sample and all the named samples
            dists = [ GetDistance( entry.Calls(), self._entryBase.GetEntry( e ).Calls() ) \
                         for e in namedEntries ]
            
            #calculate the name of the entry
            existingNames = CalcName(namedEntries, existingNames, entry.Key(), dists, self.thresholds, qcStatus)
            
            #keep track of the data
            if existingNames.HasResolvedName( entry.Key() ):
                namedEntries.append( entry.Key() )
                nameFrequencies[ existingNames.GetStrName( entry.Key() ) ] = nameFrequencies.get( existingNames.GetStrName( entry.Key() ), 0) + 1


        if iteration == 0:
            for entry in namedEntries:
                self._entryBase.GetEntry( entry ).Wgst( existingNames.GetStrName( entry ) )

        for entry in namedEntries:
            dbEntry = self._entryBase.GetEntry( entry )

            oldName = dbEntry.Wgst()

            newName = existingNames.GetStrName( dbEntry.Key() )

            if oldName == newName:
                continue
            
            if iteration > 0 and oldName and nameFrequencies[oldName] > 1:
                ob = bool( dbEntry.Outbreak() )
                self._wgstHistory.TrackChanges( oldName, newName, iteration, ob = ob, seed = entry in self._startingSet )

            dbEntry.Wgst( newName )
        
        return( existingNames )

def Main(args):

    calc = Calculator( args )
    calc.DoValidation()