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

# ======================== NAME EVENT OBJECT ================================#
class NameEvent(object):
    def __init__(self, oldname, newname, key, info):
        self._oldName = oldname
        self._newName = newname
        self._key = key
        self._info = info

# ========================= NAME HISTORY OBJECT =============================#
class NameHistory(object):
    def __init__(self):
        self._events = []
        pass

    def AddEvent(self, oldName, newName, key, info):
        self._events.append(NameEvent( oldName, newName, key, info ))

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
            addingSet = [ addingSet[i:i+150] for i in range(0, len(addingSet), 150) ]

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

#============================ NAMES OBJECT ==================================#
class Names(object):
    def __init__(self):
        #{ 'CDC_12547' : ['1.1.1.1.1.1'] }
        self._names = defaultdict(list)
        self._history = defaultdict(NameHistory)
        self._allnames = defaultdict(list)
        pass
        
    def _NameToStr(self, parts):
        return '.'.join(str(p) for p in parts)

    def _NameFromStr(self, name):
        return [int(p) for p in name.split('.')]
    
    def Save(self, oFile):
        bstring = json.dumps({
            'names': {k: self._NameToStr(v) for k, v in self._names.items()}, 
            'allnames': list(self._allnames.keys()), 
            })

        oFile.write(bstring)
        
    def Load(self, iFile):
        dataStr = iFile.read()
        data = json.loads(dataStr)

        self._allnames = defaultdict(list)
        for k in data['allnames']:
            self._allnames[k] = self._NameFromStr(k)
        self._names = defaultdict(list)
        for k, v in data['names'].items():
            self._names[k] = self._NameFromStr(v)

    @property
    def History(self):
        return self._history

    def HasName(self, key):
        return key in self._names
        
    def HasResolvedName(self, key):
        return self.HasName(key) and -1 not in self.GetName(key)

    def DropName(self, key):
        if not self.HasName(key): return
        assert(not self.HasResolvedName(key))
        del self._names[key]

    def GetName(self, key):
        return self._names.get(key, [])

    def GetStrName(self, key):
        return '.'.join(str(n) for n in self.GetName(key))

    def GetPart(self, key, i):
        if i<0: return 1
        else:
            if len(self.GetName(key))<=i:
                kkkkk=0
            return self.GetName(key)[i]

    def SetName(self, key, names):
        self._names[key] = names

    def AddToName(self, key, part):
        self._names[key].append(part)
        
    def MakeUndefined(self, key, length):
        self._names[key] = self._names[key] + [-1] * (length-len(self._names[key]))

    def GetNextIdx(self, parts):
        idx = [ p[len(parts)] for k, p in self._allnames.items() if p[:len(parts)]==parts and len(p)>len(parts)]
        if len(idx)==0:
            return 1
        return max(idx)+1 if len(idx)>0 else 1

    def IsSame(self, other):
        nameViolations = []
        nameMap = {}

        for key, name in self._names.items():
            #the real other name of the sample
            otherName = other.GetName(key)

            #the predicted other name of the sample (if any)
            strName = '.'.join(str(n) for n in name)
            strOtherName = '.'.join(str(n) for n in otherName)
            mappedName = nameMap.get(strName, None)
            if mappedName is None: nameMap[strName] = otherName
            else:
                if mappedName!=otherName:
                    #if predicted name and other name do not match, we have a problem
                    nameViolations.append([strName, strOtherName])

        return len(nameViolations) == 0, nameViolations

    def Merge(self, clusterNames, key, info):

        def NameToStr(parts):
            return '.'.join(str(p) for p in parts)

        newName = clusterNames[0]

        nextIdx  = {}
        cnt = 0
        changes = {}
        for clusterName in clusterNames:
            for key, name in self._names.items():

                if len(name)==len(clusterName):
                    if name==clusterName:
                        myOldName = NameToStr(self._names[key])
                        self._names[key][-1] = newName[-1]
                        myNewName = NameToStr(self._names[key])
                        if myOldName != NameToStr(newName):
                            changes[myOldName] = myNewName

                elif name[:len(clusterName)] == clusterName:
                    if clusterName==newName:
                        tail = name[len(clusterName)]
                        nextIdx[NameToStr(name[:len(newName)+1])]=tail

                    myOldName = NameToStr(self._names[key][:len(clusterName)])

                    idx = nextIdx.get(NameToStr(name[:len(newName)+1]), cnt+1)
                    if idx==tail and clusterName!=newName: idx+=1

                    cnt = max(cnt, idx)
                    nextIdx[NameToStr(name[:len(newName)+1])] = idx

                    self._names[key][len(newName)-1] = newName[-1]
                    self._names[key][len(newName)] = idx

                    myNewName = NameToStr(self._names[key][:len(clusterName)])
                    if myOldName != NameToStr(newName):
                        changes[myOldName] = myNewName

        for oN, nN in changes.items():
            if oN!=nN:
                self._history[oN].AddEvent(oN, nN, key, info)

        return newName[-1]
        
    def Finalize(self, key):
        self._allnames[self.GetStrName(key)] = self.GetName(key) 

#=========================== NAMING FUNCTION ================================#
def CalcName(namedEntries, existingNames, unNamedEntry, distances, 
    thresholds, qcStatus):
    """
    :param namedEntries: keys you want to use to identify the samples, can be anything but should be
                    in the same order as the list of distances
    :param existingNames: existing names and therefore clusters
    :param unNamedEntry: key for the entry to be named
    :param distances: list of distances between the unNamedEntry and the namedEntries
    :param thresholds: list of thresholds
    :return: a Names object holding the hierarchical names of all the samples
    """

    #create the object that will hold the names
    names = existingNames

    #create the initial cluster: all entries lumped together
    nEntries = len(namedEntries)+1

    #make sure the thesholds are sorted as biggest-first
    myThresholds = thresholds
    myThresholds.sort(key=lambda x: -x)

    startCluster = range(len(namedEntries)+1)
    cluster = startCluster


    for partNr, threshold in enumerate(myThresholds):

        #create the clusters at this level
        subClusters = defaultdict(list)
        for j in cluster:
            if j==len(namedEntries): continue
            assert(names.GetPart(namedEntries[j], partNr)>0)
            subClusters[names.GetPart(namedEntries[j], partNr)].append(j)

        #find the clusters that are close enough
        closestClusters = set()
        for j, cluster in subClusters.items():
            for k in cluster:
                if distances[k]<=threshold:
                    closestClusters.add(j)
        
        #assigned to just one cluster, nothing special happening here
        if len(closestClusters)==1:
            theOne = closestClusters.pop()
            cluster = subClusters[theOne]
            names.AddToName(unNamedEntry, theOne)
        #bad quality, jump out 
        elif not qcStatus:
            names.MakeUndefined(unNamedEntry, len(myThresholds))
            break
        #create new cluster
        elif len(closestClusters)==0:
            cluster = []
            idx = names.GetNextIdx(names.GetName(unNamedEntry))
            names.AddToName(unNamedEntry, idx)
        #assigned to multiple clusters, need to merge these clusters    
        else:
            #assigned to multiple clusters, need to merge these clusters
            cluster = []
            clusterSizes = []
            for c in closestClusters:
                cluster.extend(subClusters[c])
                clusterSizes.append([c, len(subClusters[c])])

            clusterSizes.sort(key=lambda x: -x[1])
            prefix = names.GetName(unNamedEntry)
            newName = names.Merge([list(prefix) + [c[0]] for c in clusterSizes], unNamedEntry, None)
            names.AddToName(unNamedEntry, newName)
    
    if qcStatus:
        names.Finalize( unNamedEntry )
    
    return names

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