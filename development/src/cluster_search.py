###########################################################
# Cluster Search Plugin
# Author: Milan Patel
# Last Modified: 08/16/2017
# Version: 0.2
# Contact: mpatel5@cdc.gov
###########################################################

import os
import sys
import re
import csv
import argparse

from datetime import datetime
from operator import itemgetter
from collections import defaultdict
from tqdm import *

# These are used to extract date, if you're seeing issues with
# dates not being extracted properly, please make sure they're
# at least in the below formats
REUPLOADDATE_1 = re.compile( r'([0-9]{1,2})/([0-9]{1,2})/([0-9]{4})' )
REUPLOADDATE_2 = re.compile( r'([0-9]{4})/([0-9]{1,2})/([0-9]{1,2})' )
REUPLOADDATE_3 = re.compile( r'([0-9]{1,2})-([0-9]{1,2})-([0-9]{4})' )
REUPLOADDATE_4 = re.compile( r'([0-9]{4})-([0-9]{1,2})-([0-9]{1,2})' )
REUPLOADDATE_5 = re.compile( r'([0-9]{4})' )

REUPLOADDATE = {
                    REUPLOADDATE_1: '%m/%d/%Y',
                    REUPLOADDATE_2: '%Y/%m/%d',
                    REUPLOADDATE_3: '%m-%d-%Y',
                    REUPLOADDATE_4: '%Y-%m-%d',
                    REUPLOADDATE_5: '%Y'
                    }

def ParseAndValidateCommandLine():

    def Usage():
        usagetxt = \
        """
        You can use this to calculate clusters outside
        of BioNumerics. Provide a query file and a historical
        database file. 

        Required:

        [--q]:  Query file, see example for necessary headers

        [--d]:  Historical database, really only needed if you
                want to include historical isolates.If you want 
                to just search a singular file: point your query 
                and database to the same file

        [--o]:  Output directory, algorithm will make the
                necessary folders inside

        [--s]:  Provide int for the level of similarity.
                This means you should provide the number
                of digits you want isolates to share
                in order to be considered a cluster

        Optional:

        [--h]:  Include historical isolates into the clusters?
                Boolean value

        [--c]:  Include nonclinical isolates
                Boolean value
        """
        return usagetxt

    parser = argparse.ArgumentParser( description = 'Cluster Search tool, requires '
        'query file and database', usage=Usage() )

    parser.add_argument( '--q', '--query', help='Query file', type=str,
        required=True )

    parser.add_argument( '--d', '--database', help='Database file', type=str,
        required=True )

    parser.add_argument( '--o', '--outputdir', help='Output directory', type=str,
        required=True )

    parser.add_argument( '--s', '--similarity', help='Higher is more similar', type=int,
        required=True )

    parser.add_argument( '--h', '--historical', help='Include historical '
        'isolates?', type=bool, default=False )

    parser.add_argument( '--c', '--nonclinical', help='Inlcude nonclinical '
        'isolates?', type=bool, default=False )

    args = parser.parse_args()

    if not os.path.isfile( args.q ) or not args.q.endswith('.csv'):
        raise IOError('Provide a valid csv formatted query file')

    if not os.path.isfile( args.d ) or not args.d.endswith('.csv'):
        raise IOError('Provide a valid csv formatted db file')

    if not os.path.isdir( args.o ):
        os.mkdir( args.o )

    args.q = os.path.abspath( args.q )
    args.d = os.path.abspath( args.d )
    args.o = os.path.abspath( args.o )

    return args


class ClusterSearch( object ):

    def __init__( self, query, db, similarity, resultsdir, inlcudeHistorical = False,
        includeNonClinical = False ):

        self._query = query
        self._database = db
        self._similarity = similarity
        self._resultsDir = resultsdir
        self._includeHistorical = inlcudeHistorical
        self._includeNonClinical = includeNonClinical

    def CalcClusters(self):

        # Get the list
        isoList = self.CreateList()
        percSimil = self._similarity
        self._clustersToCheck = {}
        start = 0

        # Do the cluster search
        while (start < len(isoList)-1):

            toCheck = isoList[start][1][:percSimil]
            i = start

            if i == len(isoList)-1:
                break

            while (True):

                if i == len(isoList)-1:
                    start = i+1
                    break

                elif toCheck == isoList[i+1][1][:percSimil]:
                    try:

                        string = '.'.join( map( str, toCheck ) )
                        self._clustersToCheck[string].AddIso(isoList[i+1][0])
                        i += 1

                    except(KeyError):

                        string = '.'.join( map( str, toCheck ) )
                        self._clustersToCheck[string] = ClusterParent(isoList[start][0])
                        self._clustersToCheck[string].AddIso(isoList[i+1][0])
                        i += 1

                else:
                    start = i + 1
                    break

        # Make sure there are at least three
        toRemove = [ cluster for cluster in self._clustersToCheck if 
            self._clustersToCheck[cluster].ClusterCount() < 3 ]

        # Remove the two-sies
        for remove in toRemove:
            del self._clustersToCheck[remove]

        # Find the historical isolates
        if self._includeHistorical:

            wgstList = defaultdict(list)

            for entry in self._database:
                wgstList[ entry.WgstPartial(self._similarity) ].append( 
                    entry.Key() )

            for cluster in self._clustersToCheck:
                
                if wgstList.get( cluster, None ) is None:
                    continue

                self._clustersToCheck[cluster].UpdateClusterList(
                    wgstList[cluster] )

    def CreateList(self):

        isoList = []

        for entry in self._query.values():

            if self._includeNonClinical:
                isoList.append( [ entry.Key(), entry.NumericalWgst() ] )

            elif 'human' in entry.SrcType():
                isoList.append( [ entry.Key(), entry.NumericalWgst() ] )

        isoList.sort( key = itemgetter( 1 ) )

        return isoList

    def Save(self):

        if len( self._clustersToCheck ) == 0:
            print('No clusters found!')
            sys.exit(0)

        OUT_DIR = os.path.join( self._resultsDir, 'cls_search_results' )

        if not os.path.isdir( OUT_DIR ):
            os.mkdir( OUT_DIR )

        OUT_NAME = 'found_clusters_{}.csv'.format( 
            datetime.now().strftime("%m-%d-%y@%H-%M-%S") )
        
        OUT_PATH = os.path.join( OUT_DIR, OUT_NAME )

        with open( OUT_PATH, 'w', newline='' ) as f:
            writer = csv.writer( f )

            for i, cluster in enumerate( self._clustersToCheck, 1 ):
                writer.writerow( ['Cluster {}'.format(i) ] )

                for iso in self._clustersToCheck[cluster].ClusterList():
                    
                    entry = self._query.get( iso, None )

                    if entry is None:
                        entry = self._database.get( iso, None )

                        if entry is None:
                            raise KeyError('Something went really wrong, should '
                                'have found this key: {}'.format( iso ) )

                    writer.writerow( [entry.Key(), entry.SrcType(),
                        entry.UploadDate(), entry.Wgst() ] )

class ClusterParent(object):

    def __init__(self, key):
        self._parent = key
        self._simIsoList = set()
        self.AddIso( self._parent )

    def AddIso(self, value):
        self._simIsoList.add( value )

    def ClusterCount(self):
        return len( self._simIsoList )

    def ClusterList(self):
        return self._simIsoList

    def UpdateClusterList(self, value):
        self._simIsoList.update( value )

class Database( object ):

    def __init__( self, searchFile, dbPath):
        self._queryPath = searchFile
        self._dbPath = dbPath
        self._query = {}
        self._database = {}
        self.Load()
    
    def GetEntries(self):
        return self._database

    def GetQueries(self):
        return self._query
    
    def Load(self):
        refHeader = ('key', 'sourcetype', 'uploaddate', 'wgst')


        with open( self._queryPath, 'r' ) as f:

            reader = csv.reader( f )
            header = tuple( map( str.lower, next( reader ) ) )

            for head in refHeader:
                if head not in header:
                    raise ValueError('Missing field: {}'.format( head ) )

            for row in tqdm( reader, desc='Loading search file'):

                fldsAndData = { h:f for h,f in zip( header, row ) }

                if len( fldsAndData['key'] ) < 1 or \
                len( fldsAndData['wgst'] ) < 1:
                    continue

                key = fldsAndData['key']

                uploadDate = None
                for matchObj in REUPLOADDATE.keys():
                    if matchObj.match( fldsAndData['uploaddate'] ):
                        uploadDate = datetime.strptime( 
                            fldsAndData['uploaddate'] , REUPLOADDATE[matchObj] ).date()
                        break

                if uploadDate is None:
                    uploadDate = datetime.strptime( '1/1/2000', "%m/%d/%Y" ).date()

                srcType = fldsAndData['sourcetype']

                wgst = fldsAndData['wgst']

                self._query[ key ] = DbaseEntry( key, uploadDate, 
                    srcType, wgst )

        with open( self._dbPath, 'r' ) as f:
            reader = csv.reader( f )
            header = tuple( map( str.lower, next( reader ) ) )

            for head in refHeader:
                if head not in header:
                    raise ValueError('Missing field: {}'.format( head ) )

            for row in tqdm( reader, desc='Loading historical'):

                fldsAndData = { h:f for h,f in zip( header, row ) }

                if len( fldsAndData['key'] ) < 1 or \
                len( fldsAndData['wgst'] ) < 1:
                    continue

                key = fldsAndData['key']

                uploadDate = None
                for matchObj in REUPLOADDATE.keys():
                    if matchObj.match( fldsAndData['uploaddate'] ):
                        uploadDate = datetime.strptime( 
                            fldsAndData['uploaddate'] , REUPLOADDATE[matchObj] ).date()
                        break

                if uploadDate is None:
                    uploadDate = datetime.strptime( '1/1/2000', "%m/%d/%Y" ).date()

                srcType = fldsAndData['sourcetype']

                wgst = fldsAndData['wgst']

                self._database[ key ] = DbaseEntry( key, uploadDate,
                    srcType, wgst )

class DbaseEntry( object ):

    def __init__( self, key, uploadDate, srcType, wgst ):

        self._key = key
        self._srcType = srcType
        self._uploadDate = uploadDate
        self._wgst = wgst
        self._numericalWgst = list( map( int, self._wgst.split('.') ) )

    def Key(self):
        return self._key

    def NumericalWgst(self):
        return self._numericalWgst
    
    def SrcType(self):
        return self._srcType

    def UploadDate(self):
        return self._uploadDate

    def Wgst(self):
        return self._wgst

    def WgstPartial(self, value):
        return '.'.join( map( str,
            self._numericalWgst[:value] ) )
        

    def __str__(self):
        return self._key

    def __repr__(self):
        return '{0} | {1} | {2}'.format( self._key, self._srcType, self._wgst )


def Main( args ):

    db = Database( args.q, args.d )
    cc = ClusterSearch( 
                        db.GetQueries(), 
                        db.GetEntries(),
                        args.s,
                        args.o,
                        args.h,
                        args.c)
    cc.CalcClusters()
    cc.Save()


if __name__ == '__main__':

    args = ParseAndValidateCommandLine()
    Main( args )
