###########################################################
# WGST SNSP Plugin
# Author: Milan Patel
# Last Modified: 08/16/2017
# Version: 0.4.1
# Contact: mpatel5@cdc.gov
###########################################################

import os
import sys
import random
import csv
import re
import math
import json
from collections import Counter, defaultdict
from datetime import datetime
import argparse

import numpy as np
from tqdm import *

# Argument parser
def ParseAndValidateCommandLine():

    def Usage():
        usagetxt = \
        """
        You can use this to calculate SNSP to monoclonal
        outbreaks. It will only do so for clinical isolates

        Required:

        [--d]:  Your database file, see example for necessary
                headers.

        [--o]:  Output directory, all the necessary folders
                will be made inside.

        Optional:

        [--p]:  Include polyclonal outbreaks in the calculation.
        """
        return usagetxt

    parser = argparse.ArgumentParser( description = 'Cluster Search tool, requires '
        'query file and database', usage=Usage() )

    parser.add_argument( '--d', '--database', help='Database file', type=str,
        required=True )

    parser.add_argument( '--o', '--outputdir', help='Output directory', type=str,
        required=True, default='./results')

    parser.add_argument( '--p', '--polyclonal', help='Include polyclonal '
        'outbreaks?', type=bool, default=False )

    args = parser.parse_args()

    args.d = os.path.abspath( args.d )
    args.o = os.path.abspath( args.o )

    if not os.path.isfile( args.d ) or not args.d.endswith('.csv'):
        raise IOError('Provide a valid csv formatted db file')

    if not os.path.isdir( args.o ):
        os.mkdir( args.o )

    OUT_DIR = os.path.join( args.o, 'snsp_results' )
    if not os.path.isdir( OUT_DIR ):
        os.mkdir( OUT_DIR )

    args.o = OUT_DIR

    return args

def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))

def unweighted_avg_std(values):
    return (np.average(values), np.std(values))


class Database(object):

    def __init__( self, path, includePoly ):
        
        self._entries = {}
        self._dataPath = path
        self._wgstList = defaultdict(list)
        self.LoadData()
        self.QC( polyclonal = includePoly )
    
    def AddEntry(self, key, mlst, outbreak, srcType, wgst):
        self._entries[ key ] = DbaseEntry( key, mlst, outbreak, srcType, wgst )

    def GetEntry(self, key):
        return self._entries[ key ]

    def GetEntries(self):
        return self._entries

    def LoadData(self):

        refHeader = ('key', 'mlst_st', 'outbreak', 'sourcetype', 'wgst')
        obQualifiers = ('NOT', 'DUP')

        with open( self._dataPath, 'r', encoding='cp1252') as f:
            reader = csv.reader(f)
            header = tuple( map(str.lower, next(reader) ) )

            for head in refHeader:
                if head not in header:
                    raise ValueError('Missing field: {}'.format( head ) )

            for row in tqdm( reader, desc = 'Loading data'):

                fldsAndData = { h:f for h,f in zip( header, row ) }
                
                # SOURCE-TYPE
                if 'human' not in fldsAndData['sourcetype'].lower():
                    continue
                else:
                    srcType = 'human'

                # Check WGST presence
                if len( fldsAndData['wgst'] ) > 1:
                    wgst = list(map( int, fldsAndData['wgst'].split('.') ))

                else:
                    continue

                # BN Key:
                key = fldsAndData['key']

                # 7-GENE MLST
                mlst = fldsAndData['mlst_st'].lower().replace( 
                    'public', '' ) if len( fldsAndData['mlst_st'] ) > 1 else None

                print(mlst)

                # OUTBREAK
                if len( fldsAndData['outbreak'] ) > 1:
                    if any( s in fldsAndData['outbreak'].upper() \
                        for s in obQualifiers ):
                        outbreak = None
                    else:
                        outbreak = fldsAndData['outbreak']
                else:
                    outbreak = None

                self.AddEntry( key, mlst, outbreak, srcType, wgst)
                
        self._numThresholds = len( self.GetEntry( key ).Wgst() )

    def QC(self, polyclonal = False):

        # Initializes a dictionary of outbreaks
        monoclonal = defaultdict(list)

        # Let's find the unique outbreak and fit them in
        for entry in self.GetEntries().values():

            if entry.Outbreak():
                monoclonal[ entry.Outbreak() ].append( entry )

        # We are going to remove outbreaks that don't meet our
        # definition of a monoclonal outbreak
        outbreaksToRemove = set()

        for outbreak in monoclonal:
            # True outbreaks have at least three clinical cases
            if len( monoclonal[ outbreak ] ) < 3:
                outbreaksToRemove.add( outbreak )
                continue

            # Filter polyclonal outbreaks based on pubMLST
            if not polyclonal:

                mlstValues = set()

                for isolate in monoclonal[ outbreak ]:
                    mlst = isolate.Mlst()

                    if mlst == None:
                        continue
                    else: mlstValues.add( mlst )

                if len( mlstValues ) > 1:
                    outbreaksToRemove.add( outbreak )

        for outbreak in outbreaksToRemove:
            for entry in monoclonal[ outbreak ]:
                del self._entries[ entry.Key() ]

            del monoclonal[ outbreak ]

        self._monoclonal = monoclonal
        
        for entry in self.GetEntries().values():
            self._wgstList[ self.StringFromWgst( entry.Wgst() ) ].append( entry )

        for key, value in monoclonal.items():
            print(key)
            for e in value:
                print( repr(e) )

    def GetMonoclonal(self):
        return self._monoclonal

    def GetThresholds(self):
        return self._numThresholds

    def GetWgstList(self):
        return self._wgstList

    def StringFromWgst(self, value):
        return '.'.join (map (str, value) )

    def GetEntryCount(self):
        return len( self._entries )

class DbaseEntry(object):

    def __init__(self, key, mlst, outbreak, srcType, wgst):

        self._key = key
        self._mlst = mlst
        self._outbreak = outbreak
        self._srcType = srcType
        self._wgst = wgst

    def Key(self):
        return self._key

    def Mlst(self):
        return self._mlst

    def Outbreak(self):
        return self._outbreak

    def SrcType(self):
        return self._srcType

    def Wgst(self, level=None):
        if isinstance(level, int):
            return '.'.join(map(str, self._wgst[:level]))
        else:
            return self._wgst

    def __str__(self):
        return self._key

    def __repr__(self):
        return '{0} | {1} | {2}'.format( self._key, self._outbreak, self._srcType )

def CalcSnSp(database):

    levels = database.GetThresholds()
    
    outbreakInfo = defaultdict(dict)
    TOTAL_ISOLATES = database.GetEntryCount()

    monoclonal = database.GetMonoclonal()
    wgstList = database.GetWgstList()

    # For each threshold, let's calculate SNSP at that level
    for level in trange( levels, desc='Calculating' ):
        
        for outbreak in monoclonal:
                        
            FPCount = 0.0
            TPCount = 0.0
            FNCount = 0.0

            modeList = Counter( [ isolate.Wgst(level + 1) for isolate in monoclonal[outbreak] ] )
            anchor = ''

            if len(modeList) == 1:
                for isolate in monoclonal[outbreak]:
                    TPCount += 1.0
                    anchor = next (iter (modeList.keys()) )

            else:
                anchor = max( modeList, key=lambda key: modeList[key] )

                for isolate in monoclonal[outbreak]:
                    if isolate.Wgst(level + 1) == anchor:
                        TPCount += 1.0

                    else:
                        FNCount += 1.0
            
            anchor = anchor.split('.')
            for isoList in wgstList:
                chain = isoList.split('.')[:level + 1]
                
                assert len(anchor) == len(chain)

                if chain == anchor:
                    for key in wgstList[isoList]:
                        if key.Outbreak() != outbreak:
                            FPCount += 1.0

            TNCount = TOTAL_ISOLATES - TPCount - FNCount - FPCount


            outbreakInfo[level + 1][outbreak] = { 'SN': (TPCount/(TPCount + FNCount)), 'SP': ((TNCount)/(TNCount + FPCount)), 'Count': len( monoclonal[outbreak] ) }

    return outbreakInfo


def Main(args):

    dbase = Database( args.d, args.p )
    obInfo = CalcSnSp(dbase)

    # The SNSP of each outbreak contributes different weights
    # to the overall SNSP. Larger outbreaks influence the SNSP
    # more than smaller ones
    weightedSnSp = []

    # Unweighted as well just to be thorough. Usually this one is
    # higher than the weighted SNSP
    unWeightedSnSp = []

    for level in obInfo:
        sn = []
        sp = []
        wght = []

        for outbreak in obInfo[ level ]:

            sn.append( obInfo[level][outbreak]['SN'] )
            sp.append( obInfo[level][outbreak]['SP'] )
            wght.append( obInfo[level][outbreak]['Count'] )


        wghtedSN, snStd = weighted_avg_and_std( sn, wght )
        wghtedSP, spStd = weighted_avg_and_std( sp, wght )

        unWghtedSN, uSnStd = unweighted_avg_std( sn )
        unWghtedSP, uSpStd = unweighted_avg_std( sp )


        weightedSnSp.append( [ '{:03}'.format( level ), '{:08f}'.format( wghtedSN ), '{:08f}'.format( snStd ), \
            '{:08f}'.format( wghtedSP ), '{:08f}'.format( spStd ) ] )
        
        unWeightedSnSp.append( [ '{:03}'.format( level ), '{:08f}'.format( unWghtedSN ), '{:08f}'.format( uSnStd ), \
            '{:08f}'.format( unWghtedSP ), '{:08f}'.format( uSpStd ) ] )

    OUT_NAME = 'snsp_results_{}.csv'.format( 
            datetime.now().strftime("%m-%d-%y@%H-%M-%S") )

    OUT_PATH = os.path.join( args.o, OUT_NAME )

    with open( OUT_PATH, 'w', newline='') as f:
        writer = csv.writer(f)

        writer.writerow([ None, 'Weighted', None, 'Weighted' ])
        writer.writerow([ None, 'Sensitivity', None, 'Specificity', None ])
        writer.writerows( weightedSnSp )
        writer.writerow([])
        writer.writerow([ None, 'Unweighted', None, 'Unweighted' ])
        writer.writerow([ None, 'Sensitivity', None, 'Specificity', None ])
        writer.writerows( unWeightedSnSp )

    OUT_NAME = 'snsp_dump_{}.json'.format( 
            datetime.now().strftime("%m-%d-%y@%H-%M-%S") )

    OUT_PATH = os.path.join( args.o, OUT_NAME)

    with open(OUT_PATH, 'w') as f:

        json.dump(obInfo, f)


if __name__ == '__main__':

    args = ParseAndValidateCommandLine()
    Main( args )