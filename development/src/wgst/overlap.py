###########################################################
# Distance Matrix calculation script
# Author: Milan Patel
# Last Modified: 08/16/2017
# Version: 0.1
# Contact: mpatel5@cdc.gov
###########################################################

import os
import sys
import json
import csv
import itertools
import numpy as np
from tqdm import *

#======================== DISTANCE FUNCTION ===============================#
def GetDistance(p1, p2):
    # Find the loci that are present in both samples
    common = np.multiply(p1>0, p2>0)

    # Count the total number of loci that are present for both
    nCommon = np.sum(common)

    # For where they both have a presence, see where they are the same 
    nSame = np.sum(p1[common]==p2[common])

    # If they have anything in common, find out how much
    if nCommon:
      return 100.0 * (float(nCommon) - float(nSame)) / float(nCommon)
    
    # This shouldn't happen too often
    else:
      return 100.0

class Database(object):

    def __init__(self, fields_path, minpresence, scheme, viewsPath):
        self._entries = {}
        self._fieldsPath = fields_path
        self._minPres = minpresence
        self._scheme = scheme
        self._viewsPath = viewsPath

        self.LoadViews()
        self.LoadCalls()

    def AddEntry( self, key, profiles):
        self._entries[ key ] = DbaseEntry( key, profiles )

    def GetEntries(self):
        return self._entries

    def GetView(self, view):
        viewValues = self._views.get( view, None )

        if viewValues:
            return viewValues
        else:
            raise ValueError('Invalid View Mode')

    def LoadCalls(self):

        with open( self._fieldsPath, 'r', encoding='cp1252') as f:
            reader = csv.reader( f )
            header = tuple( map(str.lower, next( reader ) ) )

            loci = 0
            for fld in header:
                if fld.lower().startswith('lmo'):
                    loci += 1

            if loci < 1747:
                raise ValueError("Only found {} loci, "
                    "make sure we're not missing anything".format( 
                    str( loci) ) )

            # Load the data from the file
            for row in tqdm( reader, desc='Loading data for matrix calc'):
                
                fldsAndData = { h:f for h,f in zip( header, row ) }

                """
                Make sure not to have empty keys in your files, but this
                should have you covered
                """
                if len( fldsAndData['key'] ) > 0:
                    key = fldsAndData['key']
                else:
                    continue

                if self._scheme == 'whole':
                    alleleCalls = []
                    for i, fld in enumerate( header ):
                        if fld.startswith('lmo'):
                            alleleCalls.append( row[i] )

                    alleleCalls = np.asarray( alleleCalls, dtype=int )
                    alleleCalls.flags.writeable = False

                else:
                    activeLoci = self.GetView( self._scheme )
                    calls = {}
                    for i, fld in enumerate( header ):
                        if fld.startswith('lmo'):
                            calls[fld] = row[i]

                    alleleCalls = np.asarray( [ calls[locus.lower()] for locus in activeLoci ], dtype=int )
                    alleleCalls.flags.writeable = False

                # Check how many allele calls are present
                present = sum( alleleCalls > 0)

                # Calc the % present
                if present / len( alleleCalls ) < self._minPres: continue

                alleleCalls.flags.writeable = False

                self.AddEntry( key, alleleCalls )

        print('Total isolates: {}'.format( len( self._entries ) ) )
    
    def LoadViews(self):
        self._views = {}

        with open( self._viewsPath, 'r') as f:
            data = json.load(f)

        for view in data:
            self._views[ view ] = tuple( data[ view ] )

class DbaseEntry(object):

    def __init__(self, key, alleleCalls):

        self._key = key
        self._alleleCalls = alleleCalls
        self._distances = {}

    def AddDistance(self, other, dist):
        self._distances[ other ] = dist

    def GetDistances(self):
        return self._distances

    def Key(self):
        return self._key

    def Calls( self ):
        return self._alleleCalls


def Main(args):

    # This is where we will put the distance matrix
    DIR_PATH = os.path.join( args['outdir'], 'dist_data' )
    FILE_PATH = os.path.join(DIR_PATH, \
        'distance_matrix_{}_qcd.json').format( args['scheme'] )
    
    if args['recalculate']:
        if os.path.isfile( FILE_PATH ):
            print( 'Deleting old matrix' )
            os.remove( FILE_PATH )

    # If the matrix exists for the passed scheme, let's load it
    if os.path.isfile( FILE_PATH ):
        print( 'Loading distance matrix...' )
        with open( FILE_PATH, 'r') as f:
            distance_data = json.load(f)

        # ET Phone Home
        return distance_data

    elif not os.path.isdir( DIR_PATH ):
        os.mkdir( DIR_PATH )

    # Initialize the database
    dbase = Database(
        args['fields_path'],
        args['minPres'],
        args['scheme'],
        args['views_path']
        )

    # Calculate the matrix
    for entry, other in tqdm( itertools.combinations( 
        dbase.GetEntries().values(), 2), desc = 'Creating Matrix' ):

        dist = GetDistance( entry.Calls(), other.Calls() )
        entry.AddDistance( other.Key(), dist )
        other.AddDistance( entry.Key(), dist )

    # Here we will store the matrix -> { Key: {Other Keys: %Similiarity} }...
    distance_data = {}

    # Build the matrix for shipment
    for entry in tqdm( dbase.GetEntries().values(), desc='Organizing matrix' ):
        distance_data[ entry.Key() ] = entry.GetDistances()

    # Save the matrix for future use
    print( 'Dumping matrix...this might take a while' )
    with open( FILE_PATH, 'w' ) as f:

        json.dump(distance_data, f)

    # I'm coming home
    return distance_data