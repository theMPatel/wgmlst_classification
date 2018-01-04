###############################################################################
#
# DATABASE CLASS FOR PAIRING WITH NOMENCLATURE SCRIPT. 
# AUTHOR: MILAN PATEL
# CONTACT: mpatel5@cdc.gov
#
###############################################################################

import os
import sys
import csv
from datetime import datetime
import re
import json
import numpy as np
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

#=================== DATABASE CLASS =====================================#
class DB(object):

    # Stores the views, you can access these and any custom ones
    # using the class method
    VIEWS = {}

    def __init__(self, scheme, fieldsPath, viewsPath, outdir, adding_file = None ): 
        self._entries = {}
        self._fieldsPath = fieldsPath
        self._viewsPath = viewsPath
        self._view = scheme
        self._outdir = outdir

        if adding_file is not None:
            self._addingPath = adding_file
            self._addExternal = True
        else:
            self._addExternal = False

        self.LoadViews()
        self.LoadCalls()

    def AddEntry(self, key, mlst, outbreak, srcType, uploadDate, alleleCalls):
        self._entries[ key ] = DbaseEntry( key, mlst, outbreak, srcType, uploadDate, alleleCalls )

    def GetEntry(self, key):
        return self._entries[ key ]

    def GetEntries(self):
        return self._entries

    def LoadViews(self):
        with open( self._viewsPath, 'r') as f:
            data = json.load(f)

        for view in data:
            DB.VIEWS[ view ] = tuple( data[ view ] )

    def ExtractData( self, header, reader, **flags ):

        refHeader = ( 'key', 'mlst_st', 'outbreak', 'sourcetype', 'uploaddate' )

        for head in refHeader:
            if head not in header:
                raise ValueError('Missing field: {}'.format( head ) )

        if 'externaldata' in flags:
            self._addingSet = []

        loci = 0
        for fld in header:
            if fld.lower().startswith('lmo'):
                loci += 1

        if loci < 1747:
            raise ValueError("Only found {} loci, "
                "make sure we're not missing anything".format( 
                str( loci) ) )

        for row in tqdm( reader, desc= 'Loading metadata and allele calls' ):
            
            fldsAndData = { h:f for h,f in zip( header, row ) }
            
            key = fldsAndData['key']
            
            if len( fldsAndData['outbreak'] ) > 1:
                outbreak = fldsAndData['outbreak']
            else:
                outbreak = None

            if len( fldsAndData['mlst_st'] ) > 1:
                mlst = fldsAndData['mlst_st'].lower().strip()
            else:
                mlst = fldsAndData['mlst_st']

            mlst = None

            if len( fldsAndData['sourcetype'] ) > 1:
                srcType = fldsAndData['sourcetype'].lower().strip()
            else:
                srcType = None

            uploadDate = None
            for matchObj in REUPLOADDATE.keys():
                if matchObj.match( fldsAndData['uploaddate'] ):
                    uploadDate = datetime.strptime( 
                    fldsAndData['uploaddate'] , REUPLOADDATE[matchObj] ).date()

            if uploadDate is None:
                uploadDate = datetime.strptime( '1/1/2000', "%m/%d/%Y" ).date()

            if self._view == 'whole':
                alleleCalls = []
                for i, fld in enumerate( header ):
                    if fld.startswith('lmo'):
                        alleleCalls.append( row[i] )

                alleleCalls = np.asarray( alleleCalls, dtype=int )
                alleleCalls.flags.writeable = False

            else:
                activeLoci = DB.GetView( self._view )
                calls = {}
                for i, fld in enumerate( header ):
                    if fld.startswith('lmo'):
                        calls[fld] = row[i]

                alleleCalls = np.asarray( [ calls[locus.lower()] for locus in activeLoci ], dtype=int )
                alleleCalls.flags.writeable = False

            if 'externaldata' in flags:
                if key not in self._entries:
                    self.AddEntry( key, outbreak, uploadDate, alleleCalls)
                    self._addingSet.append( self._entries[ key ] )
            else:
                self.AddEntry( key, mlst, outbreak, srcType, uploadDate, alleleCalls)

    def LoadCalls(self):

        with open(self._fieldsPath, 'r', encoding='cp1252') as myFile:
            reader = csv.reader( myFile )
            header = tuple( map( str.lower, next( reader ) ) )

            self.ExtractData( header, reader )

        if self._addExternal:

            with open( self._addingPath, 'r', encoding='cp1252' ) as f:
                reader = csv.reader( f )
                header = tuple ( map( str.lower, next( reader ) ) )

                self.ExtractData( header, reader, externaldata=True)

    def QC( self, minPresence ):
        poorQuality = 0
        for entry in tqdm( self.GetEntries().values(), desc='QC'):
            
            nonZero = len( list( filter( lambda x: x != 0, entry.Calls() ) ) )
            total = len( entry.Calls() )

            if float( nonZero ) / float( total ) < minPresence:
                entry.SetQC( False )
                poorQuality += 1
            else:
                entry.SetQC( True )

        print( 'Poor quality entries: {:d}'.format( poorQuality ) )

    def Save( self, clusterSearch = False ):

        OUT_NAME = datetime.now().strftime("%m-%d-%y@%H-%M-%S")
        OUT_DIR = os.path.join( self._outdir, 'snsp_data')

        if not os.path.isdir( OUT_DIR ):
            os.mkdir( OUT_DIR )

        PATH = os.path.join( OUT_DIR , 'database_{}.csv'.format( OUT_NAME ) )

        with open(PATH, 'w', newline='') as f:
            writer = csv.writer( f )

            writer.writerow( ['KEY', 'MLST_ST', 'OUTBREAK', 'SOURCETYPE', 'UPLOADDATE', 'WGST', \
                'WGST History'] )

            for entry in tqdm( self._entries.values(), desc = 'Saving' ):
                oArray = [ entry.Key(), entry.Mlst(), entry.Outbreak(), entry.SrcType(), \
                    entry.UploadDate(), entry.Wgst(), '|'.join( map( str, entry.WgstHistory() ) ) ]

                writer.writerow( oArray )

        # The below code is if you want to save the database to calculate
        # intra/inter cluster error

        OUT_DIR = os.path.join( self._outdir, 'clustering_validation')

        if not os.path.isdir( OUT_DIR ):
            os.mkdir( OUT_DIR )

        PATH = os.path.join( OUT_DIR , 'database_{}.csv'.format( OUT_NAME ) )

        with open( PATH, 'w', newline='') as f:

            writer = csv.writer( f )

            for entry in tqdm( self._entries.values(), desc = 'Saving For Validation' ):
                oArray = [ entry.Key(), entry.Wgst() ]
                oArray.extend( list( entry.Calls() ) )
                writer.writerow( oArray )

        if clusterSearch:
            OUT_DIR = os.path.join( self._outdir, 'cluster_search' )

            if not os.path.isdir( OUT_DIR ):
                os.mkdir( OUT_DIR )

            OUT_NAME = datetime.now().strftime("%m-%d-%y@%H-%M-%S")
            PATH = os.path.join( OUT_DIR , 'cls_search_{}.csv'.format( OUT_NAME ) )

            with open( PATH, 'w', newline='' ) as f:
                writer = csv.writer( f )

                writer.writerow( ['KEY', 'MLST_ST', 'OUTBREAK', 'SOURCETYPE', 'UPLOADDATE', 'WGST', \
                    'WGST History'] )

                for entry in tqdm( self._addingSet, desc='Saving for cluster search' ):
                    oArray = [ entry.Key(), entry.Mlst(), entry.Outbreak(), entry.SrcType(), \
                        entry.UploadDate(), entry.Wgst(), '|'.join( map( str, entry.WgstHistory() ) ) ]
                    
                    writer.writerow(oArray)

    def TotalIsolates(self):
        return len( self._entries )
    
    @classmethod
    def GetView(cls, view):
        viewValues = cls.VIEWS.get( view, None )

        if viewValues:
            return viewValues
        else:
            raise ValueError('Invalid View Mode')

class DbaseEntry(object):

    def __init__(self, key, mlst, outbreak, srcType, uploadDate, alleleCalls):

        self._key = key
        self._mlst = mlst
        self._outbreak = outbreak
        self._srcType = srcType
        self._uploadDate = uploadDate
        self._alleleCalls = alleleCalls
        self._wgst = None
        self._qcStatus = None
        self._wgstHistory = []

    def Key(self):
        return self._key
    
    def Mlst(self):
        return self._mlst

    def Outbreak(self):
        return self._outbreak

    def SetQC(self, value):
        assert isinstance( value, bool)
        self._qcStatus = value

    def SrcType(self):
        return self._srcType

    def QC(self):
        return self._qcStatus

    def UploadDate(self):
        return self._uploadDate

    def Wgst(self, value = None):
        if not value:
            return self._wgst
        else:
            assert isinstance( value, str )
            self._wgst = value
            self._wgstHistory.append( value )

    def WgstHistory(self):
        return self._wgstHistory

    def Calls( self ):
        return self._alleleCalls

    def __str__(self):
        return self._key

    def __repr__(self):
        return '{0} | {1} | {2}'.format( self._key, self._outbreak, self._uploadDate )

if __name__ == '__main__':

    Dbase = DB('core')
    Dbase.QC( 0.95 )
