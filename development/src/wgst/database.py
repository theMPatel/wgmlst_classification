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
import re
import json
import numpy as np
import traceback

from tqdm import *
from datetime import datetime
from .environment import Log

__all__ = ['Database', 'DatabaseEntry']

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
                    REUPLOADDATE_4: '%Y-%m-%d'
                    }

REQ_HEADERS = ( 'key', 'outbreak', 'sourcetype', 'uploaddate' )

#=================== DATABASE CLASS =====================================#
class Database(object):

    # Stores the views, you can access these and any custom ones
    # using the class method
    _SCHEMES = {}

    def __init__(self, organism, scheme, metadata, outdir, **kwargs): 
        self._entries = {}
        self._organism = organism
        self._metadata = metadata
        self._scheme = scheme
        self._outdir = outdir
        self._train = []
        self._test = []
        self.SetSelf( kwargs )

    def AddEntry(self, kwargs):
        entry = DatabaseEntry(kwargs)
        self._entries[ entry._key ] = entry
        return entry

    def GetEntry(self, key):
        return self._entries.get( key, None )

    def GetEntries(self):
        return self._entries.values()

    def LoadActiveLoci(self):
        # Get the views

        if hasattr(self, '_schemespath') and self._schemespath \
            is not None:

            # Check to see if it's a real file
            if not os.path.exists(self._schemespath) or not \
                self._schemespath.endswith('.json'):

                raise RuntimeError('Provided views file is not:'
                    ' a real file or not a json file.')

            # Open it
            with open( self._schemespath, 'r') as f:
                try:
                    data = json.load(f)

                except:
                    Log('Provided views file produced'
                        ' and error on reading:\n{}'.format(
                            traceback.format_exc() ))
                    Log('Using all loci')
                    self._scheme = 'WHOLE'
                    return

            # Save it for use
            for view in data:
                Database._SCHEMES[ view ] = set( data[ view ] )

            if self._scheme not in Database._SCHEMES:
                raise RuntimeError('Loci not delineated in'
                    ' json file for specified scheme: {}'.format(
                        self._scheme ))
        else:
            Log( 'Using all loci' )
            self._scheme = 'WHOLE'

    def ExtractData( self, filetype, header, reader ):

        for head in REQ_HEADERS:
            if head not in header:
                raise ValueError('Missing field: {}'.format( head ) )

        if not hasattr(self, '_locisize'):
            raise RuntimeError('Must know the amount of loci to expect')
        
        providedLoci = set()
        for fld in header:
            if fld.startswith(self._organism):
                providedLoci.add( fld )

        if self._locisize > 0 and len( providedLoci ) < self._locisize:
            raise ValueError("Only found {} loci, "
                "make sure you're not missing anything".format( 
                len(providedLoci) ) )

        if self._scheme == 'WHOLE':
            activeLoci = providedLoci
        else:
            activeLoci = Database.GetView( self._scheme )

        for row in tqdm( reader, desc= 'Loading metadata and allele calls' ):
            
            fldsAndData = {}
            alleleCalls = []

            if row[0] == 'key' or self.GetEntry( row[0] ):
                continue

            for h, fld in zip( header, row ):
                # Match the upload date
                if h == 'uploaddate':

                    # Make it a proper variable name
                    h = '_'+h

                    for matchObj in REUPLOADDATE.keys():
                        matched = matchObj.match( fld )

                        if matched:
                            fldsAndData[h] = datetime.strptime(
                                fld,
                                REUPLOADDATE[matchObj] ).date()

                            break

                    if matched is None:

                        if REUPLOADDATE_5.match( fld ):
                            fldsAndData[h] = datetime.strptime(
                                fld,
                                "%Y" ).date()
                        else:
                            fldsAndData[h] = datetime.strptime(
                                '1/1/2000',
                                '%m/%d/%Y' ).date()                       

                # Check for the allele calls
                elif h.startswith( self._organism ) and h in activeLoci:
                    if fld == '' or fld == '?':
                        alleleCalls.append( 0 )
                    else:
                        try:
                            alleleCalls.append( int(fld) )
                        except:
                            tqdm.write(row[0])

                else:
                    fldsAndData['_'+h] = fld if fld else None

            alleleCalls = np.asarray( alleleCalls, dtype=int )
            alleleCalls.flags.writeable = False

            assert len(alleleCalls) == self._locisize

            fldsAndData['_allelecalls'] = alleleCalls
            e = self.AddEntry( fldsAndData )

            # if filetype == 'train':
            #     if e not in self._test:
            #         self._train.append(e)

            # elif filetype == 'test':
            #     if e not in self._train:
            #         self._test.append(e)
            # else:
            #     raise RuntimeError('Invalid filetype provided')

    def LoadCalls(self):

        # Let's load them if present
        self.LoadActiveLoci()

        # For each file, lets add the data
        for filetype, fieldsPath in self._metadata.items():
            with open(fieldsPath, 'r', encoding='cp1252') as f:
                reader = csv.reader( f )
                header = tuple( map( str.lower, next( reader ) ) )
                self.ExtractData( filetype, header, reader )

    def QC( self, minPresence ):

        to_remove = set()
        
        poorQuality = 0
        for entry in tqdm( self.GetEntries(), desc='QC'):
            nonZero = len( entry._allelecalls[entry._allelecalls > 0] )
            total = len( entry._allelecalls )
            qcVal = True

            if float( nonZero ) / float( total ) < minPresence:
                qcVal = False
                poorQuality += 1
                to_remove.add(entry._key)

            if not hasattr( entry, '_qc' ):
                setattr( entry, '_qc', qcVal )
            else:
                entry._qc = qcVal

        for rm in to_remove:
            del self._entries[rm]

        print( 'Poor quality entries: {:d}'.format( poorQuality ) )

    def Save( self ):

        def WriteOut( path, header, attrs ):

            with open(path, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(header)

                for entry in tqdm(self._entries.values(), 
                    desc='Saving: {}'.format(path) ):

                    oArray = []
                    oArray.extend( entry.__getattribute__(x) for x in attrs if hasattr(entry, x ) )
                    writer.writerow( oArray )

        DIRS = { x : os.path.join( self._outdir, x) for x in ['database',
        'snsp_data', 'clustering_validation', 'cluster_search' ] }
        
        for dirc in DIRS.values():
            if not os.path.exists( dirc ):
                os.mkdir( dirc )

        REF_ATTRS = [ '_key', '_outbreak', '_sourcetype', '_uploaddate', '_wgst', 
                        '_wgstHistory' ]

        OUT_NAME = datetime.now().strftime("%m-%d-%y@%H-%M-%S")
        for directory, path in DIRS.items():
            if not os.path.isdir( directory ):
                os.mkdir( directory )

            path_out = os.path.join(path, directory+'_{}.csv'.format(
                OUT_NAME))

            header = [ x[1:].upper() for x in REF_ATTRS ]

            if directory != 'clustering_validation':
                WriteOut( path_out, header, REF_ATTRS )
            else:
                newHeader = header+['ALLELECALLS']
                newREFS = REF_ATTRS+['_allelecalls']
                WriteOut( path_out, newHeader, newREFS )

    def SetSelf( self, kwargs ):
        for attr, val in kwargs.items():
            if attr.startswith('_'):
                setattr( self, attr, val )
            else:
                setattr( self, '_'+attr, val )

    def TotalIsolates(self):
        return len( self._entries )
    
    @staticmethod
    def GetView(view):
        viewValues = Database._SCHEMES.get( view, None )

        if viewValues:

            return set(map(str.lower, viewValues) )
        else:
            raise ValueError('Invalid View Mode')

class DatabaseEntry(object):

    def __init__(self, kwargs):

        self.SetSelf(kwargs)
        self._wgst = None
        self._wgstHistory = []

    def SetSelf( self, kwargs ):
        for attr, val in kwargs.items():
            if attr.startswith('_'):
                setattr( self, attr, val )
            else:
                setattr( self, '_'+attr, val )

if __name__ == '__main__':

    Dbase = Database(organism = 'saen', scheme='core', metadata={'train':'../allele_calls/salm_allele_calls_all.csv'}, outdir='../salm_results/', locisize=1748)

    Dbase.LoadCalls()
    Dbase.QC( 0.30 )
    Dbase.Save()

