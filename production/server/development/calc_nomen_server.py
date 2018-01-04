##############################################################################
#
# Nomenclature naming script
# Contact: mpatel5@cdc.gov, hcarleton@cdc.gov, hannes.pouseele@applied-maths.com
# Version: 1.0
#
##############################################################################

# Stdlib imports
import os
import sys
import base64
import gzip
import json
import shutil
import traceback
import numpy as np

from datetime import datetime, date
from collections import defaultdict

# Bionumerics constants
import bns
from wgMLST_Client.wgMLSTSchema import Schema
import xml.etree.cElementTree as ET

Database = bns.Database
Fields = Database.Db.Fields
Characters = bns.Characters

# Params required for this organism
GB_PARAMS = {
    'bnparamkey': 'hierarchical_strain_nomenclature',
    'organism': 'LMO',
    'version': '1.0',
    'scheme': 'core7',
    'thresholds': (4.050815987,2.914238135,2.046511241,
        1.068159867,0.0),
    'minpres': 0.95,
    'namesgiven': 0,
    'instantiated': None,
    'wgst': None,
    'wgmlstname': None,
    'quality': {'Length':[2.8e6, 3.2e6]}
}

# Data directories
DATA_DIR = Database.Db.Info.SourceFilesDir
DATA_DIRS = [ 'tree', 'tree/current', 'allele_calls', 
    'allele_calls/current', 'logs' ]

############################# LOGGER #########################################
class Logger(object):
    """
    Logger class, will be used to log things during the naming process
    """

    current = None

    def __init__(self, logDir):
        path = [
            '{}_nomenclature_srcfiles'.format( 
            GB_PARAMS['organism'] ),
            'logs',
            'nomenclature_log_{}.txt'.format( 
                datetime.now().strftime("%m-%d-%y@%H-%M-%S"))
             ]
        self._fileName = os.path.join( logDir, *path )
        self._logStr = ''
        Logger.current = self

    def _log(self, outStr, depth=0):
        if not isinstance(outStr, basestring):
            
            try:
                self._logStr += '->' + ''.join(['\t']*depth) + str(outStr) + '\n'

            except:
                e = traceback.format_exc()
                self._logStr += '->' + ''.join(['\t']*depth) + e + '\n'

        else:
            self._logStr += '->' + ''.join(['\t']*depth) + outStr + '\n'

    def _save(self):
        with open(self._fileName, 'w') as f:
            f.write( self._logStr )

################################# CLUSTERING LOGIC ###########################
class Tree( object ):
    
    CURRENT = None
    DEPTH = -1

    def __init__(self, depth):
        
        Tree.CURRENT = self
        Tree.DEPTH = depth
        self._tree = Node( 1, 0, None )
        self._names = defaultdict(list)
        self._cdcNames = {}

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

    def CDCName(self, part):
        if not isinstance( part, list ):
            raise AssertionError('Need list to create CDC name!')

        code = GB_PARAMS['organism'] + \
                GB_PARAMS['version'] + \
                ' - ' + \
                Tree.NameToStr( part )

        return code

    def FinalizeName( self, key, name ):
        node = self.Traverse(name)
        assert( node.IsNamed( key ) )
        self._names[ key ] = name

    def FinalizeCDCNames( self ):
        for node in self._tree.GetChildrenNodes():
            for namedNode in node.DFSNamed():
                cdcname = namedNode.NTraverse()
                for i, val in enumerate( cdcname, 1):
                    if val > 0:
                        part = len(cdcname) - i
                        break

                for key in namedNode.GetChildrenKeys():
                    partialName = self.CDCName( self.GetPart( key, part ) )
                    complete = part == Tree.DEPTH
                    yield ( key, partialName, complete )

    def GetName( self, key ):
        return self._names.setdefault(key, [])

    def GetNames(self):
        return self._names

    def GetPart( self, key, level ):
        return self.GetName( key )[ :level ]

    def GetStrName( self, key ):
        return Tree.NameToStr( self.GetName(key) )

    def HasName( self, key ):
        return key in self._names and len(self._names[key]) == Tree.DEPTH

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

############################ ALLELE CALLS ####################################
class AlleleCalls(object):

    current = None

    def __init__(self, dir_path, loci_length):

        self._allele_calls = None
        self._dir_path = dir_path
        self._length = -1
        self._capacity = -1
        self._last_index = -1
        self._keys_to_index = {}
        self._invalid_indices = set()
        self._loci_length = loci_length

        AlleleCalls.current = self

    def invalidate(self, key):
        index = self._keys_to_index[key]
        self._allele_calls[index] = -1
        self._invalid_indices.add(index)

    def load(self):
        meta_data_path = os.path.join(self._dir_path, 'metadata.json')

        # Load the metadata
        if os.path.exists( meta_data_path ):
            with open(meta_data_path, 'r') as f:
                self._metadata = json.load(f)
                self._length = self._metadata['size']
                self._capacity = self._metadata['capacity']
                self._last_index = self._metadata['last']

        # Set the default metadata size
        else:
            self._metadata = { 'size':-1, 'capacity': -1, 'runs': 1, 'last': 0 }

        # Load the key mappings so we know key->index
        key_mapping_path = os.path.join(self._dir_path, 
            'keys_mapping.json')

        # Check to see if it exists
        if os.path.exists( key_mapping_path ):
            with open(key_mapping_path, 'r') as f:
                self._keys_to_index = json.load(f)

            if len(self._keys_to_index) != self._metadata['size']:
                raise RuntimeError('Invalid size for loaded mapping')
        
        # New instance, or something catastropic occurred
        else:

            if self._metadata['size'] != -1 or self._metadata['capacity'] != -1:
                raise RuntimeError( 'Missing keys_mapping file for'
                    'associated metadata' )

        array_path = os.path.join( self._dir_path, 
            'allele_calls_array.memmap' )

        if os.path.exists( array_path ):

            # In case we need to resize
            flush_flag = False

            # If we find a file but have no metadata
            if self._metadata['size'] == -1:
                raise RuntimeError('Missing metadata for '
                    'found memmap file')

            # If we need to resize
            if (self._metadata['size'] / self._metadata['capacity']) > 0.75:
                self._metadata['capacity'] = int(self._metadata['size'] * 1.25)+1
                flush_flag = True
                self._capacity = self._metadata['capacity']


            shape = ( self._metadata['capacity'], self._loci_length )
            
            try:

                self._allele_calls = np.memmap( 
                    array_path, dtype='int32', mode='r+', shape=shape )

                # If we have resized, let's flush to disk (J.I.C.)
                if flush_flag:
                    self._allele_calls.flush()

            except:

                Logger.current._log(str(self._metadata))
                Logger.current._log(array_path)
                Logger.current._log(shape)
                raise

            if self._metadata['runs']%20 == 0:
                self._validate_calls()

            # Load the invalid indices
            invalid_path = os.path.join(self._dir_path, 'invalid_indices.json')

            if os.path.exists(invalid_path):
                with open(invalid_path, 'r') as f:
                    self._invalid_indices.update( json.load(f) )

        else:
            # No file exists
            if self._metadata['size'] != -1 or self._metadata['capacity'] != -1:
                raise RuntimeError('Found array metadata but no array')

            self._allele_calls = np.memmap(
                array_path, dtype='int32', mode='w+', shape=(2000, self._loci_length) )

            self._allele_calls.flush()

            self._metadata['size'] = 0
            self._metadata['capacity'] = 2000

            self._length = 0
            self._capacity = 2000
            self._last_index = 0


        if self._allele_calls is None:
            raise RuntimeError('Why did this happen?')

    def resize(self):

        Logger.current._log('Resizing array...', depth=3)
        # Let's make sure that we save changes
        # to disk
        self._allele_calls.flush()

        # The size is always 75% of capacity
        # double it to get a capacity of 150% the old capacity
        newsize = self._last_index*2

        # Delete it so we can reload
        del self._allele_calls

        # Get the new path
        array_path = os.path.join(self._dir_path, 'allele_calls_array.memmap' )

        # Create the new shape
        shape = (newsize, self._loci_length)
        
        # Let's remap the new size
        self._allele_calls = np.memmap( 
            array_path, dtype='int32', mode='r+', shape=shape )

        # Store the new capacity:
        self._capacity = newsize

    def save(self):

        attr_mapping = {
            'size': self._length,
            'capacity': self._capacity,
            'runs': self._metadata['runs']+1,
            'last': self._last_index
        }

        meta_data_path = os.path.join(self._dir_path, 'metadata.json')

        with open(meta_data_path, 'w') as f:
            json.dump(attr_mapping, f)

        invalid_path = os.path.join( self._dir_path, 'invalid_indices.json')

        with open(invalid_path, 'w') as f:
            json.dump(list(self._invalid_indices), f)

        keys_mapping_path = os.path.join( self._dir_path, 'keys_mapping.json')

        with open(keys_mapping_path, 'w') as f:
            json.dump(self._keys_to_index, f)

        if self._allele_calls is not None:
            self._allele_calls.flush()

    def validate_calls(self):
        # Here we need to do an integrity check
        # include bionumerics logic

        logger = Logger.current

        expr = Characters.CharSetType( GB_PARAMS['wgmlstname'] )

        charNrs = [ exper.FindChar( char ) for char in \
            exper.ViewGet( GB_PARAMS['scheme'] ).GetCharacters() ]

        mishapkeys = set()

        for entry in bns.Database.Db.Entries:

            if entry.Key in self:

                entryExper = Database.Experiment( entry, 
                    GB_PARAMS['wgmlstname'] )

                charVals = []
                charPresences = []
                entryExper.LoadOrCreate().GetData(charVals, charPresences)
                eCalls = [charVals[i] for i in charNrs]

                if len(self[entry.Key]) != len(eCalls):
                    mishapkeys.add( entry.Key )
                    continue

                for x,y in zip(self[entry.Key], eCalls):
                    if x != y:
                        mishapkeys.add( entry.Key )
                        break

        if len(mishapkeys):
            logger._log('These keys have the wrong allele calls: {}'.format(
                '|'.join(mishapkeys)))

            raise RuntimeError('Stored allele calls are in discordance'
                'with allele calls in the database. Critical error!')

    @staticmethod
    def backup_srcfiles( path ):

        # logging
        logger = Logger.current

        # Log what we are doing
        logger._log('Backing up allele srcfiles', depth=1)

        # Get the files in the current directory
        PATH_DIR = os.path.join( path, 'current' )

        PATH_FILES = [ os.path.join( PATH_DIR, file) for file in \
            os.listdir( PATH_DIR ) ]        

        # Make sure we only have real files
        PATH_FILES = filter( os.path.isfile, PATH_FILES )

        missing_files = []
        # Reference files
        ref_files = ('allele_calls_array.memmap', 'invalid_indices.json',
            'keys_mapping.json', 'metadata.json'
            )

        for x in ref_files:
            if x not in PATH_FILES:
                missing_files.append(x)

        # Log and save the runtimeerror
        if len(missing_files):

            if len(missing_files)==len(ref_files):
                # Nothing to back up so do nothing
                return PATH_DIR

            debugstr = 'Missing srcfiles: {}'.format(
                '|'.join(missing_files))

            logger._log(debugstr)
            raise RuntimeError(debugstr)


        for IN_FILE in PATH_FILES:

            OUT_FILE = os.path.join( path, os.path.basename(IN_FILE) )

            logger._log('Writing to disk', depth=2)

            try:
                with open(IN_FILE, 'rb') as f_in, \
                    open( OUT_FILE, 'wb') as f_out:

                    shutil.copyfileobj(f_in, f_out)

            except:
                logger._log('Backup failed', depth=2)
                raise

            else:
                logger._log('Successfully backed up file: {}'.format(
                    IN_FILE), depth=2)

        return PATH_DIR

    def set(self, key, calls):
        if key not in self._keys_to_index:

            if not isinstance(calls, np.ndarray):
                return False

            # Check to see if we can fit into the array
            if self._last_index+1 <= int(self._capacity*0.75):

                if len(self._invalid_indices):
                    index = self._invalid_indices.pop()
                    self._allele_calls[index] = calls
                    self._keys_to_index[key] = index
                    self._length += 1

                else:
                    self._last_index += 1
                    self._allele_calls[self._last_index] = calls
                    self._keys_to_index[key] = self._last_index
                    self._length += 1

            else:
                # Resize and reload
                self.resize()

                self._last_index += 1
                self._allele_calls[self._last_index] = calls
                self._keys_to_index[key] = self._last_index
                self._length +=1

            return True

        return False

    def __len__(self):
        if self._length == -1:
            return 0

        return self._length

    def __iter__(self):
        for key, value in self._keys_to_index.items():
            yield key, self._allele_calls[value]

    def __getitem__(self, key):
        if key in self._keys_to_index:
            return self._allele_calls[ self._keys_to_index[key] ]

        else:
            raise KeyError('Do not have calls for this key: {}'.format(key))

    def __delitem__(self, key):
        if key in self._keys_to_index:
            self.invalidate(key)
            del self._keys_to_index[key]
            self._length -= 1
        else:
            Logger.current._log('Key: {} does not exist'.format(key))

    def __contains__(self, key):
        return key in self._keys_to_index

def GetDistance(p1, p2):
    common = np.multiply(p1>0, p2>0)
    nCommon = np.sum(common)
    nSame = np.sum(p1[common]==p2[common])
    if nCommon:
        return 100.0 * (float(nCommon) - float(nSame)) / float(nCommon) #nCommon-nSame
    else:
        return 100.

############################## NAMING FUNCTION ###############################
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

    # Do the clustering 
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

################################# CALCULATOR #################################
class Calculator(object):

    def __init__(self, runtimeArgs):

        self._exper = Characters.CharSetType( GB_PARAMS['wgmlstname'] )
        self._wgstfld = Database.Field( GB_PARAMS['wgst'] )
        self._scheme = GB_PARAMS['scheme']
        self._minPres = GB_PARAMS['minpres']
        self._thresholds = list( GB_PARAMS['thresholds'] )
        self._qualityRanges = GB_PARAMS['quality']['Length']
        self._datadir = runtimeArgs['dirpath']
        self._treepath = runtimeArgs['tree']
        self._callspath = runtimeArgs['allele calls']

        self._charNrs = [ self._exper.FindChar( char ) for char in \
            self._exper.ViewGet( self._scheme ).GetCharacters() ]

        self._tree = Tree( len( self._thresholds ) )
        self._alleleCalls = AlleleCalls(self._callspath, len(self._charNrs))

    def DoCalc( self ):

        def CheckCore(entryExper):
            
            charVals = []
            charPresences = []
            entryExper.LoadOrCreate().GetData(charVals, charPresences)

            # How much is present?
            present = sum( charPresences[i] for i in self._charNrs )
            if round(float(present) / float(len(self._charNrs)),2) < self._minPres:
                return None

            else:
                eCalls = np.asarray( [charVals[i] for i in self._charNrs] )
                eCalls.flags.writeable = False
                return eCalls

        def CheckLength(entryExper):

            cs = entryExper.LoadOrCreate()
            seqLength = cs.GetVal( cs.FindName( 'Length' ) )

            if seqLength > max( self._qualityRanges ) or \
                seqLength < min( self._qualityRanges ):
                
                return False

            else:
                return True

        def SaveSettings():
            settings = ET.Element( GB_PARAMS['bnparamkey'] )

            for setting in GB_PARAMS:
                ET.SubElement( settings, setting ).text = str( GB_PARAMS[setting] )

            settings = bns.Database.Db.Info.SaveSetting(
                GB_PARAMS['bnparamkey'], 
                ET.tostring( settings ), 
                True
            )

        # Log function
        logger = Logger.current

        # Tell them what we're doing
        logger._log('Retrieving cached files...', depth=1)
        
        # Load our data into memory
        if not os.path.isdir( self._treepath ):
            logger._log('Loading tree file...', depth=2)
            with open( self._treepath, 'rb' ) as f:
                self._tree.Load( f )

        # Load the allele calls
        logger._log('Loading allele calls', depth=2)
        self._alleleCalls.load()
    
        # Check to make sure the tree and allele calls file match
        logger._log('Checking data integrity...', depth=1)
        if not len( self._alleleCalls ) >= len( self._tree ):
            raise AssertionError('Fatal error, we are missing allele calls'
                ' for keys that have WGS codes in srcfiles' )

        logger._log('Success!', depth=2)

        # Build the tree
        logger._log('Restoring tree state...', depth=1)
        self._tree.BuildTree()
        logger._log('Success!', depth=2)

        # Get all the current names
        namedEntries = list( self._tree.GetNames() )
        belowQC = defaultdict(list)
        
        # Time to assign names:
        logger._log('Beginning calculation @ {}'.format(
            datetime.now()), depth=1)

        for entry in Database.Db.Entries:
            
            # Let's get the experiment
            entryExper = Database.Experiment( entry, 
                GB_PARAMS['wgmlstname'] )

            # Let's get the quality experiment
            qualityExper = Database.Experiment( entry,
                'quality' )
            
            # Let's skip if it doesn't have allele calls
            if not entryExper.IsPresent() or \
                not qualityExper.IsPresent():
                continue

            qcFail = False

            # Check to see if we have a name
            if self._tree.HasName( entry.Key ):
                if not entry.Key in self._alleleCalls:
                    seqQC = CheckLength( qualityExper )

                    if not seqQC:
                        logger._log(
                            'Entry: {} has WGS code, but does not '
                            'have an acceptable sequence length'.format(
                                entry.Key ))

                        logger._log( 'Removing entry: {}'.format(
                            entry.Key ))

                        self._tree.RemoveNamed( entry.Key )

                    else:
                        eCalls = CheckCore( entryExper )
                        if eCalls is not None:
                            retval = self._alleleCalls.set(key, eCalls)
                            if not retval:
                                logger._log('Could not add allele calls for entry:'
                                    ' {}'.format(entry.Key), depth=2)
                        else:
                            logger._log(
                                'Entry: {} has WGS code, but is below'
                                'the {} presence cutoff'.format( 
                                                    entry.Key, self._minPres ))

                            logger._log( 'Removing entry: {}'.format(
                                entry.Key ))

                        self._tree.RemoveNamed( entry.Key )

                continue

            # Check sequence length:
            seqLen = CheckLength( qualityExper )

            # Get the allele calls
            eCalls = CheckCore( entryExper )

            # Continue if we are below QC threshold
            if eCalls is None:
                belowQC[entry.Key].append('CORE')
                qcFail = True

            if not seqLen:
                belowQC[entry.Key].append('LENGTH')
                qcFail = True

            if qcFail:
                continue
            else:
                retval = self._alleleCalls.set(entry.Key, eCalls)

                if not retval:
                    logger._log('Could not add allele calls for'
                        'entry: {}'.format( entry.Key ))

            # Get the distances
            dists = []
            for i in range(len(namedEntries)):
                dist = GetDistance( eCalls, self._alleleCalls[namedEntries[i]] )

                if dist <= self._thresholds[0]:
                    dists.append((dist,i))

            # Calculate tree
            self._tree = CalcName(namedEntries, self._tree, entry.Key,
                                    dists, self._thresholds)

            # Keep track of our new data
            if self._tree.HasName( entry.Key ):
                namedEntries.append( entry.Key )

        # Number of names given
        namesGiven = 0

        # For each key
        for key, name, complete in self._tree.FinalizeCDCNames():

            # Check if the current name is different
            currName = Database.EntryField( key, GB_PARAMS['wgst'] ).Content

            if currName != name:
                attachment = '{0} -> {1}: {2}\n'.format( 
                    currName, name, datetime.now().strftime( "%m-%d-%y" ) )

                found = False
                for flat in Database.Entry(key).GetAttachmentList():
                    if flat['Name'] == 'WGS Code':
                        flat['Attach'].Content += attachment
                        found = True
                        break

                if not found:
                    Database.Entry(key).AddAttachmentFlatText(attachment, saveInRecord=True,
                        name='WGS Code', description='WGS Code Change' )

            
                Database.EntryField( key, GB_PARAMS['wgst'] ).Content = name

            if complete:
                namesGiven += 1

        for key in belowQC:
            criteria = ' , '.join( belowQC[key] )
            Database.EntryField( key, GB_PARAMS['wgst'] ).Content = \
                'FAILED QC: {}'.format( criteria )

        Database.Db.Fields.Save()
        GB_PARAMS['namesgiven'] = namesGiven

        SaveSettings()

        self._alleleCalls.save()

        if not os.path.isdir( self._treepath ):
            DIR_NAME = os.path.dirname( self._treepath )
            FILE_NAME = 'tree_{}.json'.format( 
                datetime.now().strftime("%m-%d-%y@%H-%M-%S"))
            FILE_PATH = os.path.join( DIR_NAME, FILE_NAME )

            with open( FILE_PATH, 'w' ) as f:
                self._tree.Save( f )

            os.remove( self._treepath )

        else:
            FILE_NAME = 'tree_{}.json'.format( 
                datetime.now().strftime("%m-%d-%y@%H-%M-%S"))
            FILE_PATH = os.path.join( self._treepath, FILE_NAME )

            with open( FILE_PATH, 'wb' ) as f:
                self._tree.Save( f )
        
        # Inform the user
        logger._log('Successfully assigned names @ {}'.format(
            datetime.now()))

################################# RUN ########################################

def run( runtime_args ):

    logger = Logger.current

    # If we are here, everything should be good to go
    calc = Calculator( runtime_args )
    logger._log('Successfully created calculator object', depth=1)

    # Run the calculation
    calc.DoCalc()

############################## SETUP #########################################

def setup( **kwargs ):

    RUNTIME_ARGS = {
        'dirpath': None,
        'tree': None,
        'allele calls': None
    }
    
    from ast import literal_eval as leval
    GB_PARAMS_OPS = {
        'bnparamkey': str,
        'organism': str,
        'version': str,
        'scheme': str,
        'thresholds': leval,
        'minpres': float,
        'namesgiven': int,
        'instantiated': str,
        'wgst': str,
        'wgmlstname': str,
        'quality': leval
    }

    logger = Logger.current

    def CreateSettings():
        """
        Create the settings within BN and save it
        """

        settings = ET.Element( GB_PARAMS['bnparamkey'] )

        for setting in GB_PARAMS:
            ET.SubElement( settings, setting ).text = \
                str( GB_PARAMS[setting] )

        settings = bns.Database.Db.Info.SaveSetting(
            GB_PARAMS['bnparamkey'], 
            ET.tostring( settings ), 
            True
        )

    def NewSettings( setting ):

        # Let's make sure at least something is provided:
        if GB_PARAMS.get( setting, None ) is None:
            raise RuntimeError('You cannot provide a new setting'
                ' and say nothing about what it is' )

        # Check if it's an experiment type:
        if setting in Database.Db.ExperimentTypes:
            cst = Characters.CharSetType(setting)

            if isinstance( GB_PARAMS[setting], dict ):

                for char, value in GB_PARAMS[setting].items():
                    if cst.FindChar( char ) < 0:
                        raise RuntimeError('Gave an invalid char: {0} criteria '
                            'for setting: {1}'.format( char, setting ))

            elif isinstance( GB_PARAMS[setting], str ):

                if cst.FindChar( char ) < 0:
                    raise RuntimeError('Gave an invalid char: {0} criteria '
                        'for setting: {1}'.format( char, setting ))
            else:

                raise RuntimeError( 'Do not know what to do with'
                    ' type of setting: {}'.format(type(GB_PARAMS[setting])))

    def SearchSettings():
        """
        Function to make sure we can get every
        setting stored in the database
        """

        # Get the settings from the database
        dbSettings = Database.Db.Info.LoadSetting( GB_PARAMS['bnparamkey'], 
            True )

        if dbSettings:

            logger._log('Found BioNumerics settings', depth=1)

            root = ET.fromstring( dbSettings )

            for setting in GB_PARAMS:
                val = root.find( setting )

                if val is not None:
                    GB_PARAMS[setting] = GB_PARAMS_OPS[setting](val.text)
                else:
                    NewSettings( setting )

            logger._log('Successfully loaded settings', depth=1)

        else:

            logger._log( 'Did not find settings, using default', depth=1)

            # Check to make sure the field exists:
            if len( Fields ) < 1:
                raise RuntimeError('Make sure the database field exists '
                    'to receive the WGS Code')
            
            # Make sure we have the correct field in the database
            elif 'WGST' in Fields:
                GB_PARAMS['wgst'] = 'WGST'

            # You should not be running this script if this error is raised
            else:
                raise RuntimeError('Make sure the database field exists '
                    'to receive the WGS Code')

            # Check to make sure that we have a valid wgMLST char
            # experiment in the database
            currentSchema = Schema.GetCurrent()
            if not currentSchema:
                raise AssertionError('You need to have wgMLST installed on this database'
                        'in order to issue WGS codes')
            else:
                GB_PARAMS['wgmlstname'] = currentSchema.WgMLSTExperTypeName
                cst = Characters.CharSetType( GB_PARAMS['wgmlstname'] )

                if GB_PARAMS['scheme'] not in cst.ViewGetList():
                    raise AssertionError( 'Default scheme not correct OR'
                        ' does not exist in the current database' )

            # Instantiated today!
            GB_PARAMS['instantiated'] = str( date.today() )

            CreateSettings()

            logger._log('Sucessfully created settings', depth=1)

    def ValidateDirectories():

        DIR_PATH = os.path.join( DATA_DIR,
            '{}_nomenclature_srcfiles'.format( GB_PARAMS['organism'] ) )

        if not os.path.isdir( DIR_PATH ):
            logger._log( 'Organism sourcedir does not exist', depth=1 )

            os.mkdir( DIR_PATH )

            for directory in DATA_DIRS:
                path = os.path.join( DIR_PATH,
                        directory )

                os.mkdir( path )

        else:
            for directory in DATA_DIRS:
                path = os.path.join( DIR_PATH, 
                        directory )

                if not os.path.exists( path ):
                    logger._log('Path: {} did not exist'.format(
                        path ))

                    os.mkdir( path )

        RUNTIME_ARGS['dirpath'] = DIR_PATH

        logger._log('Sucessfully validated directories', depth=1)

    def ValidateSrcFiles():

        def Backup( path ):

            # Get the paths
            PATH_DIR = os.path.join( path, 'current')

            PATH_FILES = [ os.path.join( PATH_DIR, file) for file in \
                os.listdir( PATH_DIR ) ]

            # Make sure we only have real files
            PATH_FILES = filter( os.path.isfile, PATH_FILES )

            # This should never happen
            if len( PATH_FILES ) > 1:
                raise AssertionError( 'There are too many objects in path: '
                    '{}'.format( PATH_DIR ) )

            # Make a copy of the current file
            elif len( PATH_FILES ) == 1:
                
                IN_FILE = PATH_FILES[0]
                OUT_FILE = os.path.join( path, os.path.basename( PATH_FILES[0] ))

                # Back, back, back it up!
                logger._log( 'Writing to disk', depth=2 )
                
                # This is cheating, but it works for now

                try:
                    if IN_FILE.endswith( '.gzip' ):

                        with gzip.open( IN_FILE, 'rb' ) as f_in, \
                                gzip.open( OUT_FILE, 'wb') as f_out:

                                shutil.copyfileobj( f_in, f_out )
                    else:

                        with open( IN_FILE, 'rb' ) as f_in, \
                                open( OUT_FILE, 'wb' ) as f_out:
                            
                            shutil.copyfileobj( f_in, f_out )

                except:
                    logger._log('Backup failed', depth=2)
                    raise

                else:
                    logger._log( 'Successfully backed up: {}'.format(
                        IN_FILE, depth=2 ))

                return IN_FILE

            # We might have a new instance
            else:

                OLD_FILES = os.listdir( path )
                OLD_FILES = filter( os.path.isfile, OLD_FILES )

                # Something went really wrong
                if len( OLD_FILES ) > 0:
                    raise AssertionError('We are missing our most current'
                        'file in: {}'.format( path ) )

                else:
                    logger._log( 
                                'Initializing path: {} '
                                'for organism: {}'.format(
                                    path, 
                                    GB_PARAMS['organism']
                                    ),
                                depth=2
                                )

                    return PATH_DIR

        # Add files here in the future, make sure that
        # the files follow the backup scheme in the 
        # above backup function
        files = {
            'tree': os.path.join( RUNTIME_ARGS['dirpath'], 'tree' )
        }

        # For each file, back it up
        for file, path in files.items():
            
            logger._log('Backing up if possible: {}'.format( path ), depth=1 )

            # Set the runtime arguments
            RUNTIME_ARGS[file] = Backup( path )

        RUNTIME_ARGS['allele calls'] = AlleleCalls.backup_srcfiles(
            os.path.join( RUNTIME_ARGS['dirpath'], 'allele_calls' )
            )

    SearchSettings()
    ValidateDirectories()
    ValidateSrcFiles()

    return RUNTIME_ARGS

############################### MAIN #########################################

def main():

    logger = Logger.current

    # Setup the runtime args
    logger._log('Creating habitable environment...')
    runtime_args = setup()

    # Run the calculation
    logger._log('Exceuting naming @ {}'.format(datetime.now()))
    run( runtime_args )

if __name__ == '__main__':

    # Create the logger
    logger = Logger(DATA_DIR)

    # The lockfile
    lock_file = os.path.join(DATA_DIR,'{}_nomenclature_srcfiles'.format( GB_PARAMS['organism'] ),'nomenclature.lock')

    try:

        # Check if we can go forward with the run
        if os.path.exists(lock_file):
            raise RuntimeError( 'Lockfile already exists:'
                'either a process is running, or there was an error' )
        else:
            with open(lock_file, 'w') as f:
                pass

            logger._log('Created lockfile @ {}'.format(datetime.now()))

        # Run the main program
        logger._log('Running main program @ {}'.format(datetime.now()))
        main()

    except:
        
        e = traceback.format_exc()
        logger._log( e )
        traceback.print_exc()

    else:
        # Only if there is no error do we remove the lock file
        # otherwise we will never really know if something went
        # wrong without checking the logs
        os.remove(lock_file)

    finally:
        # Always save the log file
        logger._save()

        # Save the allele calls as well
        if AlleleCalls.current is not None:
            AlleleCalls.current.save()
