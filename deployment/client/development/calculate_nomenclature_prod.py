#//$autorun;event=CreateMainWin

##############################################################
# Nomenclature Pattern Naming Script
# Author: Milan Patel
# Last Modified: 08/31/2017
# Version: 1.1
# Contact: mpatel5@cdc.gov
# Deployed: wgListeria (production)
##############################################################

# Stdlib imports
import os
import sys
import base64
import gzip
import json
import shutil
import traceback

from datetime import datetime, date
from collections import defaultdict

# External imports
import numpy as np
import bns
from wgMLST_Client.wgMLSTSchema import Schema
import xml.etree.cElementTree as ET

# BioNumerics constants
import bns
Dlg = bns.Windows.XmlDlgBuilder	
Database = bns.Database
Fields = Database.Db.Fields
MessageBox  = bns.Util.Program.MessageBox
ConfirmBox  = bns.Util.Program.ConfirmBox
Characters = bns.Characters

# Parameters that are required for this organism
GB_PARAMS = {
    'bnparamkey': 'hierarchical_strain_nomenclature',
    'organism': 'LMO',
    'version': '1.0',
    'scheme': '__Core__',
    'thresholds': (4.050815987,2.914238135,2.046511241,
        1.068159867,0.041631973),
    'minpres': 0.95,
    'namesgiven': 0,
    'instantiated': None,
    'wgst': None,
    'wgmlstname': None,
    'quality': {'Length':[2.8e6, 3.2e6]},
    'nosave': False,
}

# Data directories
DATA_DIR = Database.Db.Info.SourceFilesDir
DATA_DIRS = [ 'tree', 'tree/current', 'allele_calls', 'allele_calls/current', 'logs' ]

#========================= LOGGER CLASS =====================================#
class Logger(object):
    """
    Logger class, will be used to log things during the naming process
    """

    CURRENT = None

    def __init__(self, logDir):
        path = [
            '{}_nomenclature_srcfiles'.format( 
            GB_PARAMS['organism'] ),
            'logs',
            'wgst_log_{}.txt'.format( 
                datetime.now().strftime("%m-%d-%y@%H-%M-%S"))
             ]
        self._fileName = os.path.join( logDir, *path )
        self._logStr = ''
        Logger.CURRENT = self

    def _log(self, outStr):
        if not isinstance(outStr, str):
            
            try:
                self._logStr += str(outStr) + '\n'

            except:
                e = traceback.format_exc()
                self._logStr += e + '\n'

        else:
            self._logStr += outStr + '\n'

    def _save(self):

        with open(self._fileName, 'w') as f:
            f.write( self._logStr)

class WgstDlg( Dlg.Dialogs ):
    
    def __init__( self ):
        Dlg.Dialogs.__init__( self, 'WGSTStrainNames' )
        
        self._organism = Dlg.StaticText( str(GB_PARAMS['organism'] ) )
        self._scheme = Dlg.StaticText( str( GB_PARAMS['scheme'] ) )
        self._namesGiven = Dlg.StaticText( str(GB_PARAMS['namesgiven'] ) )
        self._dateInstantiated = Dlg.StaticText( str( GB_PARAMS['instantiated'] ) )
        self._currVersion = Dlg.StaticText( str( GB_PARAMS['version'] ) )
        self._minPres = Dlg.Input('minpres', tp='float', default=95.0, remembersettings=False)
        self._noSaveCheck = Dlg.Check('nosave', '', 
                remembersettings=False, default = 0, OnChange = self.EnableOptions )
        
        grid = [
                        ['Organism: ', self._organism ],
                        ['Scheme: ', self._scheme ],
                        ['Current Version: ', self._currVersion ],
                        #['Date initiated: ', self._dateInstantiated ],
                        #['WGS Codes issued: ', self._namesGiven ],
                        ['Name unsatisfactory: ', self._noSaveCheck ],
                        ['Core minimum presence: ', Dlg.Cell([[self._minPres, '%']])]
                    ]
        
        simple = Dlg.SimpleDialog( grid, onStart=self.OnStart, onOk=self.OnOk )
        self.AddSimpleDialog('Issue WGS Codes', simple )

    def EnableOptions(self, args):
        self._minPres.Enabled = not self._minPres.Enabled

    def OnStart( self, args ):
        self._minPres.Enabled = False
    
    def OnOk( self, args ):
        if self._minPres.Enabled:
            GB_PARAMS[ 'minpres' ] = (float(self._minPres.GetValue())/100.0)
            GB_PARAMS[ 'nosave' ] = True

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
        return key in self._names

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

class AlleleCalls( object ):

    def __init__(self):
        self._alleleCalls = {}

    def Load( self, flobj ):
        json_bytes = flobj.read()
        json_str = json_bytes.decode( 'utf-8' )
        self._alleleCalls = json.loads( json_str )
        for key in self._alleleCalls:
            self._alleleCalls[key] = np.asarray(self._alleleCalls[key], dtype=int)

    def Save( self, flobj ):
        for key in self._alleleCalls:
            self._alleleCalls[key] = list(self._alleleCalls[key])
        json_str = json.dumps(self._alleleCalls )
        json_bytes = json_str.encode( 'utf-8' )
        flobj.write(json_bytes)

    def Add( self, key, calls ):
        self._alleleCalls[ key ] = calls

    def GetCalls( self, key ):
        return self._alleleCalls.get( key, None )

    def HasKey( self, key ):
        return key in self._alleleCalls

    def __len__(self):
        return len( self._alleleCalls )

def GetDistance(p1, p2):
    common = np.multiply(p1>0, p2>0)
    nCommon = np.sum(common)
    nSame = np.sum(p1[common]==p2[common])
    if nCommon:
        return 100.0 * (float(nCommon) - float(nSame)) / float(nCommon) #nCommon-nSame
    else:
        return 100.

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
    namedEntries = list( named )
    
    for level, threshold in enumerate( thresholds, 1 ):

        # Create the clusters at this level:
        closestClusters = set()
        toRemove = []
        for i, dist in enumerate( distances ):
            
            if dist <= threshold:
                keyToGet = namedEntries[i]
                partial = tuple( tree.GetPart( keyToGet, level ) )
                closestClusters.add( partial )

            else:
                toRemove.append( i )

        # If something isn't similar at this level
        # There's no way it will be similar at the lower
        # levels, let's remove those
        for index in reversed(toRemove):
            del namedEntries[ index ]
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
        self._alleleCalls = AlleleCalls()

    def DoRefresh( self, args ):
        pass

    def DoCalc( self, args ):

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

            if seqLength >= max( self._qualityRanges ) or \
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

        # Will be used for the progress bar
        comm = args.get('communication', bns.Windows.CalcCommunication() )

        # Log function
        logger = Logger.CURRENT

        # Tell them what we're doing
        comm.SetMessage("Retrieving view: {}".format( self._scheme ) )
        
        # Load our data into memory
        if not os.path.isdir( self._treepath ):
            comm.SetMessage( "Loading current tree" )
            with open( self._treepath, 'rb' ) as f:
                self._tree.Load( f )

        if not os.path.isdir( self._callspath ):
            comm.SetMessage( "Loading historical allele calls")
            with gzip.open( self._callspath, 'rb' ) as f:
                self._alleleCalls.Load( f )

        # Check to make sure the tree and allele calls file match
        comm.SetMessage( "Checking data integrity")
        if not len( self._alleleCalls ) >= len( self._tree ):
            raise AssertionError('Fatal error, we are missing allele calls'
                'for keys that have WGS codes in srcfiles' )

        # Build the tree
        comm.SetMessage( "Building tree" )
        self._tree.BuildTree()

        # Get all the current names
        namedEntries = list( self._tree.GetNames() )
        belowQC = defaultdict(list)
        selection = set()
        
        # Time to assign names:
        comm.SetMessage( 'Calculating and assigning WGS codes' )

        nSelected = len( Database.Db.Selection )
        for i, entry in enumerate( Database.Db.Selection ):
            
            comm.SetProgress( i, nSelected )

            # Let's get the experiment
            entryExper = Database.Experiment( entry, 
                GB_PARAMS['wgmlstname'] )

            # Let's get the quality experiment
            qualityExper = Database.Experiment( entry,
                'quality' )

            # add to selection tracker
            selection.add( entry.Key )
            
            # Let's skip if it doesn't have allele calls
            if not entryExper.IsPresent() or \
                not qualityExper.IsPresent():
                continue

            qcFail = False

            # Check to see if we have a name
            if self._tree.HasName( entry.Key ):
                if not self._alleleCalls.HasKey( entry.Key ):
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
                            self._alleleCalls.Add( key, eCalls )
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
                self._alleleCalls.Add( entry.Key , eCalls )

            # Get the distances
            dists = [GetDistance( eCalls, self._alleleCalls.GetCalls(named) ) \
                        for named in namedEntries ]

            # Calculate tree
            self._tree = CalcName(namedEntries, self._tree, entry.Key,
                                    dists, self._thresholds)

            # Keep track of our new data
            if self._tree.HasName( entry.Key ):
                namedEntries.append( entry.Key )

        if GB_PARAMS['nosave']:
            
            comm.SetMessage( 'Saving codes to your desktop' )
            import csv
            
            PATH = os.path.normpath( os.path.expanduser( '~/Desktop' ) )
            FILE_NAME = 'wgs_code_calculation-{}.csv'.format( 
                datetime.now().strftime("%m-%d-%y@%H-%M-%S") )
            FILE_PATH = os.path.join( PATH, FILE_NAME )
            
            with open( FILE_PATH, 'wb' ) as f:
                writer = csv.writer( f )
                writer.writerow(['Key', 'WGST', 'UploadDate'])
                
                for key, value, complete in self._tree.FinalizeCDCNames():
                    if key in selection:
                        writer.writerow( [key, value, Database.EntryField(key, 'uploaddate').Content] )

                for key in belowQC:
                    if key in selection:
                        criteria = ' , '.join( belowQC[key] )
                        writer.writerow( [key, 'FAILED QC: {}'.format( criteria ) ] )

        else:
            comm.SetMessage( 'Saving codes to the database' )
            
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

            if not os.path.isdir( self._callspath ):
                with gzip.open( self._callspath, 'wb') as f:
                    self._alleleCalls.Save( f )
            else:
                FILE_PATH = os.path.join( self._callspath, 'calls.gzip')
                with gzip.open( FILE_PATH, 'wb' ) as f:
                    self._alleleCalls.Save( f )

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
        MessageBox('', 'Successfully assigned WGS Codes', '')

        # Open the file if there is a no save:
        if GB_PARAMS['nosave']:
            os.startfile( FILE_PATH )

def Setup( **kwargs ):

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
        'quality': leval,
        'nosave': leval
    }

    logger = Logger.CURRENT

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

            logger._log('Found BioNumerics settings')

            root = ET.fromstring( dbSettings )

            for setting in GB_PARAMS:
                val = root.find( setting )

                if val is not None:
                    GB_PARAMS[setting] = GB_PARAMS_OPS[setting](val.text)
                else:
                    NewSettings( setting )

            logger._log('Successfully loaded settings')

        else:

            logger._log( 'Did not find settings, using default')

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

            logger._log('Sucessfully created settings')

    def ValidateDirectories():

        DIR_PATH = os.path.join( DATA_DIR,
            '{}_nomenclature_srcfiles'.format( GB_PARAMS['organism'] ) )

        if not os.path.isdir( DIR_PATH ):
            logger._log( 'Organism sourcedir does not exist' )

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

        logger._log('Sucessfully validated directories')

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
                logger._log( 'Writing to disk' )
                
                # This is cheating, but it works for now
                if IN_FILE.endswith( '.gzip' ):

                    with gzip.open( IN_FILE, 'rb' ) as f_in, \
                            gzip.open( OUT_FILE, 'wb') as f_out:

                            shutil.copyfileobj( f_in, f_out )
                else:

                    with open( IN_FILE, 'rb' ) as f_in, \
                            open( OUT_FILE, 'wb' ) as f_out:
                        
                        shutil.copyfileobj( f_in, f_out )

                logger._log( 'Successfully backed up: {}'.format(
                    IN_FILE ))

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
                                    )
                                )

                    return PATH_DIR

        # Add files here in the future, make sure that
        # the files follow the backup scheme in the 
        # above backup function
        files = {
            'tree': os.path.join( RUNTIME_ARGS['dirpath'], 'tree' ),
            'allele calls': os.path.join( 
                RUNTIME_ARGS['dirpath'], 'allele_calls' )
        }

        # For each file, back it up
        for file, path in files.items():
            
            logger._log('Backing up if possible: {}'.format( path ) )

            # Set the runtime arguments
            RUNTIME_ARGS[file] = Backup( path )

    SearchSettings()
    ValidateDirectories()
    ValidateSrcFiles()

    return RUNTIME_ARGS

def Run( runtimeArgs ):

    # The main window
    winID = 1

    # If we are here, we are ready to issue WGS codes
    calc = Calculator( runtimeArgs )

    # Run the calculation
    bns.Windows.BnsWindow(winID).StartAsyncCalculation(
            calc.DoCalc, calc.DoRefresh, async=False )

def Main( args ):

    # The logger!!
    logger = Logger(DATA_DIR)

    try:

        # Let's set ourselves up for success!
        runtimeArgs = Setup()

        # Let's make the dialog box:
        DialogBox = WgstDlg()
        if not DialogBox.Show():
            bns.Stop()

        # Start the calculation
        Run( runtimeArgs )

    except:
        e = traceback.format_exc()
        logger._log( e )
        MessageBox( 'Error', e, '' )
    
    finally:
        logger._save()

if __name__ == '__main__':

    try:
        bns.Windows.BnsWindow.AddCustomSection(
                        WindowClass='main',
                        sectionID='WgMlst',
                        sectionName="wgMLST",
                        parentSectionID='WgsTools'
        )

        bns.Windows.BnsWindow.AddCustomSection(
                        WindowClass='main',
                        sectionID='WgMlstStrainNaming',
                        sectionName='wgMLST strain naming',
                        parentSectionID='WgMlst'
        )

        bns.Windows.BnsWindow.AddCustomCommand(
                        WindowClass='main',
                        sectionID='WgMlstStrainNaming',
                        commandID='assign_wgst',
                        caption='Assign WGS Codes',
                        description='Assign WGS Codes',
                        execute=Main
        )


        # bns.Windows.BnsWindow.AddCustomSection(
        #                 WindowClass='main',
        #                 sectionID='WGS codes',
        #                 sectionName='WGS codes',
        #                 parentSectionID=''
        # )

        # bns.Windows.BnsWindow.AddCustomCommand(
        #                 WindowClass='main',
        #                 sectionID='WGS codes',
        #                 commandID='assign_wgst',
        #                 caption='Assign WGST',
        #                 description='Assign WGST',
        #                 execute=Main
        # )
    except Exception:
        pass
