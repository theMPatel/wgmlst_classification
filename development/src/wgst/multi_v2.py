##############################################################
# Nomenclature naming script
# Author: Milan Patel
# Last Modified: 08/16/2017
# Version: 2.1
# Contact: mpatel5@cdc.gov
##############################################################
from multiprocessing import Process, Pipe

import os
import sys
import csv
from collections import defaultdict

from datetime import datetime, timedelta
import json
import random
import re
import numpy as np

REQUEST_GET_NEXT = '__request_next__'
RESPONSE_DATA = '__response_data__'
RESPONSE_DONE = '__response_done__'
RETURN_DATA = '__return_data__'

def DataServer( pps, ph ):

	import csv
	class OptimizeThresholds(object):

		current = None

		def __init__(self):
			self._history = []
			self._scheme = Main.args['scheme']
			self._thresholds = Main.args['thresholds']
			self._outDir = Main.args['outdir']
			OptimizeThresholds.current = self

		def AddThreshold(self, results):
			self._history.append( results )

			print('Pid: {} is sending home at {}'.format( os.getpid(), datetime.now() ))
			ph.send(( RETURN_DATA, self._history[-1] ))

		def GetThresholds(self):
			return self._thresholds

		def Save(self):
			FILE_PATH = os.path.join( self._outDir, 
				'{}_backup.csv'.format( self._scheme ) )

			with open( FILE_PATH, 'w' ) as f:
				writer = csv.writer(f)
				writer.writerows( self._history )

	def GenNextThresholds():
		import itertools

		toCheck = []
		thresh = OptimizeThresholds.current.GetThresholds()

		for threshSpace in thresh[:-1]:
			if isinstance( threshSpace, tuple ):
				toCheck.append( list( np.arange( threshSpace[0], \
					threshSpace[1], thresh[-1] ) ) )
			elif isinstance( threshSpace, list ):
				toCheck.append( threshSpace )
			else:
				raise RuntimeError('Invalid thresholds, this should have been caught earlier')
		
		wallclock = datetime.now()
		for thresholds in itertools.product( *toCheck ):

			assert len( thresholds ) == len( thresh[:-1] )

			if datetime.now() < wallclock + timedelta(days=2):
				yield( thresholds )
			else:
				raise StopIteration

	def Server( generator ):

		optimumThresh = OptimizeThresholds.current
		
		print('Beginning transmission at: {}'.format( datetime.now() ))

		while len( pps ) > 0:
			pipes_to_remove = []

			for pipe in pps:

				if pipe.poll():
					msg = pipe.recv()

					if msg[0] == REQUEST_GET_NEXT:

						try:
							pipe.send(( RESPONSE_DATA, next(generator) ))
						except StopIteration:
							pipe.send(( RESPONSE_DONE, ))

					elif msg[0] == RETURN_DATA:
						results = msg[1]
						optimumThresh.AddThreshold( results )

					elif msg[0] == RESPONSE_DONE:
						pipes_to_remove.append(pipe)

			for p in pipes_to_remove:
				pps.remove(p)

		optimumThresh.Save()
		ph.send(( RESPONSE_DONE, ))
		return

	OptimizeThresholds()
	Server( GenNextThresholds() )

	return

def DataConsumer( pipe ):
	from .tree import Tree, Node, NamedNode
	from .database import Database, DatabaseEntry

	#====================== DATABASE OBJECT =================================#
	class Database(Database):

		def AddEntry(self, **kwargs):
			entry = DatabaseEntry(kwargs)
			self._entries[ entry._key ] = entry
			return entry

		def LoadDistanceData(self):

			self._len = len( Main.distance_matrix )

			for k, v in Main.distance_matrix.items():
				
				newEntry = self.AddEntry( key=k )
				
				assert len( v ) == self._len - 1

				newEntry.AddDistanceData( v )

			print('Pid: {} has loaded distance matrix at {}'.format( os.getpid(), datetime.now() ))
		
		def CreateSubset(self):

			allSet = list( self.GetEntries() )
			allSet.sort( key=lambda x: x._uploaddate )

			startingSet = allSet[ :(len(allSet)//2) ]
			addingSet = allSet[ (len(allSet)//2): ]
			addingSet = [ addingSet[i:i+30] for i in range(0, len(addingSet), 30) ]
			return (( startingSet, addingSet ))

		def ExtractData(self, filetype, header, reader ):
			refHeader = ( 'key', 'uploaddate' )

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

			count = 0
			for row in reader:

				if row[0] == 'key':
					continue

				entry = self.GetEntry( row[0] )

				if entry is None:
					continue

				for h, fld in zip(header, row):

					if h not in refHeader:
						continue

					elif h == 'uploaddate':

						h = '_' + h

						for matchObj in REUPLOADDATE.keys():
							matched = matchObj.match( fld )

							if matched:
								upload = datetime.strptime(
									fld,
									REUPLOADDATE[matchObj] ).date()

								break

						if matched is None:

							if REUPLOADDATE_5.match( fld ):
								upload = datetime.strptime(
									fld,
									"%Y" ).date()
							else:
								upload = datetime.strptime(
									'1/1/2000',
									'%m/%d/%Y' ).date()

						setattr(entry, h, upload)

						break

				count +=1

			assert count == self._len

		def Reset(self):
			for entry in self._entries.values():
				del entry.Wgst
	
		@property
		def Threshold(self):
			return self._thresholds

		@Threshold.setter
		def Threshold( self, value ):
			self._thresholds = value

		def TotalIsolates(self):
			return len( self._entries )

		def __len__(self):
			return len(self._entries)

	#====================== DATABASE ENTRY OBJECT ===========================#
	class DatabaseEntry(DatabaseEntry):

		def AddDistanceData(self, data):
			assert isinstance(data, dict)
			self._distances = data
		
		
		def Distance(self, other):
			if isinstance(other, str):
				dist = self._distances.get(other, None)
				
				if dist == None:
					raise ValueError( 'Distance Accession Error: {}-{}'.format(self._key, other ) )
				else:
					return dist

			elif isinstance(other, type(self) ):
				dist = self._distances.get(other._key, None)
				
				if dist == None:
					raise ValueError( 'Distance Accession Error: {}-{}'.format(self._key, other.Key() ) )
				else:
					return dist

			else:
				raise TypeError( 'Invalid object for dist accession: {}'.format( repr( other ) ) )

		def GetDistances(self):
			return self._distances

		@property
		def Wgst(self):
			return self._wgst

		@Wgst.setter
		def Wgst(self, value):
			assert isinstance( value, str )
			self._wgst = value
			self._wgstHistory.append( value )

		@Wgst.deleter
		def Wgst(self):
			self._wgst = None
			self._wgstHistory = []

		@property
		def WgstHistory(self):
			return self._wgstHistory
		
		def __repr__(self):
			if not hasattr(self, '_key') or \
				not hasattr(self, '_uploaddate'):

				return 'Failed repr for key: {}'.format( self )
			else:
				return '{0} | {1}'.format( self._key, self._uploaddate )
		
		def __str__(self):
			return self._key

	#=========================== NAMING FUNCTION ============================#
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

	#========================== WGST HISTORY CLASS ==========================#
	class WgstHistory(object):

		def __init__(self, thresholds):
			
			self._thresholds = thresholds

			self._totalChanges = [0]*len( self._thresholds )
			self._totalChangesWithTime = defaultdict(list)

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

		def TrackChanges( self, oldVal, newVal, iteration ):

			theSame = self.CalculateDifference( oldVal, newVal )
			
			self._totalChanges[ theSame ] += 1
			self._totalChangesWithTime[ iteration ].append( theSame )

	#========================= CALCULATOR OBJECT ============================#
	class Calculator(object):
		
		def __init__(self, thresholds, database ):

			self._thresholds = thresholds
			self._entryBase = database

		def DoValidation(self):

			# Initialize our constants
			self._tree = Tree( len( self._thresholds ) )
			self._wgstHistory = WgstHistory( self._thresholds )
			self._startingSet, addingSet = self._entryBase.CreateSubset()
			self._namedEntries = []

			# Create starting and adding sets return the updated Names object to track history
			print('Pid: {} is initializing clusters at: {}'.format( 
					os.getpid(), datetime.now() ) )

			# Do the calculation
			self.DoCalc( 
						self._startingSet, 
						'Initializing' )

			# Notify we are starting
			print('Pid: {} is incrementing at: {}'.format( 
				os.getpid(), datetime.now() ))
			
			# Increment
			for i in range( len( addingSet) ):
				self.DoCalc( addingSet[i], 'Incrementing', iteration=i+1 )
			
		def DoCalc(self, selection, status, iteration = 0):

			nameFrequencies = {}

			for entry in self._namedEntries:
				nameFrequencies[ self._tree.GetStrName(entry) ] = \
					nameFrequencies.get( self._tree.GetStrName(entry), 0) + 1

			#for each entry, calculate the name 
			for entry in selection:
					
				if entry._key in self._namedEntries:
					continue
				
				dists = []
				for i in range(len(self._namedEntries)):
					dist = entry.Distance(self._namedEntries[i])

					if dist <= self._thresholds[0]:
						dists.append((dist,i))
				# dists = [ entry.Distance( e ) for e in self._namedEntries ]
				
				#calculate the name of the entry
				CalcName( self._namedEntries, self._tree, entry._key, dists, self._thresholds )
				
				#keep track of the data
				if self._tree.HasName( entry._key ):
					self._namedEntries.append( entry._key )
					nameFrequencies[ self._tree.GetStrName( entry._key ) ] = nameFrequencies.get( self._tree.GetStrName( entry._key ), 0) + 1

			if iteration == 0:
				for entry in self._namedEntries:
					self._entryBase.GetEntry( entry ).Wgst = \
						self._tree.GetStrName( entry )

			for entry in self._namedEntries:
				dbEntry = self._entryBase.GetEntry( entry )

				oldName = dbEntry.Wgst
				newName = self._tree.GetStrName( dbEntry._key )

				if oldName == newName:
					continue
				
				if iteration > 0 and oldName and nameFrequencies[oldName] > 1:
					self._wgstHistory.TrackChanges( oldName, newName, iteration )

				dbEntry.Wgst = newName

	db = Database( 
					Main.args['organism'],
					Main.args['scheme'],
					Main.args['fields_path'],
					Main.args['outdir'] )

	db.LoadDistanceData()
	db.LoadCalls()
	
	print( 'Pid: {} has loaded database of size {} at: {}'.format( os.getpid(), 
		len(db), datetime.now() ))
	
	print( 'Pid: {} is requesting NEXT at: {}'.format( os.getpid(), 
		datetime.now() ))
	
	pipe.send(( REQUEST_GET_NEXT, ))

	while True:

		msg = pipe.recv()

		if msg[0] == RESPONSE_DATA:

			thresholds = msg[1]
			thresholds = list( thresholds )
			thresholds.sort( key = lambda x: -x )

			db.Threshold = thresholds

			db.Reset()

			calc = Calculator( thresholds, db )
			calc.DoValidation()

			total = calc._wgstHistory.GetTotal()

			returnData = [ x for x in thresholds ]
			returnData.append( total )
			returnData = list( map( float, returnData ) )

			pipe.send(( RETURN_DATA, returnData ) )
			print('Pid: {} is requesting NEXT at {}'.format( os.getpid(), datetime.now() ))
			pipe.send(( REQUEST_GET_NEXT, ))

		elif msg[0] == RESPONSE_DONE:
			pipe.send(( RESPONSE_DONE, ))
			return

		else:
			pipe.send(( REQUEST_GET_NEXT, ))
			print( 'Pid: {} is requesting next at {}'.format( os.getpid(), datetime.now() ))

def WriteServer(pipe):
	
	scheme = Main.args['scheme']
	outPath = Main.args['outdir']
	resultsPath = os.path.join( outPath, '{}_optimized.txt'.format( scheme ) )

	print('Pid: {} is ready for writing at {}'.format( os.getpid(), datetime.now() ))
	while True:

		if pipe.poll():

			msg = pipe.recv()

			if msg[0] == RETURN_DATA:
				
				print('Pid: {} is rcving at at {}'.format( os.getpid(), datetime.now() ))
				if not os.path.isfile( resultsPath ):
					with open( resultsPath, 'w' ) as f:
						writer = csv.writer(f)
						writer.writerow( msg[1] )
				else:
					with open( resultsPath, 'a') as f:
						writer = csv.writer(f)
						writer.writerow( msg[1] )

				print('Pid: {} has written to file at {}'.format( os.getpid(), datetime.now() ))

			elif msg[0] == RESPONSE_DONE:
				return

			else:
				raise ValueError('Invalid message for this process')

def Main(args):

	# Output dir
	OUT_DIR = os.path.join( args['outdir'], 'optimized' )
	if not os.path.isdir( OUT_DIR ):
		os.mkdir( OUT_DIR )

	args['outdir'] = OUT_DIR

	# Save the args for sharing
	Main.args = args
	
	# Number of cores
	numConsumers = args['cores'] - 3

	# Assign the distance_matrix to the current namespace, children
	# will inherit it
	Main.distance_matrix = args['matrix']

	# Communication for the children
	pipes = [ Pipe( duplex=True ) for x in range( numConsumers ) ]

	# Communication between the threshold server and the write server
	pipeHome = Pipe( duplex=True )

	# Initialize the threshold server
	GenServer = Process( target = DataServer, args=(( [ p[1] for p in pipes ], 
		pipeHome[1] )))

	# Initialize the write server
	WriteOut = Process( target = WriteServer, args=(( pipeHome[0] ), ))

	# Make the children
	Children = [ Process( target = DataConsumer, args=(p[0],)) for p in pipes ]

	# Let's start our servers to they are ready to receive requests
	GenServer.start()
	WriteOut.start()

	# Start the children
	for proc in Children: proc.start()

	# Wait for everyone to finish
	for proc in Children: proc.join()
	GenServer.join()
	WriteOut.join()

if __name__ == "__main__":

	Main({
				'cores':        10,
				'matrix':       {},
				'scheme':       'whole',
				'thresholds':   [( 1.0312822275696116,  1.2031625988312136),( 0.77346167067720872, 0.94534204193881066),( 0.51564111378480582, 0.68752148504640775),( 0.25782055689240291, 0.42970092815400485 ),(0.0, 0.17188037126160194),0.00859401856308009625300790649708],
				'fields_path':  {'train':'/scicomp/home/nft9/projects/wgst/stability/dev/allele_calls/salm_allele_calls_all.csv'},
				'outdir':       '/scicomp/home/nft9/projects/wgst/stability/dev/salm_results/',
				'organism': 'saen'
				})