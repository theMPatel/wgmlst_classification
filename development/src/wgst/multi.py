##############################################################
# Nomenclature naming script
# Author: Hannes Pouseele, Milan Patel
# Last Modified: 08/16/2017
# Version: 2.0
# Contact: hannes_pouseele@applied-maths.com, mpatel5@cdc.gov
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

		for thresholds in itertools.product( *toCheck ):

			assert len( thresholds ) == len( thresh[:-1] )

			yield( thresholds )

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

	#======================== NAME EVENT OBJECT =============================#
	class NameEvent(object):
		def __init__(self, oldname, newname, key, info):
			self._oldName = oldname
			self._newName = newname
			self._key = key
			self._info = info

	#======================== NAME HISTORY OBJECT ===========================#
	class NameHistory(object):
		def __init__(self):
			self._events = []
			pass

		def AddEvent(self, oldName, newName, key, info):
			self._events.append(NameEvent( oldName, newName, key, info ))

	#====================== DATABASE OBJECT =================================#
	class Database(object):

		def __init__(self, path): 
			self._entries = {}
			self._dataPath = path
			self.LoadDistanceData()
			self.LoadCalls()

		def AddEntry(self, key ):
			self._entries[ key ] = DbaseEntry( key )
			return self._entries[ key ]

		def GetEntry(self, key):
			return self._entries.get( key, None)

		def GetEntries(self):
			return self._entries

		def LoadCalls(self):

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

			with open(self._dataPath, 'r', encoding='cp1252') as myFile:
				reader = csv.reader( myFile )
				
				# Assumes that you have headers in your file. Comment
				# the next line if that is not the case
				header = tuple( map( str.lower, next(reader) ) )
				refHeader = ( 'key', 'uploaddate' )
				
				for head in refHeader:
					if head not in header:
						raise ValueError('Missing field: {}'.format( head ) )

				for row in reader:

					fldsAndData = { h:f for h,f in zip( header, row ) }

					if len( fldsAndData['key'] ) > 0:
						key = fldsAndData['key']
					else:
						continue

					# Entries should have been loaded by now
					entry = self.GetEntry( fldsAndData['key'] )

					if entry is None: continue

					uploadDate = None
					for matchObj in REUPLOADDATE.keys():
						if matchObj.match( fldsAndData['uploaddate'] ):
							uploadDate = datetime.strptime( 
							fldsAndData['uploaddate'] , REUPLOADDATE[matchObj] ).date()

					if uploadDate is None:
						uploadDate = datetime.strptime( '1/1/2000', "%m/%d/%Y" ).date()

					entry._uploadDate = uploadDate

			print('Pid: {} has loaded qualifiers for {} isolates at {}'.format(
						os.getpid(), len( self._entries ), datetime.now() ))

		def CreateSubset(self):

			allSet = list( self.GetEntries().values() )
			allSet.sort( key=lambda x: x.UploadDate() )
			startingSet = allSet[ :(len(allSet)//2) ]
			addingSet = allSet[ (len(allSet)//2): ]
			addingSet = [ addingSet[i:i+30] for i in range(0, len(addingSet), 30) ]

			return (( startingSet, addingSet ))

		def LoadDistanceData(self):

			lenDatabase = len( Main.distance_matrix )

			for k, v in Main.distance_matrix.items():
				
				newEntry = self.AddEntry( k )
				
				assert len( v ) == lenDatabase - 1

				newEntry.AddDistanceData( v )

			print('Pid: {} has loaded distance matrix at {}'.format( os.getpid(), datetime.now() ))

		def Threshold(self, value = None):
			if not value:
				return self._thresholds

			else:
				assert isinstance( value, list )
				self._thresholds = value

		def TotalIsolates(self):
			return len( self._entries )

	#====================== DATABASE ENTRY OBJECT ===========================#
	class DbaseEntry(object):

		def __init__(self, key ):

			self._key = key
			self._uploadDate = None
			self._wgst = None
			self._qcStatus = True
			self._wgstHistory = []

		def Key(self):
			return self._key

		def AddDistanceData(self, value):
			assert isinstance( value, dict )
			self.distances = value

		def Distance(self, other):
			if isinstance(other, str):
				dist = self.distances.get(other, None)
				
				if dist == None:
					raise ValueError( 'Distance Accession Error: {}-{}'.format(self._key, other ) )
				else:
					return dist

			elif isinstance(other, type(self) ):
				dist = self.distances.get(other.Key(), None)
				
				if dist == None:
					raise ValueError( 'Distance Accession Error: {}-{}'.format(self._key, other.Key() ) )
				else:
					return dist

			else:
				raise TypeError( 'Invalid object for dist accession: {}'.format( repr( other ) ) )

		def GetDistances(self):
			return self.distances

		def QC(self):
			return self._qcStatus

		def UploadDate(self):
			return self._uploadDate

		def Wgst(self, value = None, reset = False):

			if reset:
				self._wgst = None
				return

			if not value:
				return self._wgst

			else:
				assert isinstance( value, str )
				self._wgst = value
				self._wgstHistory.append( value )

		def WgstHistory(self):
			return self._wgstHistory
		
		def __eq__(self, other):

			if isinstance( other, type(self) ):
				return self._key == other.Key()

			elif isinstance( other, str ):
				return self._key == other

			else:
				raise TypeError('Invalid comparison: {}'.format(other) )

		def __repr__(self):
			return '{0} | {1}'.format( self._key, self._uploadDate )
		
		def __str__(self):
			return self._key

	#======================= NAMES OBJECT ===================================#
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
			idx = [p[len(parts)] for k, p in self._allnames.items() if p[:len(parts)]==parts and len(p)>len(parts)]
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

	#=========================== NAMING FUNCTION ============================#
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
			if len(closestClusters) == 1:
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

			self.thresholds = thresholds
			self._entryBase = database

		def DoValidation(self):

			existingNames = Names()
			

			# Create starting and adding sets return the updated Names object to track history


			self._startingSet, addingSet = self._entryBase.CreateSubset()
			
			print('Pid: {} is initializing clusters at: {}'.format( os.getpid(), datetime.now() ) )
			existingNames = self.DoCalc( self._startingSet, existingNames, 'Initializing' )

			self._wgstHistory = WgstHistory( self.thresholds )

			print('Pid: {} is incrementing at: {}'.format( os.getpid(), datetime.now() ))
			for i in range( len( addingSet) ):
				existingNames = self.DoCalc( addingSet[i], existingNames, 'Incrementing', i+1 )
			
		def DoCalc(self, selection, existingNames, status, iteration = 0):

			namedEntries = []
			nameFrequencies = {}

			for entry in self._entryBase.GetEntries().values():
				if existingNames.HasResolvedName( entry.Key() ):
					namedEntries.append( entry.Key() )
					nameFrequencies[ entry.Wgst() ] = nameFrequencies.get( entry.Wgst(), 0) + 1

			#for each entry, calculate the name 
			for entry in selection:
				
				qcStatus = entry.QC()
					
				if qcStatus and entry.Key() in namedEntries:
					continue
					
				if existingNames.HasName( entry.Key() ):
					existingNames.DropName( entry.Key() )

				dists = [ entry.Distance( e ) for e in namedEntries ]
				
				#calculate the name of the entry
				existingNames = CalcName(namedEntries, existingNames, entry.Key(), dists, self.thresholds, qcStatus)
				
				#keep track of the data
				if existingNames.HasResolvedName( entry.Key() ):
					namedEntries.append( entry.Key() )

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
					self._wgstHistory.TrackChanges( oldName, newName, iteration )

				dbEntry.Wgst( value = newName )
			
			return( existingNames )

	db = Database( Main.args['fields_path'] )
	
	print( 'Pid: {} has loaded database at: {}'.format( os.getpid(), 
		datetime.now() ))
	
	print( 'Pid: {} is requesting NEXT at: {}'.format( os.getpid(), 
		datetime.now() ))
	
	pipe.send(( REQUEST_GET_NEXT, ))

	while True:

		msg = pipe.recv()

		if msg[0] == RESPONSE_DATA:

			thresholds = msg[1]
			thresholds = list( thresholds )
			thresholds.sort( key = lambda x: -x )

			db.Threshold( value = thresholds )

			for entry in db.GetEntries().values():
				entry.Wgst( reset = True )

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