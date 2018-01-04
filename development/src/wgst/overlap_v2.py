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
from .database import Database, DatabaseEntry

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

class Database( Database ):

	def AddEntry(self, kwargs):
		entry = DatabaseEntry(kwargs)
		self._entries[ entry._key ] = entry
		return entry

	def QC(self, minPresence):
		poorQuality = 0
		toRemove = []

		for entry in tqdm( self.GetEntries(), desc='QC'):
			nonZero = len( entry._allelecalls[entry._allelecalls > 0] )
			total = len( entry._allelecalls )
			qcVal = True

			if float( nonZero ) / float( total ) < minPresence:
				qcVal = False
				poorQuality += 1

			if not qcVal:
				toRemove.append( entry._key )

		print( 'Poor quality entries: {:d}'.format( poorQuality ) )

		for rm in tqdm( toRemove, 'Removing' ):
			del self._entries[ rm ]

class DatabaseEntry( DatabaseEntry ):
	
	def AddDistance(self, other, dist):
		if not hasattr(self, '_distances'):
			self._distances = {}

		self._distances[ other ] = dist
	
	def GetDistances(self):
		return self._distances

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
		args['organism'],
		args['scheme'],
		args['fields_path'],
		args['outdir'],
		schemespath=args['schemespath'],
		locisize=args['locisize']
		)

	# Load this thing
	dbase.LoadCalls()

	#QC Check
	dbase.QC( args['minPres'] )

	# Calculate the matrix
	for entry, other in tqdm( itertools.combinations( 
		dbase.GetEntries(), 2), desc = 'Creating Matrix' ):
		dist = GetDistance( entry._allelecalls, other._allelecalls )
		entry.AddDistance( other._key, dist )
		other.AddDistance( entry._key, dist )

	# Here we will store the matrix -> { Key: {Other Keys: %Similiarity} }...
	distance_data = {}

	# Build the matrix for shipment
	for entry in tqdm( dbase.GetEntries(), desc='Organizing matrix' ):
		distance_data[ entry._key ] = entry.GetDistances()

	# Save the matrix for future use
	print( 'Dumping matrix...this might take a while' )
	with open( FILE_PATH, 'w' ) as f:

		json.dump(distance_data, f)

	# I'm coming home
	return distance_data