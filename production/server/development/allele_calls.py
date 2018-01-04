import traceback
import csv
import os
import sys
import numpy as np
import json
import traceback
from tqdm import *

CORE_LOCI = 1748
_ORGANISM = 'lmo'

class Logger(object):

	current = None

	def __init__(self):
		Logger.current = self

	def _log(self, string):
		print(string)

class AlleleCalls(object):

	def __init__(self, dir_path):

		self._allele_calls = None
		self._dir_path = dir_path
		self._length = -1
		self._capacity = -1
		self._last_index = -1
		self._keys_to_index = {}
		self._invalid_indices = set()

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


			shape = ( self._metadata['capacity'], CORE_LOCI )
			
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
				array_path, dtype='int32', mode='w+', shape=(2000, CORE_LOCI) )

			self._allele_calls.flush()

			self._metadata['size'] = 0
			self._metadata['capacity'] = 2000

			self._length = 0
			self._capacity = 2000
			self._last_index = 0

	def resize(self):

		Logger.current._log('Resizing...')
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
		shape = (newsize, CORE_LOCI)
		
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

		pass

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
		return self._length

	def __iter__(self):
		for key, value in self._keys_to_index.items():
			yield key, self._allele_calls[value]

	def __getitem__(self, key):

		if key in self._keys_to_index:
			return self._allele_calls[ self._keys_to_index[key] ]

		else:
			return None

	def __delitem__(self, key):
		if key in self._keys_to_index:
			self.invalidate(key)
			del self._keys_to_index[key]
			self._length -= 1
		else:
			Logger.current._log('Key: {} does not exist'.format(key))

	def __contains__(self, key):
		return key in self._keys_to_index

if __name__ == '__main__':

	try:

		logging = Logger()

		metadata_path = os.path.join(os.getcwd(), 'lmo_allele_calls_all_updated.csv')

		allele_calls_dir = os.path.join(os.getcwd(), 'allele_calls', 'current')

		if not os.path.exists(allele_calls_dir):
			os.makedirs(allele_calls_dir)

		keys = set()
		all_calls = AlleleCalls(allele_calls_dir)
		all_calls.load()

		with open(metadata_path, 'rb' ) as f:

			reader = csv.reader(f)
			header = list(map(str.lower, next(reader)))

			for row in tqdm(reader):

				allele_calls = []

				for h,f in zip(header, row):
					if len(allele_calls) == 1748:
						break

					if h.startswith(_ORGANISM):

						if f == '' or f == '?':
							allele_calls.append(0)
						else:
							allele_calls.append(int(f))

				this_calls = np.asarray( allele_calls, dtype='int32' )
				
				# key = row[0].decode('cp1252').encode('utf8')
				key = unicode(row[0], errors='ignore')

				key = key.replace( u"\u00A0", ' ' )
				keys.add(key)

				inserted = all_calls.set( key, this_calls )

				# # if not inserted:
				# # 	print('Failed: {}'.format( key ))

				# if inserted:

				# 	print( key )
				# 	raise ValueError('This shouldnt have happened')

				# else:

				# 	if not np.array_equal( this_calls, all_calls[key] ):
				# 		raise ValueError('not equal')

		# randKey = keys.pop()

		# allele_calls = all_calls[randKey]

		# print(allele_calls)

		# del all_calls[randKey]

		print(all_calls._allele_calls[4117])
	except:
		traceback.print_exc()

	finally:
		all_calls.save()

		with open('keys.json', 'w') as f:
			json.dump(list(keys), f)