import os
import csv
import numpy as np
from tqdm import *
import pdb
from collections import defaultdict
from json import dumps, loads, JSONEncoder, JSONDecoder
import pickle

class PythonObjectEncoder(JSONEncoder):
    def default(self, obj):

        if isinstance(obj, set):
            return list(obj)
        return JSONEncoder.default(self, obj)

def as_python_object(dct):
    if '_python_object' in dct:
        return pickle.loads(str(dct['_python_object']))
    return dct


dir_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
allele_path = os.path.join(dir_path, 'allele_calls/lmo_allele_calls_all_updated.csv')


database = {}
rm_keys = set()
with open(allele_path, 'r', encoding='cp1252') as f:
    reader = csv.reader( f )
    header = tuple(map(str.lower, next(reader)))

    for row in tqdm( reader, desc='Loading' ):

        flds_data = {}
        allele_calls = []
        for h, f in zip(header, row):

            if h.startswith('lmo'):

                if f == '' or f=='?':
                    allele_calls.append(0)
                else:
                    try:

                        allele_calls.append(int(f))
                    except:
                        rm_keys.add(row[0])

                        break

            else:

                flds_data[h] = f


        flds_data['allele_calls'] = allele_calls

        database[row[0]] = flds_data


for rm in rm_keys:
    del database[rm]


wgst_path = os.path.join( dir_path, 'results/database/database_11-02-17@15-44-42.csv')


with open(wgst_path, 'r', encoding='cp1252') as f:

    reader = csv.reader(f)
    header = tuple(map(str.lower, next(reader)))

    for row in tqdm(reader, desc='Loading wgst'):

        flds_data = { h:f for h,f in zip(header, row) }

        if database.get(row[0], False):

            database[row[0]]['wgst'] = flds_data['wgst']

pfge_path = os.path.join(dir_path, 'allele_calls/lmo_pfge_export.csv')

with open(pfge_path, 'r', encoding='cp1252') as f:

    reader = csv.reader(f)
    header = tuple(map(str.lower, next(reader)))

    for row in tqdm(reader, desc='Loading pfge'):

        if database.get(row[0], False):
            database[row[0]].update(((h,f) for h,f in zip(header, row)))

wgst_list = defaultdict(lambda: defaultdict(set))

for entry in database.values():

    if entry.get('wgst', False):

        wgst_list[entry['wgst']]['pfge-apai-pattern'].add(entry['pfge-apai-pattern'])
        wgst_list[entry['wgst']]['pfge-asci-pattern'].add(entry['pfge-asci-pattern'])


out_path = os.path.join(dir_path, 'results/wgst_pfge.json')
with open(out_path, 'w') as f:
    f.write(dumps( wgst_list, cls=PythonObjectEncoder))


pfge_list = defaultdict(set)

for entry in database.values():

    if entry.get('pfge-apai-pattern', False) and entry.get('pfge-asci-pattern', False) and entry.get('wgst', False):

        pfge_str = '{}|{}'.format( entry['pfge-apai-pattern'], entry['pfge-asci-pattern'] )

        pfge_list[pfge_str].add( entry['wgst'] )


out_path = os.path.join(dir_path, 'results/pfge_wgst.json')

with open(out_path, 'w') as f:
    f.write(dumps( pfge_list, cls=PythonObjectEncoder))


# pdb.set_trace()





