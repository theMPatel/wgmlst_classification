import os
import sys
import csv
import numpy as np
import csv
import json

from collections import defaultdict, namedtuple
from itertools import combinations
from tqdm import *

import_file = os.path.join(os.getcwd(), 'results/clustering_validation/' 
    'database_09-13-17@15-11-39.csv')

wgst_database = defaultdict(list)
Entry = namedtuple( 'Entry', ['key', 'wgst', 'allelecalls'], verbose=True, rename=False )
thresholds = ( 4.050815987, 2.914238135, 2.046511241, 1.068159867, 0.041631973 )


def sim(p1, p2):

    return  100. - 100 * ( np.sum( p1 == p2 ) / len(p1) )

def simTri(p1, p2):

    common = np.multiply(p1>0, p2>0)
    nSame = np.sum( p1[common] == p2[common] )
    ncommon = np.sum(common)

    return (100 * (float( ncommon ) - float(nsame) ) / ( float(ncommon) ))

with open(import_file, 'r') as f:

    reader = csv.reader(f)
    header = next(reader)
    header = list(map(str.lower, header))
    for row in tqdm(reader):

        e = Entry( row[0], row[1], np.asarray( row[2:], dtype=int ) )

        wgst_database[e.wgst.split('.')[0]].append( e )


print(len(wgst_database))

overlap_test = defaultdict(list)
failed = 0
for key in tqdm( wgst_database ):
    
    for entry, other in tqdm( combinations(wgst_database[key], 2) ):
        distance = sim( entry.allelecalls, other.allelecalls )

        overlap_test[key].append(distance)

with open('newdist.json', 'w') as f:
    json.dump(overlap_test, f)

