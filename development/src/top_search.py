import os
import sys
import csv
import numpy as np
import csv
import json
import random
import time

from collections import defaultdict, namedtuple, OrderedDict
from itertools import combinations
from tqdm import *

import_file = os.path.join(os.getcwd(), 'results/clustering_validation/' 
    'database_09-13-17@15-11-39.csv')

wgst_database = defaultdict(list)
Entry = namedtuple( 'Entry', ['key', 'wgst', 'allelecalls'], verbose=True, rename=False )
thresholds = ( 4.050815987, 2.914238135, 2.046511241, 1.068159867, 0.041631973 )


def sim(p1, p2):

    present = np.multiply( p1>0, p2>0)
    ncommon = np.sum(present)
    # toReturn = np.empty(len(present))
    nsame = np.sum( p1[present] == p2[present] )
    # for i, val in enumerate( present ):

    #     if not val:
    #         toReturn[i] = -1
    #         continue

    #     if p1[i] == p2[i]:
    #         toReturn[i] = p1[i]
    #         nsame +=1
    #     else:
    #         toReturn[i] = -1

    return (100 * (float( ncommon ) - float(nsame) ) / ( float(ncommon) ))

with open(import_file, 'r') as f:

    reader = csv.reader(f)
    header = next(reader)
    header = list(map(str.lower, header))
    
    for row in tqdm(reader):

        al = np.asarray( row[2:], dtype=int)
        if (np.sum( al>0 ) / len(al)) < 0.95:
            continue 

        e = Entry( row[0], row[1], np.asarray( row[2:], dtype=int ) )

        wgst_database[e.wgst.split('.')[0]].append( e )

print(len(wgst_database))

overlap_test = {}
test_isolates = {}
failed = 0

for key in tqdm( wgst_database ):

    if len( wgst_database[key] ) == 1:
        continue

    overlap_test[key] = [ set() for _ in range(1748) ]

    random.shuffle( wgst_database[key] )
    test_isolates[key] = wgst_database[key].pop()

    for entry in wgst_database[key]:
        for i, val in enumerate(entry.allelecalls):

            if val > 0:
                overlap_test[key][i].add( int(val) )
    
    # for entry, other in tqdm( combinations(wgst_database[key], 2) ):

    #     # distance  = sim( entry.allelecalls, other.allelecalls )

    #     for i, val in enumerate( same ):
    #         if val >= 0:
    #             overlap_test[ key ][i][0].add( int(entry.allelecalls[i]) )
    #             overlap_test[ key ][i][1] += 1

testing = {}
for cluster, isolate in test_isolates.items():
    
    start = time.clock()

    foundclusters = []
    for key, identity in overlap_test.items():

        similarityHere = 0.
        totalPresent = 0.

        for isolate_call, identity_calls in zip( isolate.allelecalls, identity ):

            # similarityHere += isolate_call in identity_calls

            # totalPresent += isolate_call > 0 and len(identity_calls) > 0

            if isolate_call != 0 and len(identity_calls) > 0:

                if isolate_call in identity_calls:
                    similarityHere +=1

                totalPresent +=1

        # foundclusters.append( [key, similarityHere, totalPresent, similarityHere/totalPresent, len( wgst_database[key] )])
    print( time.clock() - start )

    # foundclusters.sort( key=lambda x: -x[3] )
    # assert cluster == foundclusters[0][0]

#     testing[cluster] = { "foundclusters": foundclusters, "query": isolate.allelecalls.tolist() }


# with open('query.json', 'w') as f:
#     json.dump( testing, f)





# for key in overlap_test:
#     toConvert = overlap_test[key]

#     for i, val in enumerate(toConvert):
#         overlap_test[key][i] = [list(val[0]), val[1]]

# with open('results.json', 'w') as f:
#     json.dump(overlap_test, f)

# count_dict = {}

# for key, val in wgst_database.items():
#     count_dict[key] = len(val)

# with open( 'counts.json', 'w') as f:
#     json.dump( count_dict, f )
