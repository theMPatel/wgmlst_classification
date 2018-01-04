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

import_file = os.path.join( os.getcwd(), 'allele_calls/salm_allele_calls_test.csv')
# import_file = os.path.join( os.getcwd(), 'allele_calls/lmo_allele_calls_all_updated.csv')
Entry = namedtuple( 'Entry', ['key', 'allelecalls'], verbose=True, rename=False )

def sim(p1, p2):

    present = np.multiply( p1>0, p2>0)
    ncommon = np.sum(present)
    nsame = np.sum( p1[present] == p2[present] )

    return (100 * (float( ncommon ) - float(nsame) ) / ( float(ncommon) ))

entryList = []
with open(import_file, 'r', encoding='cp1252') as f:

    reader = csv.reader(f)
    header = next(reader)
    header = list(map(str.lower, header))
    
    for row in tqdm(reader):


        al = list()

        for h, f in zip( header, row ):

            if h.startswith('saen'):
                if f == '':
                    al.append(0)
                else:
                    al.append(int(f))

        alleles = np.asarray( al, dtype=int)
        if (np.sum( alleles>0 ) / len(alleles)) < 0.28:
            continue

        # print(len(alleles))

        e = Entry( row[0], alleles )

        entryList.append( e )

random.shuffle( entryList )
tried = 0
averagetime = []


# for i in trange(1, len(entryList)+1):
#     start = time.clock()

#     # distances = np.empty(len(entryList[1:i]), dtype=float)
#     distances = []
#     for j in range(1,i):
#         distances.append( sim(entryList[0].allelecalls, entryList[j].allelecalls) )

#     end = time.clock() - start

#     print(end)

#     averagetime.append(end)


for i in range(100):

    start = time.clock()
    distances = []

    for j in range( 1, len(entryList)):
        distances.append( sim(entryList[0].allelecalls, entryList[j].allelecalls))


    end = time.clock() - start

    print(end)


with open('test_distance_{}.csv'.format(entryList[0].key), 'w') as f:

    writer = csv.writer(f)

    for i, val in enumerate(averagetime):
        writer.writerow([i+1, val])