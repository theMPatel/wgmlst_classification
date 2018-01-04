import os
import csv
import numpy as np
from tqdm import *

file = os.path.join(os.getcwd(),'ecoli_allele_calls_all.csv')

calls = {}

with open(file, 'r', encoding='cp1252') as f:

    reader = csv.reader(f)
    header = list(map(str.lower, next(reader)))


    for row in tqdm(reader, 'qcing'):
        allele_calls = []

        for h,f in zip(header, row):

            if h.startswith('ec'):

                if f == '' or f == '?':

                    allele_calls.append(0)

                else:
                    allele_calls.append(int(f))



        calls[row[0]] = allele_calls


def qc(calls, val):

    poorqc = []

    val = float(val)

    for entry in tqdm(calls, 'qcing'):

        total = 0

        for i in calls[entry]:
            if i >0:
                total += 1


        if 100.*( total / len(calls[entry]) ) < val:
            poorqc.append( [entry, total/len(calls[entry] ) ] )


    for result in poorqc:
        yield result


toCheck = (float("inf"),)

for val in toCheck:

    with open(str(val)+'.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(qc(calls, val))
