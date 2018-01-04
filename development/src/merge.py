import os
import sys
import csv
from datetime import datetime
from tqdm import *


# database_path = 'ecoli_results/database/database_12-12-17@09-13-30.csv'
database_path = 'ecoli_results/filtered_outbreaks/snsp_data_12-14-17@15-34-00.csv'
ecoli_export_path = 'allele_calls/export_Ecoli_ST.csv'

database = {}
ref_header = ('key', 'outbreak', 'mlst_st', 'sourcetype', 'uploaddate', 'wgst')

with open(os.path.join(os.getcwd(), database_path), 'r', encoding='cp1252') as f:
    reader = csv.reader(f)
    header = list(map(str.lower, next(reader)))

    for row in tqdm( reader, desc ='loading' ):

        flds = {h:f for h,f in zip(header, row)}

        database[row[0]] = flds

with open(os.path.join(os.getcwd(), ecoli_export_path), 'r', encoding='cp1252') as f:

    reader = csv.reader(f)
    header = list(map(str.lower, next(reader)))

    for row in tqdm( reader, desc='loading sts'):

        if not database.get(row[0], False):
            continue

        for h,f in zip(header, row):

            if h == 'sourcetype':

                database[row[0]]['sourcetype'] = f

            elif h == 'outbreak':
                database[row[0]]['outbreak'] = f

            elif h =='mlst achtman st':
                database[row[0]]['mlst_st'] = f


out_path = 'ecoli_results/snsp_data/merged_{}.csv'
out_path = out_path.format( 
            datetime.now().strftime("%m-%d-%y@%H-%M-%S") )

with open(os.path.join(os.getcwd(), out_path), 'w') as f:

    writer = csv.writer(f)
    writer.writerow(ref_header)

    for key, value in database.items():

        writer.writerow([value[x] for x in ref_header])