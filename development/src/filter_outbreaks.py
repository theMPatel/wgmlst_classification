import os
import sys
import csv

file_path = os.path.join(os.getcwd(), 'ecoli_results/snsp_data/snsp_data_12-14-17@12-40-23.csv')
out_path = os.path.join( os.getcwd(), 'ecoli_results/filtered_outbreaks/snsp_data_12-14-17@15-34-00.csv')

filtered = {'1702CAEXH-1',
'1411CAEXH-1',
'1504MLEXH-1',
'1511MTEXH-1',
'1511MLEVC-1',
'1511MLEXW-1',
'1602CAEXH-1',
'1709CAEXH-1',
'0902COEXH-1',
'1710COEVC-1',
'1609MLEXH-1',
'0710NYEXH-1mlc',
'1004MIENM-1',
'0806MLEXH-3c',
'1608MTEXH-1',
'1607NHEXH-1',
'0612mlEXH-2',
'1506SCEXH-1',
'1509TXEXD-1',
'1707UTEXH-1',
'1505WAEXH-1',
'1709MLEXH-2',
'1601MLEXK-1'}

with open(file_path, 'r') as in_file, open(out_path, 'w') as out_file:
    in_reader = csv.reader(in_file)
    out_writer = csv.writer(out_file)

    for row in in_reader:

        if row[0] == 'KEY':
            out_writer.writerow(row)
            continue

        if row[1] and row[1] not in filtered:
            row[1] = ''
            out_writer.writerow(row)
        else:
            out_writer.writerow(row)





