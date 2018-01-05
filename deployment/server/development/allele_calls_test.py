import os
import sys
import numpy as np
import json
import traceback

class Logger(object):

    current = None

    def _log(self, oStr):

        try:

            print(oStr)

        except:

            try:
                debug_str = str(oStr)
                print(debug_str)

            except:
                traceback.print_exc()

class AlleleCalls(object):

    def __init__(self, dir_path, loci_size):
        self._dir_path = dir_path
        self._entrymapping = {}
        self._length = -1
        self._array = None
        self._loci_size = loci_size

    def load(self):

        entry_mapping_path = os.path.join(self._dir_path, 'entry_mapping.json')
        
        if not os.path.exists(entry_mapping_path):
            Logger.current._log('Could not find file')
        else:
            with open(entry_mapping_path, 'r') as f:
                self._entrymapping = json.load(f)
                self._length = len(self._entrymapping)


        if self._length < 0:
            self._array = np.memmap( )

        self._array = np.memmap