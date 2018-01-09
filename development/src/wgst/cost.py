##############################################################################
# Date: 1/9/2018
# Author: Milan Patel
# Cost functions
##############################################################################

from itertools import combinations

__all__ = ['mean_error', 'mean_inter_error', 'ABSOLUTE_MARGIN_ERROR', 'SQUARED_ERROR']

# This is for Listeria, you will need to change this at runtime
ABSOLUTE_MARGIN_ERROR = 0.2860411899
SQUARED_ERROR = -1

def mean_error(unnamed, key_list, d_matrix, thresholds, level, 
    biased=True, squared=False):

    results = []
    for key in key_list:
        err = thresholds[level-1]-d_matrix[unnamed,key]

        if err >= 0.:

            if biased:
                continue
            else:
                results.append(0.)

        elif squared:
            results.append(err**2)

        else:
            results.append(err)

    if len(results):
        return abs(sum(results))/float(len(results))

    else:
        return 0.

def mean_inter_error(unnamed, key_list, d_matrix, thresholds, level,
    biased=True, squared=False):

    
    results = []

    for entry, other in combinations(key_list, 2):

        err = thresholds[level-1]-d_matrix[entry,other]

        if err >= 0.:

            if biased:
                continue
            else:
                results.append(0.)

        elif squared:
            results.append(err**2)

        else:
            results.append(err)

    for key in key_list:

        err = thresholds[level-1]-d_matrix[unnamed,key]

        if err >= 0.:

            if biased:
                continue
            else:
                results.append(0.)

        elif squared:
            results.append(err**2)

        else:
            results.append(err)

    if len(results):
        return abs(sum(results))/float(len(results))

    else:
        return 0.