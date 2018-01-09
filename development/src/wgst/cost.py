

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

