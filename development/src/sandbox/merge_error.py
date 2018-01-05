#!/usr/bin/env python

import os
import sys
from itertools import combinations


test_list = [ 1, 2, 3 ]


# Prefer the early merges, if they are within acceptable ranges
for i in range(len(test_list), 0, -1):

    for combo in combinations(test_list, i):

        print(combo)

        # if err( combo ) < max_acceptable_err:
        #   MERGE( combo )
        # else:
        #   continue