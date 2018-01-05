import os
import sys
from itertools import combinations


test_list = [ 1, 2 ,3 ,4 ,5 ]


for i in range(len(test_list), -1, -1):

    for combo in combinations(test_list, i):

        print(combo)