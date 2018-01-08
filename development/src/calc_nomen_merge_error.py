#!/usr/bin/env python
import os
import sys
import base64
import gzip
import json
import shutil
import traceback
import numpy as np
import csv


from copy import deepcopy
from tqdm import *
from nomenclature.bn_export_parser import Database, Entry
from nomenclature.tree import Tree, Node, NamedNode
from nomenclature.distance_matrix import DistanceMatrix as DM

from datetime import datetime, date
from collections import defaultdict


def get_distance( p1, p2 ):
    common = np.multiply(p1>0, p2>0)
    nCommon = np.sum(common)
    nSame = np.sum(p1[common]==p2[common])
    if nCommon:
        return 100.0 * (float(nCommon) - float(nSame)) / float(nCommon) #nCommon-nSame
    else:
        return 100.

def calculate_name(named, tree, unNamedEntry, distances, thresholds, d_matrix ):

    margin_error = 0.2860411899

    def valid_merge(level, node_list, mnode, unnamed):

        to_check = []

        for node in node_list:
            for named in mnode.GetChild(node).DFSNamed():
                to_check.extend(key for key in named.GetChildrenKeys())

        errs = [ thresholds[level-1] - d_matrix[unnamed,other] for other in to_check]
        errs = filter(lambda x: x<0, errs)

        if not len(errs):
            cost = 0.
        else:
            cost = abs(sum(errs)) / float(len(errs))

        if cost < margin_error:
            return True
        else:
            return False

    def calc_merge_cost(node, level, unnamed):

        keys_to_check = [y for x in node.DFSNamed() for y in x.GetChildrenKeys() ]
        errs = []

        for other in keys_to_check:
            t = thresholds[level-1]
            dist = d_matrix[ unnamed , other]

            errs.append( t - dist)

        # errs = [ thresholds[level-1]-d_matrix[unnamed,other] for other in keys_to_check ]

        return errs



    # Make sure thresholds are sorted biggest first
    thresholds.sort(key=lambda x: -x)

    # Let's get the head of the tree
    headNode = tree.Tree()

    # UnNamed should have no name
    patternName = tree.GetName( unNamedEntry )
    assert len( patternName ) == 0

    # Will serve as the node place holder as we drill down
    # into the proper home for a our unNamedEntry
    currentNode = headNode

    # Create a local copy of the named list
    namedEntries = list( named )
    for level, threshold in enumerate( thresholds, 1 ):

        # Create the clusters at this level:
        closestClusters = set()
        toRemove = []
        for i, dist in enumerate( distances ):
            
            if dist <= threshold:
                keyToGet = namedEntries[i]
                partial = tuple( tree.GetPart( keyToGet, level ) )
                closestClusters.add( partial )

            else:
                toRemove.append( i )

        # If something isn't similar at this level
        # There's no way it will be similar at the lower
        # levels, let's remove those
        for index in reversed(toRemove):
            del namedEntries[ index ]
            del distances[ index ]

        # Assigned to only one cluster
        if len( closestClusters ) == 1:

            pattern = closestClusters.pop()
            err = calc_merge_cost(
                currentNode.GetChild(pattern[-1]),
                level,
                unNamedEntry)
            
            cost_list = filter(lambda x: x<0, err)

            if not len(cost_list):
                cost = 0.
            else:
                cost = abs(sum(cost_list))/float(len(cost_list))

            if cost < margin_error:
                patternName.append( pattern[-1] )
                currentNode = headNode.Traverse( pattern )
            else:
                currentNode = currentNode.NewChildNode()
                patternName.append( currentNode.ID() )

                distances = []
                namedEntries = []

        # New cluster!
        elif len( closestClusters ) == 0:
            currentNode = currentNode.NewChildNode()
            patternName.append( currentNode.ID() )

        # Merge time:
        else:
            if valid_merge(level, [c[-1] for c in closestClusters], currentNode, unNamedEntry):
                toMerge =[ c[-1] for c in closestClusters ]
                currentNode = currentNode.MergeNodes( toMerge )
                patternName.append( currentNode.ID() )
            else:

                best_cluster = []

                for cnode in [ c[-1] for c in closestClusters ]:

                    child = currentNode.GetChild(cnode)
                    if child is None:
                        raise RuntimeError('NoneType child!')

                    err = calc_merge_cost(child, level, unNamedEntry)
                    cost_list = filter(lambda x: x<0, err)
                    
                    if not len(cost_list):
                        cost = 0.
                    else:
                        cost = abs(sum(cost_list))/float(len(cost_list))

                    best_cluster.append( (cnode, cost) )

                best_cluster = filter(lambda x: x[1]<margin_error, best_cluster)

                if not len(best_cluster):
                    currentNode = currentNode.NewChildNode()
                    patternName.append( currentNode.ID() )

                    distances = []
                    namedEntries = []

                else:
                    best_cluster.sort(key=lambda x: x[1])
                    patternName.append( best_cluster[0][0] )
                    currentNode = headNode.Traverse( patternName )

                    remove_others = []

                    for i, key in enumerate(namedEntries):

                        if tree.GetPart(key, level)[-1] != patternName[-1]:
                            remove_others.append(i)

                    remove_others.sort()

                    for rm in reversed( remove_others ):
                        del namedEntries[rm]
                        del distances[rm]

    assert isinstance( currentNode, NamedNode )

    currentNode.AddChild( unNamedEntry )

    assert currentNode.GetChild( unNamedEntry ) == patternName

    return tree

def main():
    data_path = os.path.expanduser('~/Documents/Work/wgst/stability/allele_calls/lmo_allele_calls_all_updated_partial.csv')
    database = Database('lmo', data_path )
    tree = Tree(5)
    thresholds = [4.050815987,2.914238135,2.046511241,1.068159867,0.041631973]
    min_pres = 0.50
    database.qc(min_pres)

    distance_holder = DM()
    
    named = []
    i = 0
    for entry in tqdm( database, desc='Clustering' ):

        dists = []
        for other in named:
            dists.append( get_distance( entry.allele_calls, database[other].allele_calls ) )


        distance_holder.add( entry.key, dists )
        calculate_name(named, tree, entry.key, deepcopy(dists), thresholds, distance_holder)

        if tree.HasName(entry.key):
            named.append(entry.key)

        i+=1


    # print(distance_holder._distance_matrix)
    out_file = 'results_{}.csv'.format(datetime.now().strftime('%m-%d-%Y@%H-%M-%S'))

    with open(out_file, 'w') as f:
        writer = csv.writer(f)

        for entry in database:

            writer.writerow([entry.key, tree.GetStrName(entry.key)])

if __name__ =='__main__':
    main()








