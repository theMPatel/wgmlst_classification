###############################################################################
#
# Initializer for pipeline. 
# AUTHOR: MILAN PATEL
# CONTACT: mpatel5@cdc.gov
#
###############################################################################

import os
import sys
import argparse
from ast import literal_eval
import json

# This function gets the command line
def ParseAndValidateCommandLine():

    def Usage():
        usagetxt = \
        """
        This pipeline is used to either calculate or optimize
        nomenclature based on your BioNumerics (cg/wg)MLST schemes


        Required:

        [--t]:  The location to your thresholds file,
                it will read the first line only

        [--f]:  The location of your flds/allele calls
                data

        [--o]:  Your output directory, pipeline will
                automatically create the necessary
                folders inside

        [--g]:  Organism name e.g. listeria; salmonella; ecoli

        [--m]:  Minimum allele presence threshold as a
                decimal value i.e.: 0.95

        Optional:
        
        [--v]:  Views file, stores the loci to look at
                for the core scheme. See scheme.json file
                ATM your only options are core/whole

        [--d]:  This flag will initialize optimization.
                Bare bones requirement: 4 threads
                For this to be useful: >=20 threads
                
                This will automatically calculate
                the distance matrix if it hasn't
                been already

        [--dm]: In case you want to only calculate
                the distance matrix and optimize 
                at another time. If you need to re-
                calculate the DM for any reason,
                use this flag. You will then need to
                launch another optimize job without this
                flag
        
        [--a]:  In case you want to manually define
                your seed database and 'new' isolates to
                add to your database. Point --f to your 
                seed file, and --a to your 'new'
                isolates file

        [--c]:  To be paired with --a; Will save a copy
                of the database with only your 'new' isolates
                inside so you can check SN/SP or cluster 
                search on new data
        """
        return usagetxt

    parser = argparse.ArgumentParser( description = 'Nomenclature calculation '
                    'or validation, default is to calculate', usage=Usage() )

    parser.add_argument('--d', '--validate', help='Provide num cores for threshold, '
                    'optimization. Minimum: 4', type=int, default=None)
    
    parser.add_argument( '--a', '--add', help='Provide an external adding file.', type=str,
        default=None, action='append')
    
    parser.add_argument( '--o', '--outdir', help='Output directory', type=str, 
        default='./results' )

    parser.add_argument( '--c', '--clustersearch', help='Search for clusters', type=bool,
        default=False )

    parser.add_argument( '--dm', '--dmonly', help='Calculate only the distance matrix', type=str,
        default=False )

    parser.add_argument( '--s', '--scheme', help='core or whole genome', required=True, 
                    type=str )

    parser.add_argument( '--t', '--thresholds', help='Path to thresholds file', required=True,
                    type=str )

    parser.add_argument( '--f', '--fields', help='Path to db csv file', required=True,
                    type=str )

    parser.add_argument( '--v', '--views', help='Path to json file for character views', 
        required=False, type=str )

    parser.add_argument( '--m', '--minPres', help='Minimum presence threshold', required=False,
                    type=float, default=-1.0 )

    parser.add_argument( '--g', '--organism', help='Organism this is for', required=True,
        type=str)

    parser = parser.parse_args()

    # Thresholds file
    if not os.path.isfile( parser.t ):
        raise RuntimeError( 'Please provide a valid thresholds file' )

    # Validate schemes
    validViews = ( 'core', 'whole' )
    if parser.s not in validViews:
        raise ValueError( 'Please provide a valid scheme type' )

    # Make sure that the metadata file is valid
    if not os.path.isfile( parser.f ) or not parser.f.endswith('.csv'):
        raise RuntimeError( 'Please provide a valid csv db' )
    
    parser.f = {'train': parser.f}

    # Import configuration file
    _CONFIG_PATH = os.path.dirname(os.path.realpath(__file__)) +'/config.json'

    if not os.path.exists(_CONFIG_PATH):
        raise RuntimeError('Missing configuration file')

    with open(_CONFIG_PATH, 'r') as f:
        config = json.load(f)

    if parser.g in config['organisms']:
        if parser.s in config['organisms'][parser.g]['schemes']:
            parser.locisize = config['organisms'][parser.g]['schemes'][parser.s]['loci']
        else:
            parser.locisize = -1
            parser.v = None

        if parser.m < 0:
            parser.m = config['organisms'][parser.g]['schemes'][parser.s]['presence']

        parser.g = config['organisms'][parser.g]['shorthand']
    else:
        raise RuntimeError('Organism not in config file'
            ' bare minimum need organism shorthand')

    # Checks the min presence, must come after the config file checking
    if not 0. <= parser.m <= 1.:
        raise ValueError( 'Please provide a valid decimal presence value' )

    # If views file is provided, check that it's real
    if parser.v:
        if not os.path.isfile( parser.v ) or not parser.v.endswith( '.json' ):
            raise ValueError( 'Please provide the character view file in json format' )
    else:
        assert parser.s == 'whole'
        parser.v = None

    # External adding set, for test data
    if parser.a is not None:
        parser.f['test'] = []
        for file in parser.a:
            if not os.path.isfile( file ):
                raise ValueError( 'Please provide a valid adding file' )
            else:
                parser.f['test'].append( file )

    # Out directory
    if not os.path.isdir( parser.o ):
        os.mkdir( parser.o )

    if parser.d is not None:

        if ( parser.d - 4 ) < 0:
            raise ValueError('You will need a minimum of 4 threads for this')
        
        else:
            thresholds = []
            with open( parser.t, 'r' ) as f:
                data = f.readline().strip().replace(' ', '').split('|')

            for threshRange in data:
                
                try:
                    thresholds.append( literal_eval( threshRange ) )
                except Exception as e:
                    print('Make sure you have formated the thresholds file properly')
                    print( e )
                    sys.exit(0)

            for value in thresholds[:-1]:
                
                if not isinstance(value, tuple) and not isinstance( 
                    value, list):

                    raise ValueError('Make sure you have provided optimization thresholds as list/tuple items')

            if not isinstance( thresholds[-1], float ) and not isinstance( thresholds[-1], int ):
                raise TypeError('Provide a valid step value as last value in threshold list not in list')
    elif parser.dm:
        thresholds = None
    else:
        with open( parser.t, 'r' ) as f:
            thresholds =  f.readline()
        try:
            thresholds = thresholds.strip().replace(' ', '').split(',')
            thresholds = list( map( float, thresholds ) )
        except Exception as e:
            print( e )
            print('Make sure that your thresholds are properly formed')
            sys.exit(0)

    for ft, fp in parser.f.items():
        if isinstance(fp, list):
            for i in range(len(fp)):
                fp[i] = os.path.abspath(fp[i])
        else:

            parser.f[ft] = os.path.abspath( parser.f[ft] )
    
    if parser.v:
        parser.v = os.path.abspath( parser.v )
    
    parser.o = os.path.abspath( parser.o )
    
    if parser.a:
        parser.a = os.path.abspath( parser.a )

    return parser, thresholds

if __name__ == '__main__':

    args, thresholds = ParseAndValidateCommandLine()

    if args.dm:
        from wgst import overlap_v2 as overlap
        distance_matrix = overlap.Main({
            'organism':         args.g,
            'fields_path':      args.f,
            'minPres':          args.m,
            'scheme':           args.s,
            'schemespath':      args.v,
            'recalculate':      args.dm,
            'outdir':           args.o,
            'locisize':         args.locisize
            })

    elif args.d:
        from wgst import overlap_v2 as overlap
        distance_matrix = overlap.Main({
            'organism':         args.g,
            'fields_path':      args.f,
            'minPres':          args.m,
            'scheme':           args.s,
            'schemespath':      args.v,
            'recalculate':      args.dm,
            'outdir':           args.o
            })

        from wgst import multi_v2 as multi
        multi.Main({
            'cores':        args.d,
            'matrix':       distance_matrix,
            'scheme':       args.s,
            'thresholds':   thresholds,
            'fields_path':  args.f,
            'outdir':       args.o,
            'organism':     args.g
            })

    else:
        from wgst import single_v2
        single.Main({
            'fields_path':      args.f,
            'minPres':          args.m,
            'scheme':           args.s,
            'organism':         args.g,
            'thresholds':       thresholds,
            'views_path':       args.v,
            'outdir':           args.o,
            'clustersearch':    args.c,
            'locisize':         args.locisize
            })