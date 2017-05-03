#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for computing the transfer factor for data. 

@file:   tf.py
@author: Andreas Søgaard
@date:   24 March 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic
import sys, os, glob, json, inspect

# Get ROOT to stop hogging the command-line options
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Scientific
try:
    
    # Numpy
    import numpy as np
    from root_numpy import *
    
except ImportError:
    print "WARNING: One or more scientific python packages were not found. If you're in lxplus, try running:"
    print "         $ source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh"
    sys.exit()
    pass

# Local include(s)
try:
    import transferfactor as tf
    from transferfactor.utils import make_directories
    from rootplotting import ap
    from rootplotting.tools import *
except ImportError:
    print "WARNING: This script uses the 'transferfactor' and 'rootplotting' packages. Clone them as e.g.:"
    print "         $ git clone git@github.com:asogaard/transferfactor.git"
    print "         $ git clone git@github.com:asogaard/rootplotting.git"
    sys.exit()
    pass


# Command-line arguments parser
import argparse

parser = argparse.ArgumentParser(description='Produce TF-estimated background files.')

#parser.add_argument('files', metavar='files', type=str, nargs='+',
#                    help='Files to combine')
#parser.add_argument('--DSID', dest='DSID', type=int,
#                    required=True,
#                    help='Pseudo dataset ID assigned to output')
parser.add_argument('--mass', dest='mass', type=int,
                    required=True,
                    help='Center of excluded mass window')
parser.add_argument('--window', dest='window', type=float,
                    default=0.2,
                    help='Relative width of excluded mass window (default: 20%%)')
parser.add_argument('--tau21DDT', dest='tau21DDT', type=float,
                    default=0.5,
                    help='Value of tau21DDT cut (default: 0.5)')
#parser.add_argument('--gamma', dest='gamma', type=float,
#                    default=1E-05,
#                    help='Kernel regression length scale (default: 1E-05)')
#parser.add_argument('--xsec_file', dest='xsec_file',
#                    default='../AnalysisTools/share/sampleInfo.csv',
#                    help='Path to cross-sections file (default: ../AnalysisTools/share/sampleInfo.csv)')
parser.add_argument('--show', dest='show',  action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save',  action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')
#parser.add_argument('--spline', dest='spline',  action='store_const',
#                    const=True, default=False,
#                    help='Use spline regression (default: False)')
#parser.add_argument('--gp', dest='gp',  action='store_const',
#                    const=True, default=False,
#                    help='Use gp regression (default: False)')
#parser.add_argument('--partialbins', dest='partialbins',  action='store_const',
#                    const=True, default=False,
#                    help='Use partially filled bins (default: False)')
#parser.add_argument('--emptybins', dest='emptybins',  action='store_const',
#                    const=True, default=False,
#                    help='Use empty bins (default: False)')
parser.add_argument('--subtractWZdata', dest='subtractWZdata',  action='store_const',
                    const=True, default=False,
                    help='Subtract estimated W/Z fail component from data (default: False)')


# Main function.
def main ():

    # Check(s)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Parse command-line arguments
    args = parser.parse_args()

    DSID = int("1%02d%03d" % (args.window * 100, args.mass))
    

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Load data
    files_data = glob.glob(tf.config['base_path'] + 'objdef_data_*.root')
    files_WZ   = glob.glob(tf.config['base_path'] + 'objdef_MC_30543*.root') + \
                 glob.glob(tf.config['base_path'] + 'objdef_MC_30544*.root')    

    if len(files_data) == 0 or len(files_WZ) == 0:
        warning("No files found.")
        return

    data      = loadData(files_data, tf.config['tree']) 
    WZ        = loadData(files_WZ,   tf.config['tree']) 
    info_data = loadData(files_data, tf.config['outputtree'], stop=1)
    info_WZ   = loadData(files_WZ,   tf.config['outputtree'], stop=1)

    
    # Rename variables
    for arr in [data, WZ]:
        arr.dtype.names = [name.replace(tf.config['prefix'], '') for name in arr.dtype.names] # @TODO: include in 'loadData(...)'?
        pass
    
    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])

    # Append new DSID field # @TODO: Make more elegant?
    WZ = append_fields(WZ, 'DSID', np.zeros((WZ.size,)), dtypes=int)
    for idx in info_WZ['id']:    
        msk = (WZ['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
        tmp_DSID = info_WZ['DSID'][idx]  # Get DSID for this file
        WZ['weight'][msk] *= xsec[tmp_DSID] # Scale by cross section x filter eff. for this DSID
        WZ['DSID']  [msk] = tmp_DSID        # Store DSID
        pass
    # @TODO: k-factors?
    WZ['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity

    # Check output.
    if data.size == 0:
        warning("No data was loaded.")
        return 

    # Compute new variables
    data = append_fields(data, 'logpt', np.log(data['pt']))
    WZ   = append_fields(WZ,   'logpt', np.log(WZ  ['pt']))
    

    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Pass/fail masks
    msk_data_pass = tf.config['pass'](data)
    msk_data_fail = ~msk_data_pass
    msk_WZ_pass   = tf.config['pass'](WZ)
    msk_WZ_fail   = ~msk_WZ_pass

    # Transfer factor calculator instance
    calc = tf.calculator(data=data, config=tf.config) # Using default configuration
    calc.mass   = args.mass
    calc.window = args.window
    # ... calc.partialbins, calc.emptybins, ...
    calc.fit() # ...(theta=0.5)
    w_nom     = calc.weights(data[msk_data_fail])
    w_up      = calc.weights(data[msk_data_fail], shift=+1)
    w_down    = calc.weights(data[msk_data_fail], shift=-1)
    w_WZ_nom  = calc.weights(WZ  [msk_WZ_fail])
    w_WZ_up   = calc.weights(WZ  [msk_WZ_fail], shift=+1)
    w_WZ_down = calc.weights(WZ  [msk_WZ_fail], shift=-1)
    if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/tf_')
    

    # Computing data-driven background estimate
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    if not os.path.exists('output/'):
        os.makedirs('output/')
        pass

    # Write TF-scaled failing data to file
    output = ROOT.TFile('output/objdef_TF_{DSID:6d}.root'.format(DSID=DSID), 'RECREATE')
    
    for shift, w, w_WZ in zip([0, 1, -1], [w_nom, w_up, w_down], [w_WZ_nom, w_WZ_up, w_WZ_down]):
        
        # -- Get branch name for current variation
        var_name = 'Nominal' if shift == 0 else ('TF_UP' if shift == 1 else 'TF_DOWN')

        # -- Prepare mass- and weight vectors
        vector_m = data['m']     [msk_data_fail]
        vector_w = data['weight'][msk_data_fail] * w
        if args.subtractWZdata:
            print "  Subtracting TF-scaled W/Z MC from background estimate"
            vector_m = np.concatenate((vector_m,   WZ['m']     [msk_WZ_fail]))
            vector_w = np.concatenate((vector_w, - WZ['weight'][msk_WZ_fail] * w_WZ))
            pass
        
        # -- Prepare DISD and isMC vectors
        vector_DSID = np.ones_like(vector_w) * DSID
        vector_isMC = np.ones_like(vector_w).astype(bool)
        
        array1 = np.array(zip(vector_m, vector_w),
                          dtype = [(tf.config['prefix'] + 'm', np.float64),
                                   ('weight',                  np.float64)])
        
        array2 = np.array(zip(vector_DSID, vector_isMC),
                          dtype = [('DSID', np.uint32),
                                   ('isMC', np.bool_)])
        
        # Mass and weight branch
        print "  Writing arrays to file " + ("(%s)" % ("up" if shift == 1 else "down") if shift != 0 else "")
        treename1 = tf.config['tree'].replace('NumLargeRadiusJets', 'Jet_tau21DDT').replace('Nominal', var_name)
        make_directories('/'.join(treename1.split('/')[:-1]), fromDir=output)
        tree1 = array2tree(array1, name=treename1.split('/')[-1])
        
        # outputTree
        treename2 = tf.config['outputtree'].replace('Nominal', var_name)
        make_directories('/'.join(treename2.split('/')[:-1]), fromDir=output)
        tree2 = array2tree(array2, name=treename2.split('/')[-1])

        output.Write()        
        pass

    output.Close()
        

    # Write TF-scaled failing W/Z MC to file
    output = ROOT.TFile('output/objdef_TF_{DSID:6d}_WZfail.root'.format(DSID=DSID), 'RECREATE')
    
    for shift, w_WZ in zip([0, 1, -1], [w_WZ_nom, w_WZ_up, w_WZ_down]):
               
        # -- Get branch name for current variation
        var_name = 'Nominal' if shift == 0 else ('TF_UP' if shift == 1 else 'TF_DOWN')

        # -- Prepare mass- and weight vectors
        vector_m = WZ['m']     [msk_WZ_fail]
        vector_w = WZ['weight'][msk_WZ_fail] * w_WZ
        
        # -- Prepare DISD and isMC vectors
        vector_DSID = np.ones_like(vector_w) * DSID
        vector_isMC = np.ones_like(vector_w).astype(bool)
        
        array1 = np.array(zip(vector_m, vector_w),
                          dtype = [(tf.config['prefix'] + 'm', np.float64),
                                   ('weight',                  np.float64)])
        
        array2 = np.array(zip(vector_DSID, vector_isMC),
                          dtype = [('DSID', np.uint32),
                                   ('isMC', np.bool_)])
        
        # Mass and weight branch
        print "  Writing arrays to file " + ("(%s)" % ("up" if shift == 1 else "down") if shift != 0 else "")
        treename1 = tf.config['tree'].replace('NumLargeRadiusJets', 'Jet_tau21DDT').replace('Nominal', var_name)
        make_directories('/'.join(treename1.split('/')[:-1]), fromDir=output)
        tree1 = array2tree(array1, name=treename1.split('/')[-1])
        
        # outputTree
        treename2 = tf.config['outputtree'].replace('Nominal', var_name)
        make_directories('/'.join(treename2.split('/')[:-1]), fromDir=output)
        tree2 = array2tree(array2, name=treename2.split('/')[-1])

        output.Write()        
        pass

    output.Close()
        
    # Save configuration
    if not os.path.exists('logs/'):
        os.makedirs('logs/')
        pass
    # -- Turn numpy arrays into lists, in order to make them JSON serializable
    def make_serializable (iterable):
        """ Turn numpy arrays into lists, in order to make them JSON serializable. """
        result = iterable
        if   type(result) == dict:
            for key, element in result.iteritems():
                result[key] = make_serializable(element)
                pass
        elif type(result) == list:
            for idx, element in enumerate(result):
                result[idx] = make_serializable(element)
                pass
        elif type(result).__module__.startswith(np.__name__):
            result = result.tolist()
        elif type(result).__name__ == 'function':
            result = inspect.getsource(result)
            pass
        return result

    cfg = make_serializable(tf.config)

    json.dump([cfg, vars(args)], open('logs/tf_config_%d.log' % DSID, 'w'))        

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
