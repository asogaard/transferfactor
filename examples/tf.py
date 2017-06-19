#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for computing the transfer factor for data. 

@file:   tf.py
@author: Andreas Søgaard
@date:   24 March 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys, os, glob, json

# Get ROOT to stop hogging the command-line options
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Scientific import(s)
try:
    # Numpy
    import numpy as np
    from root_numpy import *
except ImportError:
    print "WARNING: One or more scientific python packages were not found. If you're in lxplus, try running:"
    print "         $ source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh"
    sys.exit()
    pass

# Local import(s)
try:
    import transferfactor as tf
    from transferfactor.utils import make_directories, check_make_dir, make_serializable, get_signal_DSID
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

parser.add_argument('--mass', dest='mass', type=int,
                    required=True,
                    help='Center of excluded mass window')
parser.add_argument('--window', dest='window', type=float,
                    default=None,
                    help='Width of excluded mass window (default: None)')
parser.add_argument('--tau21DDT', dest='tau21DDT', type=float,
                    default=0.5,
                    help='Value of tau21DDT cut (default: 0.5)')
parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')
parser.add_argument('--subtractWZdata', dest='subtractWZdata', action='store_const',
                    const=True, default=False,
                    help='Subtract estimated W/Z fail component from data (default: False)')
parser.add_argument('--subtractWZMC', dest='subtractWZMC', action='store_const',
                    const=True, default=False,
                    help='Subtract W/Z MC from TF profile (default: False)')


# Main function.
def main ():

    # Parse command-line arguments
    args = parser.parse_args()

    DSID = int("1%02d%03d" % (0 if args.window is None else args.window * 100, args.mass))
    

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Get signal file
    sig_DSID = get_signal_DSID(args.mass)

    # Load data
    #files_data = glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root') 
    files_data = glob.glob(tf.config['base_path'] + 'objdef_data_*.root')
    files_WZ   = glob.glob(tf.config['base_path'] + 'objdef_MC_30543*.root') + \
                 glob.glob(tf.config['base_path'] + 'objdef_MC_30544*.root')    
    
    if sig_DSID is None:
        if args.mass < 100.:
            print "Assuming signal is W/Z"
            files_sig = files_WZ    
            files_WZ  = []
        else:
            files_sig = []
            pass
    else:
        sig_file = 'objdef_MC_{DSID:6d}.root'.format(DSID=sig_DSID)
        print "Using signal file: %s" % sig_file
        files_sig = [tf.config['base_path'] + sig_file]
        pass

    if len(files_data) == 0 or (sig_DSID and len(files_sig) == 0):
        warning("No files found.")
        return

    data      = loadData(files_data, tf.config['tree'], prefix=tf.config['prefix']) 
    signal    = loadData(files_sig,  tf.config['tree'], prefix=tf.config['prefix']) 
    WZ        = loadData(files_WZ,   tf.config['tree'], prefix=tf.config['prefix']) 
    info_data = loadData(files_data, tf.config['outputtree'], stop=1)
    info_sig  = loadData(files_sig,  tf.config['outputtree'], stop=1)
    info_WZ   = loadData(files_WZ,   tf.config['outputtree'], stop=1)
    
    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])

    # ----------------------------------------------------
    # Make more elegant!
    # ----------------------------------------------------
    # Append new DSID field # @TODO: Make more elegant?
    #for arr, info in zip([signal, WZ], [info_sig, info_WZ]):

    '''# @TEMP >>>
    if data is not None:
        data = append_fields(data, 'DSID', np.zeros((data.size,)), dtypes=int)
        for idx in info_data['id']:    
            msk = (data['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            tmp_DSID = info_data['DSID'][idx]  # Get DSID for this file
            data['weight'][msk] *= xsec[tmp_DSID] # Scale by cross section x filter eff. for this DSID
            data['DSID']  [msk] = tmp_DSID        # Store DSID
            pass
        #data['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass
    # @TEMP <<<'''

    if signal is not None:
        signal = append_fields(signal, 'DSID', np.zeros((signal.size,)), dtypes=int)
        for idx in info_sig['id']:    
            msk = (signal['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            tmp_DSID = info_sig['DSID'][idx]  # Get DSID for this file
            signal['weight'][msk] *= xsec[tmp_DSID] # Scale by cross section x filter eff. for this DSID
            signal['DSID']  [msk] = tmp_DSID        # Store DSID
            pass
        signal['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass

    if WZ is not None:
        WZ = append_fields(WZ, 'DSID', np.zeros((WZ.size,)), dtypes=int)
        for idx in info_WZ['id']:    
            msk = (WZ['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            tmp_DSID = info_WZ['DSID'][idx]  # Get DSID for this file
            WZ['weight'][msk] *= xsec[tmp_DSID] # Scale by cross section x filter eff. for this DSID
            WZ['DSID']  [msk] = tmp_DSID        # Store DSID
            pass
        # @TODO: k-factors?
        WZ['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass

    # Check output.
    if data.size == 0 or ((signal is not None) and signal.size == 0):
        warning("No data was loaded. Exiting.")
        return 

    # Compute new variables
    data   = append_fields(data,   'logpt', np.log(data  ['pt']))
    if signal is not None:
        signal = append_fields(signal, 'logpt', np.log(signal['pt']))
        pass
    if WZ is not None:
        WZ = append_fields(WZ,     'logpt', np.log(WZ    ['pt']))
        pass
    

    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Pass/fail masks
    msk_data_pass = tf.config['pass'](data)
    msk_data_fail = ~msk_data_pass
    if signal is not None:
        msk_sig_pass  = tf.config['pass'](signal)
        msk_sig_fail  = ~msk_sig_pass
        pass
    if WZ is not None:
        msk_WZ_pass   = tf.config['pass'](WZ)
        msk_WZ_fail   = ~msk_WZ_pass
        pass

    
    calc = tf.calculator(data=data, config=tf.config, subtract=WZ if args.subtractWZMC else None) # Using default configuration
    calc.mass = args.mass

    # Perform full fit
    if args.window is None:
        calc.fullfit()
       
        print "  -- Computing data weights"
        w_nom,     w_up,     w_down      = calc.fullweights(data  [msk_data_fail])
        if signal is not None:
            print "  -- Computing signal weights"
            w_sig_nom, w_sig_up, w_sig_down  = calc.fullweights(signal[msk_sig_fail])
            pass
        if WZ is not None:
            print "  -- Computing W/Z weights"
            w_WZ_nom, w_WZ_up, w_WZ_down = calc.fullweights(WZ    [msk_WZ_fail])
        else:
            w_WZ_nom, w_WZ_up, w_WZ_down = None, None, None
            pass
        print "  -- Final fit done"
        if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/tf_', MC=False)


    # Perform fit with manually-set window size
    else:    
        # @TODO: - Forcing the fit to use same length scale as 20% window fit. Improve?
        calc.window = 0.2
        calc.fit()
        theta = calc.theta()
        calc.window = args.window
        calc.fit(theta=theta)
       
        print "  -- Computing data weights"
        w_nom  = calc.weights(data[msk_data_fail])
        w_up   = calc.weights(data[msk_data_fail], shift=+1)
        w_down = calc.weights(data[msk_data_fail], shift=-1)
        if signal is not None:
            print "  -- Computing signal weights"
            w_sig_nom  = calc.weights(signal[msk_sig_fail])
            w_sig_up   = calc.weights(signal[msk_sig_fail], shift=+1)
            w_sig_down = calc.weights(signal[msk_sig_fail], shift=-1)
            pass
        if WZ is not None:
            print "  -- Computing W/Z weights"
            w_WZ_nom, = calc.weights(WZ[msk_WZ_fail])
            w_WZ_up   = calc.weights(WZ[msk_WZ_fail], shift=+1)
            w_WZ_down = calc.weights(WZ[msk_WZ_fail], shift=-1)
        else:
            w_WZ_nom, w_WZ_up, w_WZ_down = None, None, None
            pass
        print "  -- Manual fit done"
        if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/tf_', MC=False)
        pass


    
    # Computing data-driven background estimate
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    check_make_dir('output')
    
    # Write TF-scaled failing data to file
    output = ROOT.TFile('output/objdef_TF_{DSID:6d}.root'.format(DSID=DSID), 'RECREATE')
    
    for shift, w, w_WZ in zip([0,        1,       -1], 
                              [w_nom,    w_up,    w_down], 
                              [w_WZ_nom, w_WZ_up, w_WZ_down]):
        
        # -- Get branch name for current variation
        var_name = 'Nominal' if shift == 0 else ('TF_UP' if shift == 1 else 'TF_DOWN')

        # -- Prepare mass- and weight vectors
        vector_m = data['m']     [msk_data_fail]
        vector_w = data['weight'][msk_data_fail] * w
        if args.subtractWZdata:
            if WZ is not None and WZ.size > 0:
                print "  Subtracting TF-scaled W/Z MC from background estimate"
                vector_m = np.concatenate((vector_m,   WZ['m']     [msk_WZ_fail]))
                vector_w = np.concatenate((vector_w, - WZ['weight'][msk_WZ_fail] * w_WZ))
            else:
                warning("  Could not subtract failed, TF-scale W/Z MC component")
                pass
            pass
        # Note: Don't subtract the signal component; that's output as a separate histogram to be used in the simultaneous fit
        
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
        print "  Writing arrays to file: %s" % var_name
        treename1 = tf.config['tree'].replace('NumLargeRadiusJets', 'Jet_tau21DDT').replace('Nominal', var_name)
        make_directories('/'.join(treename1.split('/')[:-1]), fromDir=output)
        tree1 = ROOT.TTree(treename1.split('/')[-1], "")
        array2tree(array1, tree=tree1)
        
        # outputTree
        treename2 = tf.config['outputtree'].replace('Nominal', var_name)
        make_directories('/'.join(treename2.split('/')[:-1]), fromDir=output)
        tree2 = ROOT.TTree(treename2.split('/')[-1], "")
        array2tree(array2, tree=tree2)

        output.Write()        
        pass

    output.Close()
        

    # Write TF-scaled failing signal MC to file
    if signal is not None:
        output = ROOT.TFile('output/objdef_TF_{DSID:6d}_signalfail.root'.format(DSID=DSID), 'RECREATE')
        
        for shift, w_sig in zip([0, 1, -1], [w_sig_nom, w_sig_up, w_sig_down]):
            
            # -- Get branch name for current variation
            var_name = 'Nominal' if shift == 0 else ('TF_UP' if shift == 1 else 'TF_DOWN')
            
            # -- Prepare mass- and weight vectors
            vector_m = signal['m']     [msk_sig_fail]
            vector_w = signal['weight'][msk_sig_fail] * w_sig
            
            # -- Prepare DISD and isMC vectors
            vector_DSID = np.ones_like(vector_w) * (DSID + 1E+05)
            vector_isMC = np.ones_like(vector_w).astype(bool)
            
            array1 = np.array(zip(vector_m, vector_w),
                              dtype = [(tf.config['prefix'] + 'm', np.float64),
                                       ('weight',                  np.float64)])
            
            array2 = np.array(zip(vector_DSID, vector_isMC),
                              dtype = [('DSID', np.uint32),
                                       ('isMC', np.bool_)])
            
            # Mass and weight branch
            print "  Writing arrays to file: %s" % var_name
            treename1 = tf.config['tree'].replace('NumLargeRadiusJets', 'Jet_tau21DDT').replace('Nominal', var_name)
            make_directories('/'.join(treename1.split('/')[:-1]), fromDir=output)
            tree1 = ROOT.TTree(treename1.split('/')[-1], "")
            array2tree(array1, tree=tree1)
            
            # outputTree
            treename2 = tf.config['outputtree'].replace('Nominal', var_name)
            make_directories('/'.join(treename2.split('/')[:-1]), fromDir=output)
            tree2 = ROOT.TTree(treename2.split('/')[-1], "")
            array2tree(array2, tree=tree2)
            
            output.Write()        
            pass

        output.Close()
        pass
        
    # Save configuration
    check_make_dir('logs')
    
    # -- Turn numpy arrays into lists, in order to make them JSON serializable
    cfg = make_serializable(tf.config)

    json.dump([cfg, vars(args)], open('logs/tf_config_%d.log' % DSID, 'w'))        
    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
