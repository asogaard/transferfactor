#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for producing global background shape (GBS) inputs for the dijet + ISR analysis.

@file:   gbs.py
@author: Andreas Søgaard 
@date:   7 June 2017
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

parser = argparse.ArgumentParser(description='Compute global background shape (GBS).')

parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')
parser.add_argument('--data', dest='data', action='store_const',
                    const=True, default=False,
                    help='Use data (default: False/MC)')
parser.add_argument('--subtractWZdata', dest='subtractWZdata', action='store_const',
                    const=True, default=False,
                    help='Subtract estimated W/Z fail component from data (default: False)')
parser.add_argument('--subtractWZMC', dest='subtractWZMC', action='store_const',
                    const=True, default=False,
                    help='Subtract W/Z MC from TF profile (default: False)')


# Main function
def main ():
    
    # Parse command-line arguments
    args = parser.parse_args()

    # Check(s)
    if (not args.data) and args.subtractWZMC:
        warning("Requesting to subtract W/Z MC from MC background which contains no contamination. Exiting.")
        return

    if (not args.data) and args.subtractWZdata:
        warning("Requesting to subtract W/Z data from MC background which contains no contamination. Exiting.")
        return

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Load data
    if args.data:
        files_data = glob.glob(tf.config['base_path'] + 'objdef_data_*.root')
    else:
        files_data = glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root')
        pass
    files_WZ   = glob.glob(tf.config['base_path'] + 'objdef_MC_30543*.root') + \
                 glob.glob(tf.config['base_path'] + 'objdef_MC_30544*.root')
    
    if len(files_data) == 0:
        warning("No files found. Try to run:")
        warning(" $ source getSomeData.sh")
        return

    data      = loadData(files_data, tf.config['tree'], prefix=tf.config['prefix'])
    WZ        = loadData(files_WZ,   tf.config['tree'], prefix=tf.config['prefix'])
    info_data = loadData(files_data, tf.config['outputtree'], stop=1)
    info_WZ   = loadData(files_WZ,   tf.config['outputtree'], stop=1)

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])

    # Append new DSID field
    data = append_fields(data, 'DSID', np.zeros((data.size,)), dtypes=int)
    for idx, id in enumerate(info_data['id']):
        msk = (data['id'] == id) # Get mask of all 'data' entries with same id, i.e. from same file
        tmp_DSID = info_data['DSID'][idx]  # Get DSID for this file
        if not args.data:
            data['weight'][msk] *= xsec[tmp_DSID] # Scale by cross section x filter eff. for this DSID
            data['DSID']  [msk] = tmp_DSID        # Store DSID
            pass
        pass
    if not args.data:
        data['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass
    
    WZ = append_fields(WZ, 'DSID', np.zeros((WZ.size,)), dtypes=int)
    for idx in info_WZ['id']:
        msk = (WZ['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
        tmp_DSID = info_WZ['DSID'][idx]  # Get DSID for this file
        WZ['weight'][msk] *= xsec[tmp_DSID] # Scale by cross section x filter eff. for this DSID
        WZ['DSID']  [msk] = tmp_DSID        # Store DSID
        pass
    WZ['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
    
    # Compute new variables
    data = append_fields(data, 'logpt', np.log(data['pt']))
    WZ   = append_fields(WZ  , 'logpt', np.log(WZ  ['pt']))

 
    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    # Pass/fail masks
    msk_pass    = tf.config['pass'](data)
    msk_fail    = ~msk_pass
    msk_WZ_pass = tf.config['pass'](WZ)
    msk_WZ_fail = ~msk_WZ_pass

    # Transfer factor calculator instance
    calc = tf.calculator(data=data, config=tf.config, subtract=WZ if (args.subtractWZMC and args.data) else None)

    # GBS mass bins
    masses = np.linspace(100, 270, 34 + 1, endpoint=True) # GBS mass bins

    # Weight and counter arrays
    weights_bkg_nom  = np.zeros((np.sum(msk_fail),), dtype=float)
    weights_bkg_up   = np.zeros((np.sum(msk_fail),), dtype=float)
    weights_bkg_down = np.zeros((np.sum(msk_fail),), dtype=float)
    counter_bkg      = np.zeros((np.sum(msk_fail),), dtype=float)

    weights_WZ_nom  = np.zeros((np.sum(msk_WZ_fail),), dtype=float)
    weights_WZ_up   = np.zeros((np.sum(msk_WZ_fail),), dtype=float)
    weights_WZ_down = np.zeros((np.sum(msk_WZ_fail),), dtype=float)
    counter_WZ      = np.zeros((np.sum(msk_WZ_fail),), dtype=float)

    #ctemp = ap.canvas(batch=True)
    for mass in masses:
        print " --", mass

        # Fit TF profile
        calc.mass = mass
        calc.fullfit()

        if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/tf_gbs_%s_%dGeV_' % ('data' if args.data else 'MC', mass), MC=not args.data)

        # Get TF weights
        w_nom,    w_up,    w_down    = calc.fullweights(data[msk_fail])
        w_WZ_nom, w_WZ_up, w_WZ_down = calc.fullweights(WZ  [msk_WZ_fail])

        # Compute mask for which jets to use in GBS computation
        msk_gbs    = ~(np.abs(data[msk_fail]   ['m'] - mass) < 0.2 * mass)
        msk_gbs_WZ = ~(np.abs(WZ  [msk_WZ_fail]['m'] - mass) < 0.2 * mass)

        # Store weights and increment counter for masked jets
        weights_bkg_nom [msk_gbs] += w_nom [msk_gbs]
        weights_bkg_up  [msk_gbs] += w_up  [msk_gbs]
        weights_bkg_down[msk_gbs] += w_down[msk_gbs]
        counter_bkg     [msk_gbs] += 1.

        weights_WZ_nom [msk_gbs_WZ] += w_WZ_nom [msk_gbs_WZ]
        weights_WZ_up  [msk_gbs_WZ] += w_WZ_up  [msk_gbs_WZ]
        weights_WZ_down[msk_gbs_WZ] += w_WZ_down[msk_gbs_WZ]
        counter_WZ     [msk_gbs_WZ] += 1.
        pass

    # Take average of jets in signal regions
    msk = (counter_bkg > 0)
    weights_bkg_nom [msk] /= counter_bkg[msk]
    weights_bkg_up  [msk] /= counter_bkg[msk]
    weights_bkg_down[msk] /= counter_bkg[msk]

    msk = (counter_WZ > 0)
    weights_WZ_nom [msk] /= counter_WZ[msk]
    weights_WZ_up  [msk] /= counter_WZ[msk]
    weights_WZ_down[msk] /= counter_WZ[msk]


    # Computing data-driven background estimate
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    check_make_dir('output')

    DSID = 400000 if args.data else 400001

    # Write TF-scaled failing data to file
    output = ROOT.TFile('output/objdef_GBS{MC}_{DSID}.root'.format(DSID=DSID, MC='' if args.data else 'MC'), 'RECREATE')

    for shift, w, w_WZ in zip([0,               1,             -1],
                              [weights_bkg_nom, weights_bkg_up, weights_bkg_down],
                              [weights_WZ_nom,  weights_WZ_up,  weights_WZ_down]):

        # -- Get branch name for current variation
        var_name = 'Nominal' if shift == 0 else ('TF_UP' if shift == 1 else 'TF_DOWN')

        # -- Prepare mass- and weight vectors
        vector_m = data['m']     [msk_fail]
        vector_w = data['weight'][msk_fail] * w
        if args.subtractWZdata and args.data:
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

    # Save configuration
    check_make_dir('logs')

    # -- Turn numpy arrays into lists, in order to make them JSON serializable
    cfg = make_serializable(tf.config)

    json.dump([cfg, vars(args)], open('logs/gbs_config_%s_%d.log' % ('data' if args.data else 'MC', DSID), 'w'))

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
