#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for studying the jet mass peak at different pT-thresholds.

@file:   masspeak.py
@author: Andreas Søgaard
@date:   27 September 2017
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

parser = argparse.ArgumentParser(description='Study jet mass peak(s) at different pT.')

parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')


# Main function.
def main ():

    # Parse command-line arguments
    args = parser.parse_args()

    new_pt = 450


    # W/Z MC
    # --------------------------------------------------------------------------

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Load data (W/Z MC)
    files = glob.glob(tf.config['base_path'] + 'objdef_MC_3054*.root')

    if len(files) == 0:
        warning("No files found.")
        return

    data = loadData(files, tf.config['finaltree'], prefix=tf.config['prefix']) 
    info = loadData(files, tf.config['outputtree'], stop=1)
    
    # Check output.
    if data.size == 0:
        warning("No data was loaded. Exiting.")
        return 

    # Scale by cross-section, generator filter efficiency, and luminosity
    xsec = loadXsec(tf.config['xsec_file'])
    data = scale_weights(data, info, xsec)


    # Show normalised W/Z jet mass peak for different pT thresholds
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
    msk = data['pt'] > new_pt
  
    # Create canvas
    c = ap.canvas(batch=not args.show)
    bins = np.linspace(50, 110, 30 + 1, endpoint=True)

    # Plot histograms
    h1 = c.hist(data['m'],      bins=bins, weights=data['weight'],      normalise=True, 
                linecolor=ROOT.kRed   - 4, label='p_{T} > 200 GeV')
    h2 = c.hist(data['m'][msk], bins=bins, weights=data['weight'][msk], normalise=True, 
                linecolor=ROOT.kAzure + 7, label='p_{T} > %d GeV' % new_pt)

    # Decorations
    c.legend()
    c.xlabel("Large-radius jet mass [GeV]")
    c.ylabel("Jets / {:.0f} GeV (normalised)".format(np.diff(bins)[0]))

    c.text(["#sqrt{s} = 13 TeV,  L = %s fb^{-1}" % tf.config['lumi'],
            "Photon channel",
            "W/Z + #gamma Monte Carlo",
       ], qualifier='Internal')
        
    # Show/Save
    if args.show: c.show()
    if args.save: c.save('plots/masspeak_WZ.pdf')


    # Signal(220) MC
    # --------------------------------------------------------------------------

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Load data (W/Z MC)
    files = [sorted(glob.glob(tf.config['base_path'] + 'objdef_MC_30836*.root'))[-1]]

    if len(files) == 0:
        warning("No files found.")
        return

    data = loadData(files, tf.config['finaltree'], prefix=tf.config['prefix']) 
    info = loadData(files, tf.config['outputtree'], stop=1)
    
    # Check output.
    if data.size == 0:
        warning("No data was loaded. Exiting.")
        return 

    # Scale by cross-section, generator filter efficiency, and luminosity
    xsec = loadXsec(tf.config['xsec_file'])
    data = scale_weights(data, info, xsec)

    # Get jet channel plot
    path = '/afs/cern.ch/user/l/lkaplan/public/forAndreas/signals_fatjetunc/sig_220.root'
    f = ROOT.TFile(path, 'READ')
    t = f.Get('Signal_ISRjet')
    data_isrjet = tree2array(t)
    print data_isrjet, data_isrjet.shape


    # Show normalised signal jet mass peak for different pT thresholds
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
    msk = data['pt'] > new_pt
  
    # Create canvas
    c = ap.canvas(batch=not args.show)
    bins = np.linspace(50, 300, 25 + 1, endpoint=True)

    # Plot histograms
    h0 = c.hist(data_isrjet['mJ'], bins=bins, weights=data_isrjet['weight'], normalise=True, 
                linecolor=ROOT.kRed - 4, label='p_{T} > 450 GeV (jet ch.)')
    h0 = c.hist(h0, option='E2', fillcolor=ROOT.kRed - 4, alpha=0.3)

    h1 = c.hist(data['m'],      bins=bins, weights=data['weight'],      normalise=True, 
                linecolor=ROOT.kAzure + 4, label='p_{T} > 200 GeV (#gamma ch.)')
    h1 = c.hist(h1, option='E2', fillcolor=ROOT.kAzure + 4, alpha=0.3)

    h2 = c.hist(data['m'][msk], bins=bins, weights=data['weight'][msk], normalise=True, 
                linecolor=ROOT.kAzure + 7, label='p_{T} > %d GeV (#gamma ch.)' % new_pt)
    h2 = c.hist(h2, option='E2', fillcolor=ROOT.kAzure + 7, alpha=0.3)

    # Decorations
    c.legend(xmin=0.5)
    c.xlabel("Large-radius jet mass [GeV]")
    c.ylabel("Jets / {:.0f} GeV (normalised)".format(np.diff(bins)[0]))

    c.text(["#sqrt{s} = 13 TeV",
            "Z'(220 GeV) + #gamma/jet",
       ], qualifier='Internal')
        
    # Show/Save
    if args.show: c.show()
    if args.save: c.save('plots/masspeak_sig220.pdf')

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
