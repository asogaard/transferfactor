#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for validating the transfer factor fitting and systematics assignement procedure.

@file:   validation.py
@author: Andreas Søgaard
@date:   18 July 2017
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

#parser.add_argument('--mass', dest='mass', type=int,
#                    required=True,
#                    help='Center of excluded mass window')
parser.add_argument('--window', dest='window', type=float,
                    default=None,
                    help='Width of excluded mass window (default: None)')
parser.add_argument('--N', dest='N', type=int,
                    default=10,
                    help='Number of toy experiments (default: 10)')
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


    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Load data
    files = glob.glob(tf.config['base_path'] + 'objdef_data_*.root')

    if len(files) == 0:
        warning("No files found.")
        return

    data = loadData(files, tf.config['tree'], prefix=tf.config['prefix']) 
    info = loadData(files, tf.config['outputtree'], stop=1)
    
    # Check output.
    if data.size == 0:
        warning("No data was loaded. Exiting.")
        return 

    # Compute new variables
    data = append_fields(data, 'logpt', np.log(data['pt']))
    

    # Pass/fail masks
    msk_pass = tf.config['pass'](data)
    msk_fail = ~msk_pass
    
    
    # Validating transfer factor fit using toys
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    #for mass in [85] + list(np.linspace(100, 250, 15 + 1, endpoint=True)):
    for mass in list(np.linspace(110, 250, 14 + 1, endpoint=True)):

        print "-------- MASS: %d GeV" % mass
        
        # Set up transfer factor calculator instance
        calc = tf.calculator(data=data, config=tf.config, verbose=False) # Using default configuration
        calc.mass   = mass
        calc.window = 0.2 if (args.window is None) else args.window
        
        # Get nomnial best-fit theta
        calc.fit()
        theta = calc.theta()
        nominal_weights = calc.weights(data[msk_fail], shift=0), \
            calc.weights(data[msk_fail], shift=+1), \
            calc.weights(data[msk_fail], shift=-1)
        
        # "Throw toys" from TF profile, fit N times
        calc.toysfit(N=args.N, theta=theta)
        
        # Get weights for each toys experiment fit
        toys_weights = calc.toysweights(data[msk_fail])
        
        # Plot variations
        bins = tf.config['massbins']
        c = ap.canvas(num_pads=2, batch=not args.show)
        
        # -- Nominal background(s)
        hist_nom  = c.hist(data[msk_fail]['m'], bins=bins, weights=nominal_weights[0], fillcolor=ROOT.kAzure + 7, label='Nominal bkg.')
        h_sum = c.hist(hist_nom,
                       fillstyle=3245, fillcolor=ROOT.kGray + 2, linecolor=ROOT.kGray + 3, option='E2',
                       label='Stat. uncert.')
        
        
        # -- Toys backgrounds
        toys_hists = list()
        for idx, weights in enumerate(toys_weights):
            h = c.hist(data[msk_fail]['m'], bins=bins, weights=weights[0], fillstyle=0, linecolor=ROOT.kRed + idx % 5, linestyle = 1 + idx // 5, label='Toys %d' % (idx + 1) if idx < 5 else None)
            toys_hists.append(h)
            pass
        
        # -- Nominal variations
        hist_up   = c.hist(data[msk_fail]['m'], bins=bins, weights=nominal_weights[1], fillstyle=0, linecolor=ROOT.kGreen, linestyle=2, label='Syst. uncert.')
        hist_down = c.hist(data[msk_fail]['m'], bins=bins, weights=nominal_weights[2], fillstyle=0, linecolor=ROOT.kGreen, linestyle=2)
        
        
        # -- Data
        hist_data = c.plot(data[msk_pass]['m'], bins=bins, label='Data')
        
        # -- Ratio plots
        c.ratio_plot((h_sum,     hist_nom), option='E2')
        for idx, h in enumerate(toys_hists):
            c.ratio_plot((h, hist_nom), option='HIST')
            pass
        c.ratio_plot((hist_up,   hist_nom), option='HIST')
        c.ratio_plot((hist_down, hist_nom), option='HIST')
        c.ratio_plot((hist_data, hist_nom), oob=True)
        
        # -- Decorations
        c.xlabel('Large-#it{R} jet mass [GeV]')
        c.ylabel('Events / 5 GeV')
        c.pads()[1].ylabel('Ratio wrt. nominal')
        c.pads()[1].ylim(0.8, 1.2)
        c.pads()[1].yline(1.)
        c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
                "Photon channel"],
                qualifier="Internal")

        c.region("SR", 0.8 * mass, 1.2*mass)
        
        c.legend()
        c.log()
        if args.show: c.show()
        if args.save: c.save('plots/validation_%dGeV_N%d.pdf' % (mass, args.N))

        pass
        
    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
