#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for studying the closure of the transfer factor method on MC.

@file:   closure.py
@author: Andreas Søgaard
@date:   4 April 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys, glob

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

parser = argparse.ArgumentParser(description='Perform closure test.')

parser.add_argument('--mass', dest='mass', type=int,
                    required=True,
                    help='Center of excluded mass window')
parser.add_argument('--window', dest='window', type=float,
                    default=0.2,
                    help='Relative width of excluded mass window (default: 0.2)')
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
    files  = glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root')

    if len(files) == 0:
        warning("No files found.")
        return

    data = loadData(files, tf.config['tree'], prefix=tf.config['prefix']) 
    info = loadData(files, tf.config['outputtree'], stop=1)

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    data = append_fields(data, 'DSID', np.zeros((data.size,)), dtypes=int)
    for idx in info['id']:    
        msk = (data['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
        DSID = info['DSID'][idx]  # Get DSID for this file
        data['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
        data['DSID']  [msk] = DSID        # Store DSID
        pass
    data['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity

    # Check output.
    if data.size == 0:
        warning("No data was loaded.")
        return 

    # Compute new variables
    data = append_fields(data, 'logpt', np.log(data['pt']))


    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Pass/fail masks
    msk_pass = tf.config['pass'](data)
    msk_fail = ~msk_pass

    # Transfer factor calculator instance
    calc = tf.calculator(data=data, config=tf.config) # Using default configuration
    calc.mass   = args.mass
    calc.window = args.window
    # ... calc.partialbins, calc.emptybins, ...
    calc.fit() # ...(theta=0.5)
    w_nom  = calc.weights(data[msk_fail])
    w_up   = calc.weights(data[msk_fail], shift=+1)
    w_down = calc.weights(data[msk_fail], shift=-1)
    if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/new_closure_')


    # Comparing jet mass distrbutions (closure)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    if args.show or args.save:

        c = ap.canvas(num_pads=2, batch=not args.show)
        p0, p1 = c.pads()

        bins = tf.config['massbins']

        h_bkg  = c.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w_nom,  display=False)
        h_up   = c.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w_up,   display=False)
        h_down = c.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w_down, display=False)
        h_data = c.plot(data['m'][msk_pass], bins=bins, weights=data['weight'][msk_pass],          display=False)

        for bin in range(1, h_bkg.GetXaxis().GetNbins() + 1):
            width = float(h_bkg.GetBinWidth(bin))
            h_bkg .SetBinContent(bin, h_bkg .GetBinContent(bin) / width)
            h_bkg .SetBinError  (bin, h_bkg .GetBinError  (bin) / width)
            h_up  .SetBinContent(bin, h_up  .GetBinContent(bin) / width)
            h_up  .SetBinError  (bin, h_up  .GetBinError  (bin) / width)
            h_down.SetBinContent(bin, h_down.GetBinContent(bin) / width)
            h_down.SetBinError  (bin, h_down.GetBinError  (bin) / width)
            h_data.SetBinContent(bin, h_data.GetBinContent(bin) / width)
            h_data.SetBinError  (bin, h_data.GetBinError  (bin) / width)
            pass

        h_bkg  = c.hist(h_bkg, fillcolor=ROOT.kAzure + 7, label='Background est.')
        h_err  = c.hist(h_bkg, fillstyle=3245, fillcolor=ROOT.kGray+2, linecolor=ROOT.kGray + 3, label='Stat. uncert.', option='E2')
        h_up   = c.hist(h_up, linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST', label='Syst. uncert.')
        h_down = c.hist(h_down, linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST')
        h_data = c.plot(h_data, label='Pseudo-data')
        
        c.ratio_plot((h_err,  h_bkg), option='E2')
        c.ratio_plot((h_up,   h_bkg), option='HIST')
        c.ratio_plot((h_down, h_bkg), option='HIST')
        c.ratio_plot((h_data, h_bkg))

        c.xlabel('Large-#it{R} jet mass [GeV]')
        c.ylabel('Events / GeV')
        p1.ylabel('Data / Est.')

        c.ylim(1E+00, 1E+06)
        p1.ylim(0.80, 1.20)
        p1.yline(1.0)
        c.region("SR", 0.8 * args.mass, 1.2 * args.mass)

        #for x in [args.mass * (1 - args.window), args.mass * (1 + args.window)]:
        #    p0.line(x, 1E+01, x, 2E+04)
        #    pass
        #p1.xlines([args.mass * (1 - args.window), args.mass * (1 + args.window)])

        c.text(["#sqrt{s} = 13 TeV,  %s fb^{-1}" % tf.config['lumi'],
                "Incl. #gamma Monte Carlo",
                "Photon channel",
               ], qualifier='Simulation Internal')

        c.log()
        c.legend()

        if args.save: c.save('plots/new_closure_%dGeV_pm%d.pdf' % (args.mass, args.window * 100.))
        if args.show: c.show()
        pass

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
