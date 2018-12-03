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
    calc.fit()
    w_nom  = calc.weights(data[msk_fail])
    w_up   = calc.weights(data[msk_fail], shift=+1.0)
    w_down = calc.weights(data[msk_fail], shift=-1.0)
    # @TEMP
    #calc.fullfit()
    #w_nom, w_up, w_down = calc.fullweights(data[msk_fail])
    if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/closure_temp_restricted_')


    # Comparing jet mass distrbutions (closure)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    if args.show or args.save:

        c = ap.canvas(num_pads=2, batch=not args.show)
        p0, p1 = c.pads()

        bins = tf.config['massbins']
        
        h_bkg  = c.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w_nom,
                        fillcolor=ROOT.kAzure + 7,
                        label='Background pred.')
        h_err  = c.hist(h_bkg, 
                        fillstyle=3245, fillcolor=ROOT.kGray+2, linecolor=ROOT.kGray + 3,
                        label='Stat. uncert.', option='E2')
        h_up   = c.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w_up,
                        linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST',
                        label='Syst. uncert.')
        h_down = c.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w_down,
                        linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST')
        h_data = c.plot(data['m'][msk_pass], bins=bins, weights=data['weight'][msk_pass],
                        label='Pseudo-data')


        print "== Mass: {:.0f} GeV".format(args.mass)
        def get_contents (h):
            return np.asarray([h.GetBinContent(ibin) for ibin in range(1, h.GetXaxis().GetNbins() + 1)])
        def get_errors (h):
            return np.asarray([h.GetBinError  (ibin) for ibin in range(1, h.GetXaxis().GetNbins() + 1)])
        bincentres = bins[:-1] + 0.5 * np.diff(bins)
        
        cd = get_contents(h_data)
        ed = get_errors  (h_data)
        cb = get_contents(h_bkg)
        eb = np.maximum(np.abs(get_contents(h_up)   - get_contents(h_bkg)),
                        np.abs(get_contents(h_down) - get_contents(h_bkg)))
        e = np.sqrt(ed**2 + eb**2)

        msk = np.abs(bincentres - args.mass) / args.mass < args.window

        pulls = (cd - cb) / e

        mean = pulls[msk].mean()
        std  = pulls[msk].std()
        eom  = std / np.sqrt(msk.sum())
        eom1 = 1. / np.sqrt(msk.sum())
        sign  = abs(mean) / eom
        sign1 = abs(mean) / eom1
        print "Pull mean:", mean
        print "Pull mean(abs):", np.abs(pulls[msk]).mean()
        print "Pull std.:", std
        print "Pull eom.:", eom
        print "Pull sign:", sign
        print "Pull sign (wrt. std=1.0):", sign1

        #exit()
        c.ratio_plot((h_err,  h_bkg), option='E2')
        c.ratio_plot((h_up,   h_bkg), option='HIST')
        c.ratio_plot((h_down, h_bkg), option='HIST')
        c.ratio_plot((h_data, h_bkg))

        c.xlabel('Signal jet mass [GeV]')
        c.ylabel('Events')
        p1.ylabel('Data / Est.')

        c.ylim(1E+01, 1E+07)
        p1.ylim(0.8, 1.2)
        p1.yline(1.0)
        c.region("SR", 0.8 * args.mass, 1.2 * args.mass)

        
        #for x in [args.mass * (1 - args.window), args.mass * (1 + args.window)]:
        #    p0.line(x, 1E+01, x, 2E+04)
        #    pass
        #p1.xlines([args.mass * (1 - args.window), args.mass * (1 + args.window)])

        c.text(["#sqrt{s} = 13 TeV,  L = %s fb^{-1}" % tf.config['lumi'],
                "Sherpa incl. #gamma MC",
                "Trimmed anti-k_{t}^{R=1.0} jets",
                "ISR #gamma selection",
                "Window: %d GeV #pm %d %%" % (args.mass, args.window * 100.)
               ], qualifier='Simulation Internal')

        c.log()
        c.legend()

        if args.save: c.save('plots/closure_%dGeV_pm%d.pdf' % (args.mass, args.window * 100.))
        if args.show: c.show()
        pass

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
