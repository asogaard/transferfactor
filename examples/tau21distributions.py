#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for producing publication-level plots of tau21 and tau21DDT distributions.

@file:   tau21distributions.py
@author: Andreas Søgaard
@date:   28 June 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys, glob, itertools

# Get ROOT to stop hogging the command-line options
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Scientific import(s)
try:
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
    from rootplotting.style import *
    from snippets.functions import dict_product, displayName, displayUnit, displayNameUnit
except ImportError:
    print "WARNING: This script uses the 'snippets' package. Clone it as e.g.:"
    print "         $ git clone git@github.com:asogaard/snippets.git"
    sys.exit()
    pass

# Command-line arguments parser
import argparse

parser = argparse.ArgumentParser(description='Make tau21 distributions')

parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')


# Main function definition.
def main ():

    # Parse command-line arguments
    args = parser.parse_args()

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Input file paths
    paths = glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root')

    # Load data
    xsec = loadXsec(tf.config['xsec_file'])

    # Standard definitions
    axes = {
        'tau21':    (20, -0.1, 0.9),
        'tau21DDT': (20,  0.0, 1.0), #  0.1, 1.1),
        }
        
    # Getting data
    data = loadData(paths, tf.config['tree'], prefix=tf.config['prefix'])
    info = loadData(paths, tf.config['outputtree'], stop=1)
    
    data = scale_weights(data, info, xsec, lumi=tf.config['lumi'])


    # Plotting tau21(DDT) profiles
    # ----------------------------------------------------------------
       
    slices = {
        'pt': [
            #( 200,  2000), # For inclusive distribution (test)
            ( 200,  300),
            ( 300,  500),
            ( 500, 1000),
            ],
        'm': [
            #(0,1000), # For inclusive distribution (test)
            ( 50, 100),
            (100, 150),
            (150, 200),
            #(200, 250),
            ],
        #'rhoDDT': [
        #    (1.5, 2.5),
        #    (2.5, 4.0),
        #    (4.0, 5.5),
        #    ],
        }

    colours = [ROOT.kRed + i for i in np.arange(0,5,2)] # Overwriting

    keys = slices.keys()
    key = keys[0]
    category_names = ["[%.0f, %.0f] %s" % (slices[key][i][0], slices[key][i][1], displayUnit(key)) for i in range(len(slices[key]))]

    for xvar, axis in axes.iteritems():
        print "Plotting %s:" % xvar

        c = ap.canvas(batch=not args.show)
        bins = np.linspace(axis[1], axis[2], axis[0] + 1, endpoint=True)

        if xvar == 'tau21DDT' and args.save:
            f = ROOT.TFile('output/hists_isrgamma_tau21distributions.root', 'RECREATE')
        else:
            f = None
            pass

        # Fill sliced histograms
        histograms = list()
        for i, sl in enumerate(dict_product(slices)):
            print "  %d:" % i, sl

            # Create full mask for current slice
            msk = np.ones_like(data['weight']).astype(bool)
            for key, limits in sl.iteritems():
                msk &= (data[key] >= limits[0]) & (data[key] < limits[1])
                pass

            # Create distribution for current slice
            key = keys[1]
            label = "[%.0f, %.0f] %s" % (sl[key][0], sl[key][1], displayUnit(key))
            hist = c.hist(data[xvar][msk], bins=bins, weights=data['weight'][msk], linecolor=colours[i % 3], linestyle=1 + (i//3), label=label if i < 3 else None, normalise=True)
            if f:
                f.cd()
                hist.SetName('h_%s_%d_%d_%s_%d_%s' % (keys[0], sl[keys[0]][0], sl[keys[0]][1], keys[1], sl[keys[1]][0], sl[keys[1]][1]))
                #hist.Write()
                pass
            pass

        if f:
            f.Write()
            f.Close()
            pass

        # Decorations
        c.xlabel('Signal jet %s' % displayNameUnit(xvar))
        c.ylabel('Jets (a.u.)')
        c.text(["#sqrt{s} = 13 TeV", "ISR #gamma selection"], qualifier="Simulation Internal")
        c.padding(0.50)

        ymax = 0.735
        c.legend(header='Jet %s in:' % displayName(keys[1]),
                 ymax=ymax)
        c.legend(header='Jet %s in:' % displayName(keys[0]), categories=[
                (category_names[idx], {'linestyle': 1 + idx}) for idx in range(len(category_names))
                ],
                 ymax=ymax, xmin=0.19)

        if args.show: c.show()
        if args.save: c.save('plots/tau21distributions_%s__%s_x_%s.pdf' % (xvar, keys[0], keys[1]))
        pass
        
    return


# Main function call.
if __name__ == '__main__':
   main()
   pass

