#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for producing cutflow numbers.

@file:   tau21distributions.py
@author: Andreas Søgaard
@date:   30 September 2017
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

parser = argparse.ArgumentParser(description='Compute cutflow numbers')

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
    paths = glob.glob('/afs/cern.ch/user/a/asogaard/Analysis/2016/BoostedJetISR/AnalysisTools/outputObjdef/objdef_MC_30836*.root')

    base = 'BoostedJet+ISRgamma/Nominal/'
    trees = [
        'PreSelection/Nominal/GRL/Postcut',
        'PreSelection/Nominal/HLT_g140_loose/Postcut',
        'LargeRadiusJets/Nominal/eta/Postcut',
        'LargeRadiusJets/Nominal/pt/Postcut',
        'LargeRadiusJets/Nominal/dPhiPhoton/Postcut',
        'LargeRadiusJets/Nominal/BoostedRegime/Postcut',
        'LargeRadiusJets/Nominal/rhoDDT/Postcut',
        'EventSelection/Pass/NumPhotons/Postcut',
        'EventSelection/Pass/NumLargeRadiusJets/Postcut',
        'EventSelection/Pass/Jet_tau21DDT/Postcut'
        ]
    trees = [base + tree for tree in trees]

    names = [
        "GRL",
        "Trigger",
        "Jet eta",
        "Jet pT",
        "Jet dPhi(y)",
        "Jet pT > 2m",
        "Jet rhoDDT",
        "Photon pT",
        "(Jet count)",
        "Jet tau21DDT"
        ]

    DSIDs = [308363, 308364, 308365, 308366, 308367]

    basecounts = base = {
        308363: 50000,
        308364: 48000,
        308365: 50000,
        308366: 50000,
        308367: 48000,
        }

    genfilteff = {
        308363: 2.9141e-01,
        308364: 1.3795e-01,
        308365: 7.5920e-02,
        308366: 4.5985e-02,
        308367: 2.9914e-02,
        }

    for idx, DSID in enumerate(DSIDs):
        print "\nDSID: {}".format(DSID)

        passingEvents = None
        for name, tree in zip(names,trees):
            
            # Getting data
            data = loadData([paths[idx]], tree)
            
            if passingEvents is None:
                passingEvents = sorted(list(set(data['eventNumber'])))
                assert len(passingEvents) == data.shape[0]
            else:
                passingEvents = sorted(list(set(passingEvents) & set(data['eventNumber'])))
                pass
            print "    {:15s}: {:6d} | {:4.1f}% | {:4.1f}%".format(
                name, 
                len(passingEvents), 
                len(passingEvents) / float(basecounts[DSID]) * 100., 
                len(passingEvents) / float(basecounts[DSID]) * genfilteff[DSID] * 100., 
                )
            pass
        pass
    return 

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
            #( 50, 100),
            (100, 150),
            (150, 200),
            (200, 250),
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

        if xvar in ['tau21DDT/Postcut', 'tau21'] and args.save:
            f = ROOT.TFile('output/hists_isrgamma_%sdistributions.root' % xvar, 'RECREATE')
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

