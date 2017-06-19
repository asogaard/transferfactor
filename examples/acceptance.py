#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for ...

@file:   acceptance.py
@author: Andreas Søgaard
@date:   16 June 2017
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
    files  = glob.glob(tf.config['base_path'] + 'objdef_MC_30836*.root')

    if len(files) == 0:
        warning("No files found.")
        return

    cats = ['Pass', 'Fail']
    data = {cat: loadData(files, tf.config['finaltree'].replace('Pass', cat), prefix=tf.config['prefix']) for cat in cats}
    info = loadData(files, tf.config['outputtree'], stop=1)

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    for cat in cats:
        data[cat] = append_fields(data[cat], 'DSID', np.zeros((data[cat].size,)), dtypes=int)
        for idx in info['id']:    
            msk = (data[cat]['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            DSID = info['DSID'][idx]  # Get DSID for this file
            data[cat]['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
            data[cat]['DSID']  [msk] = DSID        # Store DSID
            pass
        data[cat]['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass

    # Check output.
    if 0 in [data[cat].size for cat in cats]:
        warning("No data was loaded.")
        return 

    # Compute new variables
    for cat in cats:
        data[cat] = append_fields(data[cat], 'logpt', np.log(data[cat]['pt']))
        pass


    # Plot acceptance versus signal mass point
    # --------------------------------------------------------------------------
    
    DSIDs = sorted(list(set(data['Pass']['DSID'])))
    base_events     = np.array([50000., 48000., 50000., 50000., 48000.], dtype=float)
    accepted_events = np.zeros_like(base_events)
    passed_events   = np.zeros_like(base_events)

    for idx, DSID in enumerate(DSIDs):
        num_pass = np.sum(data['Pass']['DSID'] == DSID)
        num_fail = np.sum(data['Fail']['DSID'] == DSID)

        accepted_events[idx] = num_pass + num_fail
        passed_events  [idx] = num_pass
        pass

    acceptance  = accepted_events / base_events
    efficiency  = passed_events / accepted_events


    # Binomial uncertainty of fraction 'f' given 'base_events' trials.
    acceptance_uncertainties = np.sqrt(acceptance * (1. - acceptance) / base_events)
    efficiency_uncertainties = np.sqrt(efficiency * (1. - efficiency) / accepted_events)
    acc_eff_uncertainties    = np.sqrt(acceptance * efficiency * (1. - acceptance * efficiency)  / accepted_events)

    masses = np.array([100, 130, 160, 190, 220], dtype=float)

    acc_e = ROOT.TGraphErrors(len(DSIDs), masses, acceptance, np.zeros_like(acceptance_uncertainties), acceptance_uncertainties)
    acc   = ROOT.TGraphErrors(len(DSIDs), masses, acceptance)

    eff_e = ROOT.TGraphErrors(len(DSIDs), masses, efficiency, np.zeros_like(efficiency_uncertainties), efficiency_uncertainties)
    eff   = ROOT.TGraphErrors(len(DSIDs), masses, efficiency)

    acceff_e = ROOT.TGraphErrors(len(DSIDs), masses, acceptance * efficiency, np.zeros_like(acc_eff_uncertainties), acc_eff_uncertainties)
    acceff   = ROOT.TGraphErrors(len(DSIDs), masses, acceptance * efficiency)

    # Plot: Acceptance
    c = ap.canvas(batch=not args.show)
    c.graph(acc_e, linecolor=ROOT.kAzure+2, linewidth=2, fillcolor=ROOT.kAzure+7, option='A3', label='ISR #gamma channel', legend_option='FL') 
    c.graph(acc,   linecolor=ROOT.kAzure+2, linewidth=2, option='L')
    c.xlabel("Signal mass point [GeV]")
    c.ylabel("Acceptance (Z' #rightarrow large-#it{R} jet + #gamma)")
    c.text(["#sqrt{s} = 13 TeV"], qualifier="Simulation Internal")
    c.legend()
    if args.save: c.save('plots/acceptance.pdf')
    if args.show: c.show()

    # Plot: Efficiency
    c = ap.canvas(batch=not args.show)
    c.graph(eff_e, linecolor=ROOT.kAzure+2, linewidth=2, fillcolor=ROOT.kAzure+7, option='A3', label='ISR #gamma channel', legend_option='FL') 
    c.graph(eff,   linecolor=ROOT.kAzure+2, linewidth=2, option='L')
    c.xlabel("Signal mass point [GeV]")
    c.ylabel("Efficiency (Z' #rightarrow large-#it{R} jet + #gamma)")
    c.text(["#sqrt{s} = 13 TeV"], qualifier="Simulation Internal")
    c.legend()
    if args.save: c.save('plots/efficiency.pdf')
    if args.show: c.show()

    # Plot: Acceptance x efficiency
    c = ap.canvas(batch=not args.show)
    c.graph(acceff_e, linecolor=ROOT.kAzure+2, linewidth=2, fillcolor=ROOT.kAzure+7, option='A3', label='ISR #gamma channel', legend_option='FL') 
    c.graph(acceff,   linecolor=ROOT.kAzure+2, linewidth=2, option='L')
    c.xlabel("Signal mass point [GeV]")
    c.ylabel("Acc #times eff (Z' #rightarrow large-#it{R} jet + #gamma)")
    c.text(["#sqrt{s} = 13 TeV"], qualifier="Simulation Internal")
    c.legend()
    if args.save: c.save('plots/acceptance_times_efficiency.pdf')
    if args.show: c.show()

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
