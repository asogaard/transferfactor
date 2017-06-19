#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for ...

@file:   acceptance.py
@author: Andreas Søgaard
@date:   16 June 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys, glob, re

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
    files = glob.glob(tf.config['base_path'] + 'objdef_MC_30836*.root')

    if len(files) == 0:
        warning("No files found.")
        return

    colours = [ROOT.kViolet + 7, ROOT.kAzure + 7, ROOT.kTeal, ROOT.kSpring - 2, ROOT.kOrange - 3, ROOT.kPink]
    masses = np.array([100, 130, 160, 190, 220], dtype=float)
    DSIDs  = [int(re.search('objdef_MC_(\d{6})\.root', f).group(1)) for f in files]


    # Plot jet-level acceptance versus signal mass point
    # --------------------------------------------------------------------------

    counts = dict()
    labels = list()
    base_events = dict()
    for DSID, filename in zip(DSIDs, files):
        f = ROOT.TFile(filename, 'READ')
        cutflow = f.Get('BoostedJet+ISRgamma/Nominal/LargeRadiusJets/Nominal/Cutflow')
        ax = cutflow.GetXaxis()
        for bin in range(1, ax.GetNbins() + 1):
            label = ax.GetBinLabel(bin)
            if label not in counts: counts[label] = dict()
            if len(labels) < bin:   labels.append(label)
            counts[label][DSID] = cutflow.GetBinContent(bin)
            pass

        # Get number of events passing pre-selection
        cutflow = f.Get('BoostedJet+ISRgamma/Nominal/PreSelection/Nominal/Cutflow')
        base_events[DSID] = cutflow.GetBinContent(cutflow.GetXaxis().GetNbins())
        pass
    
    acceptances = {DSID: np.zeros((len(labels),)) for DSID in DSIDs}
    for idx, label in enumerate(labels):
        for DSID in DSIDs:
            acceptances[DSID][idx] = counts[label][DSID] / float(base_events[DSID]) #float(counts[labels[0]][DSID])
            pass
        pass

    c = ap.canvas(batch=not args.show)
    names = ["Pre-seleciton", "|#eta| < 2.0", "p_{T} > 200 GeV", "#Delta(#phi,J) > #pi/2", "p_{T} > 2 #times m", "#rho^{DDT} > 1.5"]
    for idx, (label, name) in enumerate(zip(labels, names)):
        acceptance_values = np.array([acceptances[DSID][idx] for DSID in DSIDs])
        acc = ROOT.TGraphErrors(len(DSIDs), masses, acceptance_values)
        c.graph(acc, linecolor=colours[idx], markercolor=colours[idx], linewidth=2, option=('A' if idx == 0 else '') + 'PL', label=name, legend_option='PL') 
        pass

    c.padding(0.5)
    c.xlabel("Signal mass point [GeV]")
    c.ylabel("#LTNum. large-#it{R} jets / event#GT")
    c.text(["#sqrt{s} = 13 TeV", "ISR #gamma channel", "Jet selection"], qualifier="Simulation Internal")
    c.legend()
    if args.save: c.save('plots/acceptance_jet.pdf')
    if args.show: c.show()


    # Plot event-level acceptance versus signal mass point
    # --------------------------------------------------------------------------
    
    treenames = ['BoostedJet+ISRgamma/Nominal/PreSelection/Nominal/HLT_g140_loose/Postcut',
                 tf.config['finaltree'].replace('Jet_tau21DDT', 'NumPhotons'), 
                 tf.config['finaltree'].replace('Jet_tau21DDT', 'NumLargeRadiusJets'), 
                 tf.config['finaltree']]
    data = {treename: loadData(files, treename, prefix=tf.config['prefix']) for treename in treenames}
    info = loadData(files, tf.config['outputtree'], stop=1)

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    for key in data:
        data[key] = append_fields(data[key], 'DSID', np.zeros((data[key].size,)), dtypes=int)
        for idx in info['id']:    
            msk = (data[key]['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            DSID = info['DSID'][idx]  # Get DSID for this file
            data[key]['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
            data[key]['DSID']  [msk] = DSID        # Store DSID
            pass
        data[key]['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass

    # Check output.
    if 0 in [data[key].size for key in data]:
        warning("No data was loaded.")
        return 

    base_events     = np.array([50000., 48000., 50000., 50000., 48000.], dtype=float)
    acceptances = list()

    for key in treenames:
        accepted_events = np.zeros_like(base_events)

        for idx, DSID in enumerate(DSIDs):
            accepted_events[idx] = np.sum(data[key]['DSID'] == DSID)
            pass

        acceptances.append(accepted_events / base_events)
        pass

    # Plot: Event-levle ecceptance
    c = ap.canvas(batch=not args.show)
    names = ["Pre-selection", "= 1 #gamma", "#geq 1 jet", "#tau_{21}^{DDT} < 0.5"]
    for icut, (name, acceptance) in enumerate(zip(names,acceptances)):
        acc = ROOT.TGraphErrors(len(DSIDs), masses, acceptances[icut])
        c.graph(acc, linecolor=colours[icut], markercolor=colours[icut], linewidth=2, option=('A' if icut == 0 else '') + 'PL', label=name, legend_option='PL') 
        pass
    c.xlabel("Signal mass point [GeV]")
    c.ylabel("Acceptance (Z' #rightarrow large-#it{R} jet + #gamma)")
    c.text(["#sqrt{s} = 13 TeV", "ISR #gamma channel", "Event selection"], qualifier="Simulation Internal")
    c.legend()
    if args.save: c.save('plots/acceptance_event.pdf')
    if args.show: c.show()

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
