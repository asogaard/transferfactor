#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for optimising the tau21DDT cut and producing publication-level plots.

@file:   cutoptimisation.py
@author: Andreas Søgaard
@date:   26 June 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys, glob, itertools, re

# Get ROOT to stop hogging the command-line options
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

# Scientific import(s)
try:
    import ROOT
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
    print "WARNING: This script uses the 'snippets' package. Clone it as e.g.:"
    print "         $ git clone git@github.com:asogaard/snippets.git"
    sys.exit()
    pass

# Command-line arguments parser
import argparse

parser = argparse.ArgumentParser(description='Perform substructure cut optimisation.')

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
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Input file paths
    paths = {
        'bkg': glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root'),
        'sig': glob.glob(tf.config['base_path'] + 'objdef_MC_30836*.root'),
        }
    
    sig_names = [
        "Z' + #gamma (100 GeV)", 
        "Z' + #gamma (130 GeV)", 
        "Z' + #gamma (160 GeV)", 
        "Z' + #gamma (190 GeV)", 
        "Z' + #gamma (220 GeV)", 
        ]


    # Load data
    treename = 'BoostedJet+ISRgamma/Nominal/EventSelection/Pass/NumLargeRadiusJets/Postcut'

    # Standard definitions
    # -- Data
    prefix = 'Jet_'

    xvars = ['m']
    axis = {
        'm':        (40, 50, 250),
        'tau21DDT': (20, 0, 1),
        }

    # Getting data
    bkg      = loadData(paths['bkg'], treename, prefix=prefix)
    sig      = loadData(paths['sig'], treename, prefix=prefix)
    info_bkg = loadData(paths['bkg'], tf.config['outputtree'], stop=1)
    info_sig = loadData(paths['sig'], tf.config['outputtree'], stop=1)

    # Check output
    if bkg.size == 0 or sig.size == 0:
        warning("No data was loaded.")
        return

    # Scaling by cross-section
    xsec = loadXsec(tf.config['xsec_file'])

    bkg = scale_weights(bkg, info_bkg, xsec, lumi=tf.config['lumi'])
    sig = scale_weights(sig, info_sig, xsec, lumi=tf.config['lumi'], verbose=True)


    # Compute significance improvements
    # --------------------------------------------------------------------------

    # Loop tau21DDT bins
    cuts = np.linspace(axis['tau21DDT'][1], axis['tau21DDT'][2], 10 * axis['tau21DDT'][0] + 1, endpoint=True)
    bins = np.linspace(axis['m']       [1], axis['m']       [2],      axis['m']       [0] + 1, endpoint=True)
    significances = [list() for _ in range(len(set(sig['DSID'])))]
    plot_cut=36

    c_temp = ap.canvas(batch=True)
    for icut, cut in enumerate(cuts):
        print "Bin %2d: tau21DDT < %.3f" % (icut, cut)

        # Create histograms
        msk = bkg['tau21DDT'] < cut
        h_bkg = c_temp.hist(bkg['m'][msk], bins=bins, weights=bkg['weight'][msk], fillcolor=ROOT.kAzure+7, name='Background', display=(icut == plot_cut))
        h_sigs = list()
        for idx, (DSID, name) in enumerate(zip(sorted(list(set(sig['DSID']))), sig_names)):
            msk = (sig['DSID'] == DSID) & (sig['tau21DDT'] < cut)
            h_sigs.append(c_temp.hist(sig['m'][msk], bins=bins, weights=sig['weight'][msk], linestyle=1+idx, label=name, display=(icut == plot_cut)))
            pass

        # Fill histograms
        for i,h_sig in enumerate(h_sigs):
            print " -- Signal %d |" % (i + 1), 
            sign = 0.
            for bin in range(1, h_bkg.GetXaxis().GetNbins() + 1):
                s = h_sig.GetBinContent(bin)
                b = max(h_bkg.GetBinContent(bin), 0.5)
                if icut == 35:
                    print s,b
                    pass
                #sign += np.square(s/np.sqrt(b))
                sign += 2 * ((s + b) * np.log(1 + s/b) - s) # Slide 29 in [http://www2.warwick.ac.uk/fac/sci/physics/research/epp/events/seminars/cowan_warwick_2011.pdf] -- here using the _square_ of the per-bin significance, in order to add them in quadrature
                pass
            print "Significance:", sign
            sign = np.sqrt(sign)
            significances[i].append(sign)
            pass

        pass

    # Fix for icut == 35, isig == 2 # @TEMP!!!
    significances[2][35] = 0.5 * (significances[2][34] + significances[2][36])

    c_temp.xlabel("Large-#it{R} jet mass [GeV]")
    c_temp.logy()
    c_temp.legend()
    if args.save: c_temp.save("plots/cutoptimisation_massspectrum.pdf")

    # Compute significance improvements
    improvements     = [np.array(sign) / sign[-1] for sign in significances]
    improvements_avg = np.mean(improvements, axis=0)
    idx_max = np.argmax(improvements_avg)

    print "Optimal cut: %.2f (avg. improvement: %.3f)" % (cuts[idx_max], improvements_avg[idx_max])
    msk = bkg['tau21DDT'] < cuts[idx_max]
    print "Background efficiency: %.2f%%" % (np.sum(bkg['weight'][msk])/np.sum(bkg['weight']) * 100.)
    print "Signal efficiencies:"
    for name, DSID in zip(sig_names, sorted(list(set(sig['DSID'])))):
        msk_num = (sig['DSID'] == DSID) & (sig['tau21DDT'] < cuts[idx_max])
        msk_den = (sig['DSID'] == DSID)
        print "  %s: %.2f%%" % (name, np.sum(sig['weight'][msk_num])/np.sum(sig['weight'][msk_den]) * 100.)
        pass

    print ""
    cut_man = 0.50
    print "For manual cut: %.2f (avg. improvement: %.3f)" % (cut_man, float(improvements_avg[np.where(cuts == cut_man)]))
    msk = bkg['tau21DDT'] < cut_man
    print "Background efficiency: %.2f%%" % (np.sum(bkg['weight'][msk])/np.sum(bkg['weight']) * 100.)
    print "Signal efficiencies:"
    for name, DSID in zip(sig_names, sorted(list(set(sig['DSID'])))): 
        msk_num = (sig['DSID'] == DSID) & (sig['tau21DDT'] < cut_man)
        msk_den = (sig['DSID'] == DSID)
        print "  %s: %.2f%%" % (name, np.sum(sig['weight'][msk_num])/np.sum(sig['weight'][msk_den]) * 100.)
        pass


    # Plot improvements
    # --------------------------------------------------------------------------

    # Create canvas
    c = ap.canvas(batch=not args.show)

    bins = np.linspace(axis['tau21DDT'][1], axis['tau21DDT'][2], axis['tau21DDT'][0] + 1, endpoint=True)

    # Draw histograms
    msk = bkg['tau21DDT'] < cut
    h_bkg = c.hist(bkg['tau21DDT'][msk], bins=bins, weights=bkg['weight'][msk], fillcolor=ROOT.kAzure+7, name='Incl. #gamma MC', normalise=True)
    h_sigs = list()
    for idx, (DSID, name) in enumerate(zip(sorted(list(set(sig['DSID']))), sig_names)):
        msk = (sig['DSID'] == DSID) & (sig['tau21DDT'] < cut)
        h_sigs.append(c.hist(sig['tau21DDT'][msk], bins=bins, weights=sig['weight'][msk], linecolor=ROOT.kRed + idx, label=name, normalise=True))
        pass

    # Overlay
    o = ap.overlay(c, color=ROOT.kViolet)
    graphs = list()
    for idx, impr in enumerate(improvements):
        graphs.append(o.graph(impr, bins=cuts, display=None))
        pass
    gr = o.graph(improvements_avg, bins=cuts, linecolor=ROOT.kViolet, linewidth=2, markerstyle=0, option='L')
    o.padding(0.50)
    o.label("Average significance improvement")
    
    # Text
    c.padding(0.40)
    c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
            "ISR #gamma selection",
            ], qualifier='Simulation Internal')

    c.xlabel('Signal jet #tau_{21}^{DDT}')
    c.ylabel('Events (normalised)')

    # Lines
    o.line(cuts[idx_max], 0, cuts[idx_max], improvements_avg[idx_max])

    # Legend
    c.legend(width=0.28)

    # Show
    if args.save: c.save('plots/cutoptimisation.pdf')
    if args.show: c.show()


    # Write output to file
    f = ROOT.TFile('output/hists_isrgamma_cutoptimisation.root', 'RECREATE')
    h_bkg.SetName('h_tau21DDT_bkg')
    h_bkg.Write()
    for idx in range(len(h_sigs)):
        name = sig_names[idx]
        m = re.search('\(([0-9]+) GeV\)', name)
        h_sigs[idx].SetName('h_tau21DDT_sig%s' % m.group(1))
        h_sigs[idx].Write()
        graphs[idx].SetName('gr_improvements_sig%s' % m.group(1))
        graphs[idx].Write()
        pass
    gr.SetName('gr_improvements_avg')
    gr.Write()
    f.Write()
    f.Close()

    return


# Main function call.
if __name__ == '__main__':
   main()
   pass
