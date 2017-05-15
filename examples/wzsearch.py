#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for fitting final jet mass spectra in search of W/Z peak and producing publication-level plots.

@file:   wzsearch.py
@author: Andreas Søgaard
@date:   3 May 2017
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

parser = argparse.ArgumentParser(description='Perform W/Z test search.')

parser.add_argument('--window', dest='window', type=float,
                    default=None,
                    help='Width of excluded mass window (default: None)')
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

    signal_DSID = int("1%02d085" % (0 if args.window is None else args.window * 100))

    # Load data
    files = {
        'data': glob.glob(tf.config['base_path'] + 'objdef_data_*.root'),
        'bkg':  glob.glob(tf.config['base_path'] + 'objdef_TF_%d.root' % signal_DSID),
        'sig':  glob.glob(tf.config['base_path'] + 'objdef_MC_3054*.root'),
        'sfl':  glob.glob(tf.config['base_path'] + 'objdef_TF_%d_signalfail.root' % signal_DSID),
        }

    if 0 in map(len, files.values()):
        warning("No files found.")
        return

    data      = {key: loadData(files[key], tf.config['finaltree'],                               prefix=tf.config['prefix']) for key in files}   
    data_up   = {key: loadData(files[key], tf.config['finaltree'].replace('Nominal', 'TF_UP'),   prefix=tf.config['prefix']) for key in ['bkg']} 
    data_down = {key: loadData(files[key], tf.config['finaltree'].replace('Nominal', 'TF_DOWN'), prefix=tf.config['prefix']) for key in ['bkg']} 
    info      = {key: loadData(files[key], tf.config['outputtree'], stop=1)                      for key in files}
    
    # Rename variables
    # ...

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    # Only needs to be done for 'signal', for 'nominal' syst. variation
    data['sig'] = append_fields(data['sig'], 'DSID', np.zeros((data['sig'].size,)), dtypes=int)
    for idx in info['sig']['id']:    
        msk = (data['sig']['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
        DSID = info['sig']['DSID'][idx]  # Get DSID for this file
        data['sig']['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
        data['sig']['DSID']  [msk] = DSID        # Store DSID
        pass
    # @TODO: k-factors?
    data['sig']['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity

    # Check output.
    if data['data'].size == 0:
        warning("No data was loaded.")
        return 


    # Plotting pre-fit jet mass spectrum
    # ----------------------------------------------------------------
       
    bestfit_mu = None

    for mu, fit, prefit, subtract in zip([0, 1, 1, None], 
                                         [False, False, True,  False],
                                         [True,  True,  True,  False], 
                                         [True,  True,  False, True]):

        if not prefit:
            mu = bestfit_mu[0]
            pass


        # Plotting
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
        c = ap.canvas(num_pads=2, batch=not args.show)
        p0, p1 = c.pads()

        # -- Histograms: Main pad
        bins = np.linspace(50, 110, 30 + 1, endpoint=True)#tf.config['massbins']
        
        h_bkg      = c.hist(data     ['bkg']['m'], bins=bins, weights=data     ['bkg']['weight'], display=False)
        h_bkg_up   = c.hist(data_up  ['bkg']['m'], bins=bins, weights=data_up  ['bkg']['weight'], display=False)
        h_bkg_down = c.hist(data_down['bkg']['m'], bins=bins, weights=data_down['bkg']['weight'], display=False)

        h_sig  = c.hist(data['sig']['m'], bins=bins, weights=data['sig']['weight'], scale=mu, display=False)
        h_sfl  = c.hist(data['sfl']['m'], bins=bins, weights=data['sfl']['weight'], scale=mu, display=False)
        
        if not fit:
            h_bkg     .Add(h_sfl, -1) # Subtracting signal
            h_bkg_up  .Add(h_sfl, -1) # --
            h_bkg_down.Add(h_sfl, -1) # --
            pass

        
        h_bkg  = c.stack(h_bkg,
                         fillcolor=ROOT.kAzure + 7, 
                         label='Background pred.')
        h_sig  = c.stack(h_sig,
                         fillcolor=ROOT.kRed + 1,
                         label="W/Z (#mu = %s)" % ("%.0f" % mu if prefit else "%.2f #pm %.2f" % (mu, bestfit_mu[1])))
        
        h_sum = h_bkg#c.getStackSum()
        h_sum = c.hist(h_sum, 
                       fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                       label='Stat. uncert.')

        h_bkg_up   = c.hist(h_bkg_up,
                            linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST',
                            label='Syst. uncert.')
        h_bkg_down = c.hist(h_bkg_down,
                            linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST')
    
        h_data = c.plot(data['data']['m'], bins=bins, weights=data['data']['weight'],
                        label='Data')

        # -- Histograms: Ratio pad
        c.ratio_plot((h_sig,      h_sum), option='HIST', offset=1, fillcolor=ROOT.kRed + 1)
        c.ratio_plot((h_sum,      h_sum), option='E2')
        c.ratio_plot((h_bkg_up,   h_sum), option='HIST')
        c.ratio_plot((h_bkg_down, h_sum), option='HIST')
        c.ratio_plot((h_data,     h_sum))

        # -- Text
        c.text(["#sqrt{s} = 13 TeV,  L = %s fb^{-1}" % tf.config['lumi'],
                "Trimmed anti-k_{t}^{R=1.0} jets",
                "ISR #gamma selection",
                "Window: %d GeV #pm %d %%" % (85, (0.2 if args.window is None else args.window) * 100.), # @TODO: Improve
               ], qualifier='Internal')

        # -- Axis labels
        c.xlabel('Signal jet mass [GeV]')
        c.ylabel('Events')
        p1.ylabel('Data / Est.')

        # -- Axis limits
        #c.ylim(2E+01, 2E+06)
        #c.ylim(0, 14000)
        p1.ylim(0.8, 1.2)

        # -- Line(s)
        p1.yline(1.0)
        for d in [0.8, 1.2] + ([] if args.window is None else [1 - args.window, 1 + args.window]):
            p0.xline(d * 85.)
            p1.xline(d * 85.)
            pass
        
        #c.log()            
        c.legend()
        if args.show and not fit: c.show()
        if args.save and not fit: c.save('plots/wzsearch_%s%s.pdf' % ('' if args.window is None else 'pm%d_' % (args.window * 100), 'prefit_mu%d' % mu if prefit else 'postfit'))

        
        # Fit
        if fit:

            hs_save = [
                    h_bkg_down.Clone('h_save_down'),
                    h_bkg     .Clone('h_save_nom'),
                    h_bkg_up  .Clone('h_save_up'),
                    ]

            bestfit_mu = list()

            for variation in range(3):

                print "Variation: " + ("Nominal" if variation == 1 else ("Up" if variation == 0 else "Down"))

                # Get correct histogram for this variation
                h_bkg_use = hs_save[variation]
                
                # -- Define jet mass variable
                mJ = ROOT.RooRealVar('mJ', 'mJ', 50, 300)
                mJ.setBins(50)
                
                # -- Define histograms
                rdh_bkg = ROOT.RooDataHist('rdh_bkg', 'rdh_bkg', ROOT.RooArgList(mJ), h_bkg_use)
                rdh_sig = ROOT.RooDataHist('rdh_sig', 'rdh_sig', ROOT.RooArgList(mJ), h_sig)
                rdh_sfl = ROOT.RooDataHist('rdh_sfl', 'rdh_sfl', ROOT.RooArgList(mJ), h_sfl)
                
                # -- Turn histograms into pdf's
                rhp_bkg = ROOT.RooHistPdf('rhp_bkg', 'rhp_bkg', ROOT.RooArgSet(mJ), rdh_bkg)
                rhp_sig = ROOT.RooHistPdf('rhp_sig', 'rhp_sig', ROOT.RooArgSet(mJ), rdh_sig)
                rhp_sfl = ROOT.RooHistPdf('rhp_sfl', 'rhp_sfl', ROOT.RooArgSet(mJ), rdh_sfl)
                
                # -- Define integrals as constants
                n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', h_bkg_use.Integral())
                n_sig = ROOT.RooRealVar('n_sig', 'n_sig', h_sig.Integral())
                n_sfl = ROOT.RooRealVar('n_sfl', 'n_sfl', h_sfl.Integral())
                
                # -- Define signal strength and constant(s)
                mu   = ROOT.RooRealVar('mu',   'mu', 1, 0, 5)
                neg1 = ROOT.RooRealVar('neg1', 'neg1', -1)
                
                # -- Define fittable normalisation factors
                c_bkg = ROOT.RooFormulaVar('c_bkg', 'c_bkg', '@0',           ROOT.RooArgList(n_bkg))
                c_sig = ROOT.RooFormulaVar('c_sig', 'c_sig', '@0 * @1',      ROOT.RooArgList(mu, n_sig))
                c_sfl = ROOT.RooFormulaVar('c_sfl', 'c_sfl', '@0 * @1 * @2', ROOT.RooArgList(neg1, mu, n_sfl))
                
                # -- Construct combined pdf
                pdf = ROOT.RooAddPdf('pdf', 'pdf', ROOT.RooArgList(rhp_bkg, rhp_sig, rhp_sfl), ROOT.RooArgList(c_bkg, c_sig, c_sfl))
                
                # -- Construct data histogram
                rdh_data = ROOT.RooDataHist('rdh_data', 'rdh_data', ROOT.RooArgList(mJ), h_data)
                
                # -- Fit pdf to data histogram
                pdf.chi2FitTo(rdh_data, ROOT.RooLinkedList())
                
                print "Best fit mu: %.3f +/- %.3f" % (mu.getValV(), mu.getError())
                bestfit_mu.append( (mu.getValV(), mu.getError()) )
                pass

            bestfit_mu = bestfit_mu[1][0], np.sqrt(np.power(abs(bestfit_mu[0][0] - bestfit_mu[2][0]) / 2., 2.) + np.power(bestfit_mu[1][1], 2.))
            pass

        pass
    
    return


# Main function call.
if __name__ == '__main__':
   main()
   pass
