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

    # Set correct weighting
    ROOT.TH1.SetDefaultSumw2(True)


    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    signal_DSID = int("1%02d085" % (0 if args.window is None else args.window * 100))

    # Load data
    files = {
        'data': glob.glob(tf.config['base_path'] + 'objdef_data_*.root'),
        'bkg':  glob.glob(tf.config['base_path'] + 'objdef_TF_%d.root' % signal_DSID),
        #'sig':  glob.glob(tf.config['base_path'] + 'objdef_MC_3054*.root'),
        'W':    glob.glob(tf.config['base_path'] + 'objdef_MC_30543*.root'),
        'Z':    glob.glob(tf.config['base_path'] + 'objdef_MC_30544*.root'),
        'sfl':  glob.glob(tf.config['base_path'] + 'objdef_TF_%d_signalfail.root' % signal_DSID),
        }

    if 0 in map(len, files.values()):
        warning("No files found.")
        return
                 
    data      = {key: loadData(files[key], tf.config['finaltree'],                               prefix=tf.config['prefix']) for key in files}   
    data_up   = {key: loadData(files[key], tf.config['finaltree'].replace('Nominal', 'TF_UP'),   prefix=tf.config['prefix']) for key in ['bkg', 'sfl']} 
    data_down = {key: loadData(files[key], tf.config['finaltree'].replace('Nominal', 'TF_DOWN'), prefix=tf.config['prefix']) for key in ['bkg', 'sfl']} 
    info      = {key: loadData(files[key], tf.config['outputtree'], stop=1)                      for key in files}
    
    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    # Only needs to be done for 'signal', for 'nominal' syst. variation
    for comp in ['W', 'Z']:  # 'sig'
        data[comp] = append_fields(data[comp], 'DSID', np.zeros((data[comp].size,)), dtypes=int)
        for idx in info[comp]['id']:    
            msk = (data[comp]['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            DSID = info[comp]['DSID'][idx]  # Get DSID for this file
            data[comp]['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
            data[comp]['DSID']  [msk] = DSID        # Store DSID
            pass
        data[comp]['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass

    # Check output.
    if data['data'].size == 0:
        warning("No data was loaded.")
        return 

    # k-factors
    kfactor = { 'W': 0.7838,
                'Z': 1.3160 }

    data['W']['weight'] *= kfactor['W']
    data['Z']['weight'] *= kfactor['Z']
    data['sig'] = np.concatenate((data['W'], data['Z']))


    # Plotting pre-fit jet mass spectrum
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
       
    bestfit_mu = None

    for mu, fit, prefit, subtract in zip([0, 1, 1, None], 
                                         [False, False, True,  False],
                                         [True,  True,  True,  False], 
                                         [True,  True,  False, True]):

        if not prefit:
            print "#"* 20 , bestfit_mu, "#"*20
            mu = bestfit_mu[0]
            pass


        # Plotting
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        
        c = ap.canvas(num_pads=2, batch=not args.show)
        p0, p1 = c.pads()

        # -- Histograms: Main pad
        bins = np.linspace(50, 112, 31 + 1, endpoint=True) # tf.config['massbins']

        h_bkg      = c.hist(data     ['bkg']['m'], bins=bins, weights=data     ['bkg']['weight'], display=False)
        h_bkg_up   = c.hist(data_up  ['bkg']['m'], bins=bins, weights=data_up  ['bkg']['weight'], display=False)
        h_bkg_down = c.hist(data_down['bkg']['m'], bins=bins, weights=data_down['bkg']['weight'], display=False)

        h_sig  = c.hist(data['sig']['m'], bins=bins, weights=data['sig']['weight'],
                        scale=mu,
                        display=False)

        h_sfl      = c.hist(data     ['sfl']['m'], bins=bins, weights=data     ['sfl']['weight'], scale=mu, display=False)
        h_sfl_up   = c.hist(data_up  ['sfl']['m'], bins=bins, weights=data_up  ['sfl']['weight'], scale=mu, display=False)
        h_sfl_down = c.hist(data_down['sfl']['m'], bins=bins, weights=data_down['sfl']['weight'], scale=mu, display=False)
        
        if not fit:
            h_bkg     .Add(h_sfl,      -1) # Subtracting signal
            h_bkg_up  .Add(h_sfl_up,   -1) # --
            h_bkg_down.Add(h_sfl_down, -1) # --
            pass

        
        h_bkg  = c.stack(h_bkg,
                         fillcolor=ROOT.kAzure + 7, 
                         label='Background pred.')
        
        h_Z = c.stack(data['Z']['m'], bins=bins, weights=data['Z']['weight'],
                      scale=mu,
                      fillcolor=ROOT.kAzure + 3,
                      label="Z + #gamma  (#times #mu)")
        h_W = c.stack(data['W']['m'], bins=bins, weights=data['W']['weight'],
                      scale=mu,
                      fillcolor=ROOT.kAzure + 2,
                      label="W + #gamma (#times #mu)")
        
        h_sum = h_bkg
        h_sum = c.hist(h_sum, 
                       fillstyle=3245, fillcolor=ROOT.kGray + 2, linecolor=ROOT.kGray + 3, option='E2',
                       label='Stat. uncert.')

        h_bkg_up   = c.hist(h_bkg_up,
                            linecolor=ROOT.kRed + 1, linestyle=2, option='HIST',
                            label='TF syst. uncert.')
        h_bkg_down = c.hist(h_bkg_down,
                            linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST')
    
        h_data = c.plot(data['data']['m'], bins=bins, weights=data['data']['weight'],
                        label='Data')


        # Manuall check chi2 for nominal systematic
        """
        if not fit:
            chi2 = 0.
            chi2_data = 0.
            for ibin in range(1, h_data.GetXaxis().GetNbins() + 1):
                stat_data = h_data.GetBinError(ibin)
                cont_data = h_data.GetBinContent(ibin)
                cont_bkg  = h_bkg .GetBinContent(ibin)
                cont_sig  = h_sig .GetBinContent(ibin)
                diff = np.abs(cont_data - (cont_bkg + cont_sig))
                chi2 +=  np.square(diff / stat_data)
                pass
            print "[INFO] mu = {:.2f} | chi2 = {:.1f}".format(mu, chi2)
            pass
            """

        # -- Histograms: Ratio pad
        c.ratio_plot((h_sig,      h_sum), option='HIST', offset=1, fillcolor=ROOT.kAzure + 2)
        c.ratio_plot((h_Z,        h_sum), option='HIST', offset=1, fillcolor=ROOT.kAzure + 3)
        c.ratio_plot((h_sum,      h_sum), option='E2')
        c.ratio_plot((h_bkg_up,   h_sum), option='HIST')
        c.ratio_plot((h_bkg_down, h_sum), option='HIST')
        c.ratio_plot((h_data,     h_sum))

        # -- Text
        c.text(["#sqrt{s} = 13 TeV,  L = %s fb^{-1}" % tf.config['lumi'],
                "Trimmed anti-k_{t}^{R=1.0} jets",
                "ISR #gamma selection",
                "Window: %d GeV #pm %d %%" % (85, (0.2 if args.window is None else args.window) * 100.), # @TODO: Improve
                "%s-fit: #mu = %s" % ("Pre" if prefit else "Post", "%.0f" % mu if prefit else "%.2f #pm %.2f" % (mu, bestfit_mu[1])),
               ], qualifier='Internal')

        # -- Axis labels
        c.xlabel('Signal jet mass [GeV]')
        c.ylabel('Events')
        p1.ylabel('Data / Est.')

        # -- Axis limits
        c.padding(0.45)

        p1.ylim(0.95, 1.15)

        # -- Line(s)
        p1.yline(1.0)
        
        # -- Region labels
        c.region("SR", 0.8 * 85, 1.2 * 85)
        if args.window is not None:
            c.region("VR", (1 - args.window) * 85, 0.8 * 85)
            c.region("VR", 1.2 * 85, (1 + args.window) * 85)
            pass
        
        #c.log()            
        c.legend()
        if args.show and not fit: c.show()
        if args.save and not fit: c.save('plots/wzsearch_%s%s.pdf' % ('' if args.window is None else 'pm%d_' % (args.window * 100), 'prefit_mu%d' % mu if prefit else 'postfit'))

        
        # Saving output for harmonised paper plots
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if not prefit:
            outfile = ROOT.TFile('output/hists_isrgamma_WZ.root', 'RECREATE')
            h_bkg.SetName('h_mJ_QCD')
            h_bkg_up.SetName('h_mJ_QCD_up')
            h_bkg_down.SetName('h_mJ_QCD_down')
            h_data.SetName('h_mJ_data')
            h_W.SetName('h_mJ_WHad')
            h_Z.SetName('h_mJ_ZHad')
            for hist in [h_bkg, h_bkg_up, h_bkg_down, h_data, h_W, h_Z]:
                hist.Write()
                pass
            outfile.Close()
            pass

        # Fitting
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

        if fit:

            hs_save = [
                    h_bkg_down.Clone('h_save_down'),
                    h_bkg     .Clone('h_save_nom'),
                    h_bkg_up  .Clone('h_save_up'),
                    ]

            hdata_save = [
                    h_data.Clone('hdata_save_down'),
                    h_data.Clone('hdata_save_nom'),
                    h_data.Clone('hdata_save_up'),
                    ]

            hsfl_save = [
                    h_sfl_down.Clone('h_sfl_down'),
                    h_sfl     .Clone('h_sfl_nom'),
                    h_sfl_up  .Clone('h_sfl_up'),
                    ]

            bestfit_mu = list()

            for variation in range(3):

                print "Variation: " + ("Nominal" if variation == 1 else ("Down" if variation == 0 else "Up"))

                # Get correct histogram for this variation
                h_bkg_use  = hs_save   [variation]
                h_sfl_use  = hsfl_save [variation]
                h_data_use = hdata_save[variation]

                # -- Define jet mass variable
                mJ = ROOT.RooRealVar('mJ', 'mJ', 50, 112)
                mJ.setBins(31)
                
                # -- Define histograms
                rdh_bkg = ROOT.RooDataHist('rdh_bkg', 'rdh_bkg', ROOT.RooArgList(mJ), h_bkg_use)
                rdh_sig = ROOT.RooDataHist('rdh_sig', 'rdh_sig', ROOT.RooArgList(mJ), h_sig)
                rdh_sfl = ROOT.RooDataHist('rdh_sfl', 'rdh_sfl', ROOT.RooArgList(mJ), h_sfl_use)
                
                # -- Turn histograms into pdf's
                rhp_bkg = ROOT.RooHistPdf('rhp_bkg', 'rhp_bkg', ROOT.RooArgSet(mJ), rdh_bkg)
                rhp_sig = ROOT.RooHistPdf('rhp_sig', 'rhp_sig', ROOT.RooArgSet(mJ), rdh_sig)
                rhp_sfl = ROOT.RooHistPdf('rhp_sfl', 'rhp_sfl', ROOT.RooArgSet(mJ), rdh_sfl)
                
                # -- Define integrals as constants
                n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', h_bkg_use.Integral())
                n_sig = ROOT.RooRealVar('n_sig', 'n_sig', h_sig    .Integral())
                n_sfl = ROOT.RooRealVar('n_sfl', 'n_sfl', h_sfl_use.Integral())
                
                # -- Define signal strength and constant(s)
                v_mu   = ROOT.RooRealVar('v_mu',   'v_mu', 1, 0, 5)
                v_neg1 = ROOT.RooRealVar('v_neg1', 'v_neg1', -1)
                
                # -- Define fittable normalisation factors
                c_bkg = ROOT.RooFormulaVar('c_bkg', 'c_bkg', '@0',           ROOT.RooArgList(n_bkg))
                c_sig = ROOT.RooFormulaVar('c_sig', 'c_sig', '@0 * @1',      ROOT.RooArgList(v_mu, n_sig))
                c_sfl = ROOT.RooFormulaVar('c_sfl', 'c_sfl', '@0 * @1 * @2', ROOT.RooArgList(v_neg1, v_mu, n_sfl))
                
                # -- Construct combined pdf
                pdf = ROOT.RooAddPdf('pdf', 'pdf', ROOT.RooArgList(rhp_bkg, rhp_sig, rhp_sfl), ROOT.RooArgList(c_bkg, c_sig, c_sfl))
                
                # -- Construct data histogram
                rdh_data = ROOT.RooDataHist('rdh_data', 'rdh_data', ROOT.RooArgList(mJ), h_data_use)
                
                # -- Fit pdf to data histogram
                #pdf.fitTo(rdh_data)
                chi2Fit = ROOT.RooChi2Var("chi2","chi2", pdf, rdh_data, ROOT.RooFit.Extended())
                minuit = ROOT.RooMinuit (chi2Fit)
                minuit.migrad()
                minuit.hesse()
                r2 = minuit.save()

                bestfit_mu.append( (v_mu.getValV(), v_mu.getError()) )
                pass

            print "("* 10, bestfit_mu, ")" * 10
            bestfit_mu = bestfit_mu[1][0], np.sqrt(np.power(abs(bestfit_mu[0][0] - bestfit_mu[2][0]) / 2., 2.) + np.power(bestfit_mu[1][1], 2.))
            pass

        pass
    
    return


# Main function call.
if __name__ == '__main__':
   main()
   pass
