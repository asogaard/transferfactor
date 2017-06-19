#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for plotting unblinded data distributions and producing publication-level plots.

@file:   datadistributions.py
@author: Andreas Søgaard
@date:   2 June 2017
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
    from transferfactor.utils import *
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

parser.add_argument('--mass', dest='mass', type=float,
                    required=True,
                    help='Center of excluded mass window (default: None)')
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

    TF_DSID = int("100%03d" % (args.mass))
    signal_DSID = get_signal_DSID(args.mass)
    signal = bool(signal_DSID)
    #if not signal_DSID:
    #    return

    # Load data
    files = {
        'data': glob.glob(tf.config['base_path'] + 'objdef_data_*.root'),
        'bkg':  glob.glob(tf.config['base_path'] + 'objdef_TF_%d.root' % TF_DSID),
        'gbs':  glob.glob(tf.config['base_path'] + 'objdef_GBS_400000.root'),
        'sig':  glob.glob(tf.config['base_path'] + 'objdef_MC_%d.root' % signal_DSID) if signal else [],
        'W':    glob.glob(tf.config['base_path'] + 'objdef_MC_30543*.root'),
        'Z':    glob.glob(tf.config['base_path'] + 'objdef_MC_30544*.root'),
        'sfl':  glob.glob(tf.config['base_path'] + 'objdef_TF_%d_signalfail.root' % TF_DSID),
        }
    
    if 0 in map(len, files.values()) and signal:
        warning("No files found.")
        return
    
    data      = {key: loadData(files[key], tf.config['finaltree'],                               prefix=tf.config['prefix']) for key in files}   
    data_up   = {key: loadData(files[key], tf.config['finaltree'].replace('Nominal', 'TF_UP'),   prefix=tf.config['prefix']) for key in ['bkg', 'sfl', 'gbs']} 
    data_down = {key: loadData(files[key], tf.config['finaltree'].replace('Nominal', 'TF_DOWN'), prefix=tf.config['prefix']) for key in ['bkg', 'sfl', 'gbs']} 
    info      = {key: loadData(files[key], tf.config['outputtree'], stop=1)                      for key in files}
    
    
    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    # Only needs to be done for 'signal', for 'nominal' syst. variation
    for comp in ['sig', 'W', 'Z']:
        if comp == 'sig' and not signal: continue
        data[comp] = append_fields(data[comp], 'DSID', np.zeros((data[comp].size,)), dtypes=int)
        for idx in info[comp]['id']:    
            msk = (data[comp]['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
            DSID = info[comp]['DSID'][idx]  # Get DSID for this file
            data[comp]['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
            data[comp]['DSID']  [msk] = DSID        # Store DSID
            pass
        # @TODO: k-factors?
        pass
    
    # Failing signal component is scaled by cross section, but not luminosity
    for comp in ['sig', 'W', 'Z']: # 'sfl'
        if comp in ['sig', 'sfl'] and not signal: continue
        data[comp]['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
        pass
    
    # Check output.
    if data['data'].size == 0:
        warning("No data was loaded.")
        return 
    
    
    # Plotting
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    for mu, gbs, log in itertools.product([None, 0, 0.5, 1.], [False, True], [False, True]):
        
        signal = bool(signal_DSID) and bool(mu is not None)
        
        c = ap.canvas(num_pads=2, batch=not args.show)
        p0, p1 = c.pads()
        
        # -- Histograms: Main pad
        bins = tf.config['massbins']
        
        bkgvar = 'gbs' if gbs else 'bkg'
        h_bkg      = c.hist(data     [bkgvar]['m'], bins=bins, weights=data     [bkgvar]['weight'], display=False)
        h_bkg_up   = c.hist(data_up  [bkgvar]['m'], bins=bins, weights=data_up  [bkgvar]['weight'], display=False)
        h_bkg_down = c.hist(data_down[bkgvar]['m'], bins=bins, weights=data_down[bkgvar]['weight'], display=False)
        
        # -- Subtract failing signal component (if any)
        if signal:
            h_sig  = c.hist(data['sig']['m'], bins=bins, weights=data['sig']['weight'], scale=mu, display=False)
            pass


        h_Z = c.stack(data['Z']['m'], bins=bins, weights=data['Z']['weight'],
                      fillcolor=ROOT.kAzure + 3, 
                      label="Z + #gamma")
        
        h_W = c.stack(data['W']['m'], bins=bins, weights=data['W']['weight'],
                      fillcolor=ROOT.kAzure + 2, 
                      label="W + #gamma")
        
        # ---------------------------------------------------------------------
        # Saving output for harmonised paper plots
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if (not gbs) and log:
            outfile = ROOT.TFile('output/hists_isrgamma_datadistributions_%dGeV%s.root' % (args.mass, '' if mu is None else '_mu%d' % (100 * mu)), 'RECREATE')
            h_bkg.SetName('h_mJ_QCD')
            h_bkg_up.SetName('h_mJ_QCD_up')
            h_bkg_down.SetName('h_mJ_QCD_down')
            h_data.SetName('h_mJ_data')
            h_W.SetName('h_mJ_WHad')
            h_Z.SetName('h_mJ_ZHad')
            if signal:
                h_sig.SetName('h_mJ_sig')
                pass
            for hist in [h_bkg, h_bkg_up, h_bkg_down, h_data, h_W, h_Z] + ([h_sig] if signal else []):
                hist.Write()
                pass
            outfile.Close()
            pass
        # ---------------------------------------------------------------------        

        # -- Subtract failing signal component (if any)
        if signal:
            h_sfl      = c.hist(data     ['sfl']['m'], bins=bins, weights=data     ['sfl']['weight'], scale=mu, display=False)
            h_sfl_up   = c.hist(data_up  ['sfl']['m'], bins=bins, weights=data_up  ['sfl']['weight'], scale=mu, display=False)
            h_sfl_down = c.hist(data_down['sfl']['m'], bins=bins, weights=data_down['sfl']['weight'], scale=mu, display=False)
            
            h_bkg     .Add(h_sfl, -1) # Subtracting signal
            h_bkg_up  .Add(h_sfl, -1) # --
            h_bkg_down.Add(h_sfl, -1) # --
            pass                    
        
        h_bkg = c.stack(h_bkg,
                        fillcolor=ROOT.kAzure + 7, 
                        label='Bkg. pred. (GBS)' if gbs else 'Background pred.')

        
        h_sum = c.getStackSum()
        
        if signal:
            h_sig = c.stack(h_sig,
                            fillcolor=ROOT.kRed - 4,
                            label="Z' (%d GeV) + #gamma" % args.mass)#  (#times #mu)")
            pass
        
        h_sum = c.hist(h_sum, 
                       fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                       label='Stat. uncert.')

        # Add h_sum errors in quadrature to TF systematics (?)
        if True:
            for idx in range(1, h_sum.GetXaxis().GetNbins() + 1):
                nom       = h_bkg.GetBinContent(idx)
                s         = h_sum.GetBinContent(idx)
                stat      = h_sum.GetBinError  (idx)
                syst_up   = np.abs(h_bkg_up  .GetBinContent(idx) - nom)
                syst_down = np.abs(h_bkg_down.GetBinContent(idx) - nom)
                h_bkg_up  .SetBinContent(idx, s + np.sqrt(np.square(stat) + np.square(syst_up)))
                h_bkg_down.SetBinContent(idx, s - np.sqrt(np.square(stat) + np.square(syst_down)))
                pass
            pass
        
        h_bkg_up   = c.hist(h_bkg_up,
                            #linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST',
                            linecolor=ROOT.kGray + 3, option='HIST', display=False)
                            #label='TF syst. uncert.') 
                            #label='Stat. #oplus syst.')
        c.hist(h_sum, linecolor=ROOT.kGray + 3, fillstyle=0, option='HIST',  label='Stat. #oplus syst.')
        c.hist(h_sum, linecolor=ROOT.kBlack,    fillstyle=0, option='HIST')
        
        h_bkg_down = c.hist(h_bkg_down,
                            #linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST')
                            linecolor=ROOT.kGray + 3, option='HIST', display=False)
        
        h_data = c.plot(data['data']['m'], bins=bins, weights=data['data']['weight'],
                        label='Data')
        
        # -- Axis limits
        #c.padding(0.40) # 0.45
        p1.ylim(0.8, 1.2)
        
        # -- Histograms: Ratio pad
        if signal:
            c.ratio_plot((h_sig,      h_sum), option='HIST', offset=1, fillcolor=ROOT.kRed - 4)
            pass
        c.ratio_plot((h_sum,      h_sum), option='E2')
        c.ratio_plot((h_bkg_up,   h_sum), option='HIST')
        c.ratio_plot((h_bkg_down, h_sum), option='HIST')
        c.ratio_plot((h_data,     h_sum), oob=True)
        
        # -- Text
        c.text(["#sqrt{s} = 13 TeV,  L = %s fb^{-1}" % tf.config['lumi'],
                "Trimmed anti-k_{t}^{R=1.0} jets",
                "ISR #gamma selection",
                #"Window: %d GeV #pm %d %%" % (args.mass, 20.),
                ] + (["Pre-fit: #mu = %.2f" % mu] if signal else []), 
               qualifier='Internal')
        
        # -- Axis labels
        c.xlabel('Signal jet mass [GeV]')
        c.ylabel('Events')
        p1.ylabel('Data / Est.')
        
        # -- Line(s)
        p1.yline(1.0)
        
        # -- Region labels
        c.region("SR", 0.8 * args.mass, 1.2 * args.mass)
        
        c.log(log)
        c.legend()
        if args.show: c.show()
        if args.save: c.save('plots/datadistribution_%dGeV%s%s%s.pdf' % (args.mass, '_mu%d' % (mu * 100)if signal else '', '_gbs' if gbs else '', '_log' if log else ''))
        pass

    return


# Main function call.
if __name__ == '__main__':
   main()
   pass
