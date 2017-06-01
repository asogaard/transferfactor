#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for computing the global background shape (GBS) for the dijet + ISR analysis.

@file:   globalbackground.py
@author: Andreas Søgaard 
@date:   1 May 2017
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
    from transferfactor.utils import get_signal_DSID
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

parser = argparse.ArgumentParser(description='Compute global background shape (GBS).')

parser.add_argument('--mass', dest='mass', type=int,
                    required=True,
                    help='Center of excluded mass window')
parser.add_argument('--inject', dest='inject', action='store_const',
                    const=True, default=False,
                    help='Inject signal (default: False)')
parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')


# Main function
def main ():
    
    # Parse command-line arguments
    args = parser.parse_args()


    # Setup.
    # – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – 

    # Get signal file
    sig_DSID = get_signal_DSID(args.mass, tolerance=10)
    if sig_DSID is None:
        warning ("No signal file was found")
        return
    sig_file = 'objdef_MC_{DSID:6d}.root'.format(DSID=sig_DSID)

    # Load data
    files = glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root')
    if args.inject:
        files += [tf.config['base_path'] + sig_file]
        pass

    if len(files) == 0:
        warning("No files found. Try to run:")
        warning(" $ source getSomeData.sh")
        return

    data = loadData(files, tf.config['tree'], prefix=tf.config['prefix'])
    info = loadData(files, tf.config['outputtree'], stop=1)

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    lumi = 36.1

    # Append new DSID field
    data = append_fields(data, 'DSID', np.zeros((data.size,)), dtypes=int)
    for idx in info['id']:    
        msk = (data['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
        DSID = info['DSID'][idx]  # Get DSID for this file
        data['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
        data['DSID']  [msk] = DSID # Store DSID
        pass
    data['weight'] *= lumi # Scale all events (MC) by luminosity
    # k-factors?
    data = append_fields(data, 'logpt', np.log(data['pt']))
    
    msk_sig  = (data['DSID'] == sig_DSID)
    msk_data = ~msk_sig

    signal = data[msk_sig]

 
    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Pass/fail masks
    # -- Data (incl. signal)
    msk_pass = tf.config['pass'](data)
    msk_fail = ~msk_pass

    # -- Signal
    msk_sig_pass = tf.config['pass'](signal)

    # Transfer factor calculator instance
    calc = tf.calculator(data=data, config=tf.config)
    
    # Nominal fit
    calc.fit()
    w_nom = calc.weights(data[msk_fail])
    if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/globalbackground_%s_' % ('injected' if args.inject else 'notinjected'))

    # mass +/- 20% stripe fit
    calc.mass   = args.mass
    calc.window = 0.2
    calc.fit()
    w_stripe = calc.weights(data[msk_fail])
    
    
    # Get GBS
    bins = tf.config['massbins']

    masses = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220]
    dict_backgrounds = dict()
    ctemp = ap.canvas(batch=True)
    for mass in masses:
        print " --", mass

        # Fit TF profile
        calc.mass = mass
        calc.fit()
        if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/globalbackground_%dGeV_%s_' % (args.mass, 'injected' if args.inject else 'notinjected'))

        # Get TF weights
        w = calc.weights(data[msk_fail])

        # Compute background distribution for this mass point
        dict_backgrounds[mass] = hist2array( ctemp.hist(data['m'][msk_fail], bins=bins, weights=data['weight'][msk_fail] * w, display=False) )

        # Store SR bin values (bin center within 20% if mass point)
        bincentres = bins[:-1] + (bins[1] - bins[0]) / 2.
        dict_backgrounds[mass] = np.ma.array(dict_backgrounds[mass], mask=~(np.abs(bincentres - mass) < 0.2 * mass))
        pass

    # Convert backgrounds from dict to numpy.ma.array
    backgrounds = np.ma.zeros((len(masses), len(bins) - 1))
    for idx, mass in enumerate(masses):
        backgrounds[idx,:] = dict_backgrounds[mass]
        pass

    
    # Plotting
    # – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – 
    
    # Setup canvas
    c = ap.canvas(num_pads=2, batch=not args.show)
    p0, p1 = c.pads()

    # Add stacked backgrounds
    h_bkg_nom    = c.stack(data  ['m'][msk_fail],     bins=bins, weights=data  ['weight'][msk_fail] * w_nom,    fillcolor=ROOT.kAzure  + 7, label="Bkg. (Nom.)")
    if args.inject:
        h_sig    = c.stack(signal['m'][msk_sig_pass], bins=bins, weights=signal['weight'][msk_sig_pass],        fillcolor=ROOT.kOrange + 1, label="Z' (%d GeV)" % args.mass)
        pass
    h_bkg_stripe = c.hist (data  ['m'][msk_fail],     bins=bins, weights=data  ['weight'][msk_fail] * w_stripe, linecolor=ROOT.kGreen  + 1, label="Bkg. (%d GeV #pm 20%%)" % args.mass)
    h_gbs        = c.hist(backgrounds.mean(axis=0),   bins=bins,                                                linecolor=ROOT.kViolet + 1, label="Bkg. (GBS)")
    
    # Draw stats. error of stacked sum
    h_sum = h_bkg_nom.Clone('h_sum') #c.getStackSum()
    h_sum = c.hist(h_sum, fillstyle=3245, fillcolor=ROOT.kGray+2, linecolor=ROOT.kGray + 3, label='Stats. uncert.', option='E2')
    
    # Add (pseudo-) data
    np.random.seed(21)
    h_data = c.plot (data['m'][msk_pass], bins=bins, weights=data['weight'][msk_pass], markersize=0.8, label='Pseudo-data', scale=1)

    # Draw error- and ratio plots
    if args.inject:
        hr_sig  = c.ratio_plot((h_sig,       h_bkg_nom), option='HIST', offset=1)
        pass
    h_err   = c.ratio_plot((h_sum,       h_bkg_nom), option='E2')
    h_ratio = c.ratio_plot((h_data,      h_bkg_nom), markersize=0.8)
    h_rgbs  = c.ratio_plot((h_gbs,       h_bkg_nom), linecolor=ROOT.kViolet + 1, option='HIST ][')
    h_rgbs  = c.ratio_plot((h_bkg_stripe,h_bkg_nom), linecolor=ROOT.kGreen  + 1, option='HIST ][')
    
    # Add labels and text
    c.xlabel('Signal jet mass [GeV]')
    c.ylabel('Events')
    p1.ylabel('Data / Nom.')
    c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
            "Sherpa incl. #gamma MC",
            "Trimmed anti-k_{t}^{R=1.0} jets",
            "ISR #gamma selection",
            ("Signal" if args.inject else "No signal") + " injected",
            ], 
           qualifier='Simulation Internal')

    # Add line(s) and limit(s)
    c.ylim(5E+01, 5E+06) # @TEMP
    p1.yline(1.0)
    p1.ylim(0.8, 1.2)

    # Configure axis padding and -scale
    c.padding(0.35)
    c.log(True)

    # Draw legend
    c.legend()

    # Save and show plot
    if args.save: c.save('plots/globalbackground_spectrum_%dGeV_%s.pdf' % (args.mass, 'injected' if args.inject else 'notinjected'))
    if args.show: c.show()


    # p0-plot
    # – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – 
    
    # Setup canvas
    c2 = ap.canvas(batch=not args.show)

    p_local  = h_data.Clone('p_local')
    p_global = h_data.Clone('p_global')

    for bin in range(1, h_data.GetXaxis().GetNbins() + 1):
        c_data = h_data      .GetBinContent(bin)
        e_data = h_data      .GetBinError  (bin)
        c_loc  = h_bkg_stripe.GetBinContent(bin)
        e_loc  = h_bkg_stripe.GetBinError  (bin)
        c_glb  = h_gbs       .GetBinContent(bin)
        e_glb  = e_loc # h_gbs    .GetBinError  (bin)

        z_loc = (c_data - c_loc) / np.sqrt( np.square(e_data) + np.square(e_loc) )
        z_glb = (c_data - c_glb) / np.sqrt( np.square(e_data) + np.square(e_glb) ) if c_glb > 0 else 0

        p_loc = min(ROOT.TMath.Erfc(z_loc / np.sqrt(2)), 1)
        p_glb = min(ROOT.TMath.Erfc(z_glb / np.sqrt(2)), 1)

        p_local .SetBinContent(bin, p_loc)
        p_global.SetBinContent(bin, p_glb)
        p_local .SetBinError  (bin, 0)
        p_global.SetBinError  (bin, 0)
        pass

    c2.plot(p_local,  markercolor=ROOT.kGreen  + 1, linecolor=ROOT.kGreen  + 1, option='PL', label="Local (20% window)")
    c2.plot(p_global, markercolor=ROOT.kViolet + 1, linecolor=ROOT.kViolet + 1, option='PL', label="Global (GBS)")
    c2.xlabel("Signal jet mass [GeV]")
    c2.ylabel("p_{0}")
    c2.log()

    c2.ylim(1E-04, 1E+04)
    for sigma in range(4):
        c2.yline(ROOT.TMath.Erfc(sigma / np.sqrt(2)))
        pass

    c2.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
            "Sherpa incl. #gamma MC",
            "Trimmed anti-k_{t}^{R=1.0} jets",
            "ISR #gamma selection",
            ("Signal" if args.inject else "No signal") + " injected" + (" at m = %d GeV" % args.mass if args.inject else ""),
            ], 
           qualifier='Simulation Internal')

    c2.legend()
    if args.save: c2.save('plots/globalbackground_p0_%dGeV_%s.pdf' % (args.mass, 'injected' if args.inject else 'notinjected'))
    if args.show: c2.show()
    
    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
