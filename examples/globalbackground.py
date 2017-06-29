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
parser.add_argument('--data', dest='data', action='store_const',
                    const=True, default=False,
                    help='Use data (default: False/MC)')
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

    DSID = int("100%03d" % args.mass)


    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Get signal file
    sig_DSID = get_signal_DSID(args.mass, tolerance=10)
    if sig_DSID is None:
        warning ("No signal file was found")
        return
    sig_file = 'objdef_MC_{DSID:6d}.root'.format(DSID=sig_DSID)

    # Load data
    files = {
        'data': glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root'),
        'gbs':  glob.glob(tf.config['base_path'] + 'objdef_GBSMC_400001.root'),
        'WZ':   glob.glob(tf.config['base_path'] + 'objdef_MC_3054*.root')
        }
    if args.inject:
        files['sig'] = glob.glob(tf.config['base_path'] + sig_file)
        pass

    if len(files) == 0:
        warning("No files found. Try to run:")
        warning(" $ source getSomeData.sh")
        return

    data   = loadData(files['data'], tf.config['tree'],      prefix=tf.config['prefix'])
    gbs    = loadData(files['gbs'],  tf.config['finaltree'], prefix=tf.config['prefix'])
    WZ     = loadData(files['WZ'],   tf.config['tree'],      prefix=tf.config['prefix'])
    if args.inject:
        signal = loadData(files['sig'], tf.config['tree'],   prefix=tf.config['prefix'])
    else:
        signal = None
        pass
    info   = {key: loadData(files[key], tf.config['outputtree'], stop=1) for key in files}

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])

    # Append new DSID field
    if args.inject:
        signal = append_fields(signal, 'DSID', np.zeros((signal.size,)), dtypes=int)
        for idx, id in enumerate(info['sig']['id']):
            msk = (signal['id'] == id) # Get mask of all 'signal' entries with same id, i.e. from same file
            DSID = info['sig']['DSID'][idx]  # Get DSID for this file
            signal['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
            signal['DSID']  [msk] = DSID # Store DSID
            pass
        signal['weight'] *= tf.config['lumi']
        pass

    WZ = append_fields(WZ, 'DSID', np.zeros((WZ.size,)), dtypes=int)
    for idx, id in enumerate(info['WZ']['id']):
        msk = (WZ['id'] == id) # Get mask of all 'WZ' entries with same id, i.e. from same file
        DSID = info['WZ']['DSID'][idx]  # Get DSID for this file
        WZ['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
        WZ['DSID']  [msk] = DSID # Store DSID
        pass
    WZ['weight'] *= tf.config['lumi']

    #if not args.data:
    data = append_fields(data, 'DSID', np.zeros((data.size,)), dtypes=int)
    for idx, id in enumerate(info['data']['id']):
        msk = (data['id'] == id) # Get mask of all 'data' entries with same id, i.e. from same file
        DSID = info['data']['DSID'][idx]  # Get DSID for this file
        data['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
        data['DSID']  [msk] = DSID # Store DSID
        pass
    data['weight'] *= tf.config['lumi']
    #pass

    # Compute new variables
    data       = append_fields(data,   'logpt', np.log(data  ['pt']))
    WZ         = append_fields(WZ,     'logpt', np.log(WZ    ['pt']))
    if signal is not None:
        signal = append_fields(signal, 'logpt', np.log(signal['pt']))
        pass

    # Inject signal into data
    if args.inject:
        data = np.array(np.concatenate((data,signal)), dtype=data.dtype)
        pass
    #if not args.data:
    data = np.array(np.concatenate((data,WZ)), dtype=data.dtype)
    #pass

    """ @TODO: Not sure this script works for data input... But it's not used anyway. """

    
    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Pass/fail masks
    # -- Data (incl. signal)
    msk_pass = tf.config['pass'](data)
    msk_fail = ~msk_pass

    # -- W/Z
    msk_WZ_pass = tf.config['pass'](WZ)
    msk_WZ_fail = ~msk_WZ_pass

    # -- Signal
    if args.inject:
        msk_sig_pass = tf.config['pass'](signal)
        msk_sig_fail = ~msk_sig_pass
        pass

    # Transfer factor calculator instance
    calc = tf.calculator(data=data, config=tf.config, subtract=WZ)
    
    # Nominal fit
    calc.fit()
    w_nom    = calc.weights(data[msk_fail])
    w_nom_WZ = calc.weights(WZ  [msk_WZ_fail])
    if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/globalbackground_%s_%s_' % ('injected' if args.inject else 'notinjected', 'data' if args.data else 'MC'))

    # mass +/- 20% stripe fit
    calc.mass   = args.mass
    calc.window = 0.2
    calc.fit()
    w_stripe         = calc.weights(data[msk_fail])
    w_stripe_WZ      = calc.weights(WZ  [msk_WZ_fail])
    if args.inject:
        w_stripe_sig = calc.weights(signal[msk_sig_fail])
        pass
   
    
    # Plotting
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    bins = np.linspace(100, 250, 30 + 1, endpoint=True) 
    
    # Setup canvas
    c = ap.canvas(num_pads=2, batch=not args.show)
    p0, p1 = c.pads()

    # Add stacked backgrounds
    h_bkg_nom     = c.hist(data  ['m'][msk_fail],     bins=bins, weights=data  ['weight'][msk_fail] * w_nom,            display=False)
    h_bkg_stripe  = c.hist(data  ['m'][msk_fail],     bins=bins, weights=data  ['weight'][msk_fail] * w_stripe,         display=False)
    h_WZfl_nom    = c.hist(WZ    ['m'][msk_WZ_fail],  bins=bins, weights=WZ    ['weight'][msk_WZ_fail] * w_nom_WZ,      display=False)
    h_WZfl_stripe = c.hist(WZ    ['m'][msk_WZ_fail],  bins=bins, weights=WZ    ['weight'][msk_WZ_fail] * w_stripe_WZ,   display=False)
    if args.inject:
        h_sig     = c.hist(signal['m'][msk_sig_pass], bins=bins, weights=signal['weight'][msk_sig_pass],                display=False)
        h_sfl     = c.hist(signal['m'][msk_sig_fail], bins=bins, weights=signal['weight'][msk_sig_fail] * w_stripe_sig, display=False)
        pass
    h_gbs         = c.hist(gbs   ['m'],               bins=bins, weights=gbs   ['weight'],                              display=False)

    # -- Subtract (opt.)
    if args.inject:
        h_bkg_stripe.Add(h_sfl, -1)
        h_gbs       .Add(h_sfl, -1)
        pass
    h_bkg_nom   .Add(h_WZfl_nom,    -1)
    h_bkg_stripe.Add(h_WZfl_stripe, -1)

    # -- Actually draw
    #if not args.data:
    h_WZ         = c.stack(WZ['m'][msk_WZ_pass], bins=bins, weights=WZ['weight'][msk_WZ_pass], fillcolor=ROOT.kRed - 4, label='W/Z + #gamma')
    #pass

    h_bkg_nom    = c.stack(h_bkg_nom,    fillcolor=ROOT.kAzure  + 7, label="Bkg. (full)")
    h_sum = c.getStackSum()
    h_bkg_stripe.Add(h_WZ)
    h_gbs       .Add(h_WZ)
    if args.inject:
        h_sig    = c.stack(h_sig,        fillcolor=ROOT.kViolet - 4,    label="Z' (%d GeV)" % args.mass)
        pass
    h_bkg_stripe = c.hist (h_bkg_stripe, linecolor=ROOT.kGreen  + 1, label="Bkg. (window)")# % args.mass)
    h_gbs        = c.hist (h_gbs,        linecolor=ROOT.kViolet + 1, label="Bkg. (GBS)")
    


    # Draw stats. error of stacked sum
    h_sum = c.hist(h_sum, fillstyle=3245, fillcolor=ROOT.kGray+2, linecolor=ROOT.kGray + 3, label='Stats. uncert.', option='E2')
    
    # Add (pseudo-) data
    h_data = c.plot (data['m'][msk_pass], bins=bins, weights=data['weight'][msk_pass], markersize=0.8, label='Data' if args.data else 'Pseudo-data')

    # Axis limits
    p1.ylim(0.8, 1.2)
    c.padding(0.45)
    c.log(True)

    # Draw error- and ratio plots
    if args.inject:
        hr_sig  = c.ratio_plot((h_sig,       h_sum), option='HIST', offset=1)
        pass
    h_err   = c.ratio_plot((h_sum,       h_sum), option='E2')
    h_ratio = c.ratio_plot((h_data,      h_sum), oob=True)
    h_rgbs  = c.ratio_plot((h_gbs,       h_sum), linecolor=ROOT.kViolet + 1, option='HIST ][')
    h_rgbs  = c.ratio_plot((h_bkg_stripe,h_sum), linecolor=ROOT.kGreen  + 1, option='HIST ][')
    
    # Add labels and text
    c.xlabel('Signal jet mass [GeV]')
    c.ylabel('Events')
    p1.ylabel('Data / Nom.')
    c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",] +
           (["Sherpa incl. #gamma MC",] if not args.data else [])+
           ["Trimmed anti-k_{t}^{R=1.0} jets",
            "ISR #gamma selection",] +
           (["Signal injected"] if args.inject else []),
           qualifier='%sInternal' % ("Simulation " if not args.data else ""))

    # Add line(s)
    p1.yline(1.0)

    # Draw legend
    c.legend()
    c.region("SR", 0.8 * args.mass, 1.2 * args.mass)

    # Save and show plot
    if args.save: c.save('plots/globalbackground_spectrum_%dGeV_%s_%s.pdf' % (args.mass, 'injected' if args.inject else 'notinjected', 'data' if args.data else 'MC'))
    if args.show: c.show()


    # p0-plot
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
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

    c2.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",] + 
            (["Sherpa incl. #gamma MC",] if not args.data else []) + 
            ["Trimmed anti-k_{t}^{R=1.0} jets",
            "ISR #gamma selection",
            ("Signal" if args.inject else "No signal") + " injected" + (" at m = %d GeV" % args.mass if args.inject else ""),
            ], 
           qualifier='Simulation Internal')

    c2.region("SR", 0.8 * args.mass, 1.2 * args.mass)
    c2.legend()
    if args.save: c2.save('plots/globalbackground_p0_%dGeV_%s_%s.pdf' % (args.mass, 'injected' if args.inject else 'notinjected', 'data' if args.data else 'MC'))
    if args.show: c2.show()
    
    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
