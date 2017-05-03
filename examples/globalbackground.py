#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for computing the global background shape (GBS) for the dijet + ISR analysis.

@file:   globalbackground.py
@date:   1 May 2017
@author: Andreas Søgaard 
@email:  andreas.sogaard@cern.ch
"""

# Basic
import sys, glob

# Scientific import(s)
import ROOT
import numpy as np
from root_numpy import *

# Local import(s)
from rootplotting.tools import *
from rootplotting import ap

import transferfactor as tf

# Main function
def main ():
    

    # Setup.
    # – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – 

    # Load data
    files  = glob.glob('data/objdef_MC_3610*.root') + ['data/objdef_MC_308365.root'] #  + glob.glob('data/objdef_MC_3054*.root') 

    if len(files) == 0:
        warning("No files found. Try to run:")
        warning(" $ source getSomeData.sh")
        return

    tree   = 'BoostedJet+ISRgamma/Nominal/EventSelection/Pass/NumLargeRadiusJets/Postcut'
    prefix = 'Jet_'

    data = loadData(files, tree)
    info = loadData(files, 'BoostedJet+ISRgamma/Nominal/outputTree', stop=1)

    # Rename variables
    data.dtype.names = [name.replace(prefix, '') for name in data.dtype.names]

    # Scaling by cross section
    xsec = loadXsec('sampleInfo.csv')
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
    
    msk_sig  = (data['DSID'] == 308365)
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
    calc = tf.calculator(data=data, config=tf.config) # incl
    
    # Nominal fit
    calc.fit()
    w_nom  = calc.weights(data[msk_fail])
    #calc.plot()

    # 160 GeV +/- 20% fit
    calc.mass = 160
    calc.window = 0.2
    calc.fit()
    w_160GeV = calc.weights(data[msk_fail])
    
    
    # Get GBS
    bins = np.linspace(50, 300, 50 + 1, endpoint=True)

    masses = [120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220]
    dict_backgrounds = dict()
    ctemp = ap.canvas(batch=True)
    for mass in masses:
        print " --", mass
        # Fit TF profile
        calc.mass = mass
        calc.fit()

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

    print backgrounds
    print backgrounds.mean(axis=0), len(backgrounds.mean(axis=0))
    print backgrounds.mean(axis=1), len(backgrounds.mean(axis=1))

    
    # Plotting
    # – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – – 
    #bins = np.linspace(50, 300, 25 + 1, endpoint=True)
    
    # Setup canvas
    c = ap.canvas(num_pads=2)
    p0, p1 = c.pads()

    # Add stacked backgrounds
    h_bkg_nom    = c.stack(data  ['m'][msk_fail],     bins=bins, weights=data  ['weight'][msk_fail] * w_nom,    fillcolor=ROOT.kAzure  + 7, label="Bkg. (Nom.)")
    h_sig        = c.stack(signal['m'][msk_sig_pass], bins=bins, weights=signal['weight'][msk_sig_pass],        fillcolor=ROOT.kRed    + 1, label="Z' (160 GeV)")
    h_bkg_160GeV = c.hist (data  ['m'][msk_fail],     bins=bins, weights=data  ['weight'][msk_fail] * w_160GeV, linecolor=ROOT.kGreen  + 1, label="Bkg. (160 GeV #pm 20%)")
    h_gbs        = c.hist(backgrounds.mean(axis=0),   bins=bins,                                                linecolor=ROOT.kViolet + 1, label="Bkg. (GBS)")
    
    # Draw stats. error of stacked sum
    h_sum  = h_bkg_nom.Clone('h_sum') #c.getStackSum()
    c.hist(h_sum, fillstyle=3245, fillcolor=ROOT.kGray+2, linecolor=ROOT.kGray + 3, label='Stats. uncert.', option='E2')
    
    # Add (pseudo-) data
    np.random.seed(21)
    h_data = c.plot (data['m'][msk_pass], bins=bins, weights=data['weight'][msk_pass], markersize=0.8, label='Pseudo-data', scale=1)

    # Draw error- and ratio plots
    h_err   = c.ratio_plot((h_sum,       h_bkg_nom), option='E2')
    h_ratio = c.ratio_plot((h_data,      h_bkg_nom), markersize=0.8)
    h_rgbs  = c.ratio_plot((h_gbs,       h_bkg_nom), linecolor=ROOT.kViolet + 1, option='HIST ][')
    h_rgbs  = c.ratio_plot((h_bkg_160GeV,h_bkg_nom), linecolor=ROOT.kGreen  + 1, option='HIST ][')
    
    # Add labels and text
    c.xlabel('Signal jet mass [GeV]')
    c.ylabel('Events')
    p1.ylabel('Data / Nom.')
    c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
            "Sherpa MC 15",
            "Trimmed anti-k_{t}^{R=1.0} jets",
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
    c.save('gbs.pdf')
    c.show()
    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
