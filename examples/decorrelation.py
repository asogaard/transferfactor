#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Script for performing DDT transform and producing publication-level plots.

@file:   decorrelation.py
@author: Andreas Søgaard
@date:   13 April 2017
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
    from snippets.functions import loadDataFast, loadXsec, displayName, displayNameUnit
    from rootplotting import ap
except ImportError:
    print "WARNING: This script uses the 'snippets' package. Clone it as e.g.:"
    print "         $ git clone git@github.com:asogaard/snippets.git"
    sys.exit()
    pass

# Command-line arguments parser
import argparse

parser = argparse.ArgumentParser(description='Perform decorrelation study.')

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

    paths = {
        'incl': glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root'),
        'W':    glob.glob(tf.config['base_path'] + 'objdef_MC_30543*.root'),
        'Z':    glob.glob(tf.config['base_path'] + 'objdef_MC_30544*.root'),
        'sig':  glob.glob(tf.config['base_path'] + 'objdef_MC_308365.root'),
        }

    # Load data
    xsec = loadXsec('../AnalysisTools/share/sampleInfo.csv')

    treenames = ['BoostedJet+ISRgamma/Nominal/LargeRadiusJets/Nominal/dPhiPhoton/Postcut',
                 'BoostedJet+ISRgamma/Nominal/LargeRadiusJets/Nominal/BoostedRegime/Postcut',
                 'BoostedJet+ISRgamma/Nominal/LargeRadiusJets/Nominal/rhoDDT/Postcut',
                 ]

    # Standard definitions
    # -- Data
    prefix = ''
    getvars = ['m', 'pt', 'pt_ungroomed', 'tau21', 'tau21_ungroomed']# , 'tau21DDT']

    xvars = ['rho', 'rhoDDT', 'm']
    yvars = ['tau21', 'tau21DDT', 'tau21_ungroomed']
    axis = {
        'rho':    (30, -9, 0),
        'rhoDDT': (30, -2, 7),
        'm':      (30, 0, 300),
        }
        
    # -- Plotting
    colours = [ROOT.kRed + ((i * 2) % 5) for i in range(3)]

    qualifier = 'Simulation Internal'

    lines = ["#sqrt{s} = 13 TeV,  Incl. #gamma MC",
             "Trimmed anti-k_{t}^{R=1.0} jets",
             ]

    xmin, xmax, ymax = 0.57, 0.9, 0.85 - ROOT.gROOT.GetStyle("AStyle").GetTextSize() * 1.45 # 1.4


    # Fitting transform
    # ----------------------------------------------------------------
    data = {key: loadDataFast(paths[key], treenames[1], ['m', 'pt', 'tau21'], prefix, xsec, quiet=True) for key in paths}
        
    # Compute new variables
    for key in data:
        data[key]['rho']    = np.log(np.square(data[key]['m']) / np.square(data[key]['pt']))
        data[key]['rhoDDT'] = np.log(np.square(data[key]['m']) / (data[key]['pt'] * 1.))
        pass

    # Compute profile
    xvar, yvar = 'rhoDDT', 'tau21'
    profile_before = get_profile(data['incl'], xvar, yvar, 
                                 bins=(28, 0, axis[xvar][2]), 
                                 weights=data['incl']['weight'])

    fit = ROOT.TF1('fit', 'pol1', 1.5, 7)
    fit.SetLineColor(ROOT.kGray + 3)
    fit.SetLineWidth(2)
    profile_before.Fit('fit', 'RQ0')

    print "Fit:"
    print "  ", fit.GetParameter(0), "+/-", fit.GetParError(0)
    print " ",  fit.GetParameter(1), "+/-", fit.GetParError(1)
    
    for key in data:
        data[key]['tau21DDT'] = data[key]['tau21'] - (fit.GetParameter(1) * (data[key]['rhoDDT'] - 1.5))
        pass

    xvar, yvar = 'rhoDDT', 'tau21DDT'
    profile_after = get_profile(data['incl'], xvar, yvar,
                                 bins=(28, 0, axis[xvar][2]),
                                 weights=data['incl']['weight'])

    # Draw
    c = ap.canvas(batch=not args.show)

    c.plot(profile_before, linecolor=colours[0], markercolor=colours[0], label='Original, #tau_{21}',          option='PE', legend_option='LP')
    c.plot(profile_after,  linecolor=colours[1], markercolor=colours[1], label='Transformed, #tau_{21}^{DDT}', option='PE', legend_option='LP')
    fit.Draw('SAME')

    c.legend()
    c.xlabel('Large-#it{R} jet %s' % displayNameUnit(xvar))
    c.ylabel("#LT%s#GT" % displayName('tau21') + ", #LT%s#GT" % displayName(yvar))
    c.ylim(0, 1.2)
    
    c.text(lines + ["Jet p_{T} > 2 #times M, 200 GeV"],
           qualifier=qualifier)
    c.latex("Fit: %.3f + %s #times (%.4f)" % (fit.GetParameter(0), displayName(xvar), fit.GetParameter(1)), 0.20, 0.20, NDC=True, align=11)
    
    c.line(1.5, fit.Eval(1.5), axis[xvar][2], fit.Eval(1.5), linecolor=ROOT.kGray + 3, linestyle=2, linewidth=2)
    c.line(1.5, 0.14, 1.5, fit.Eval(1.5),                    linecolor=ROOT.kGray,     linestyle=1, linewidth=1)
    c.line(1.5, 0.0,  1.5, 0.05,                             linecolor=ROOT.kGray,     linestyle=1, linewidth=1)
    
    savename = 'plots/decorrelation_fit_%s_vs_%s.pdf' % ('tau21', xvar)
    if args.save: c.save(savename)
    if args.show: c.show()


    # Performing pT-slicing
    # --------------------------------------------------------------------------
    # Loop trees (succesive cuts)
    for itreename, treename in enumerate(treenames):
        
        data = {key: loadDataFast(paths[key], treename, getvars, prefix, xsec, quiet=True) for key in paths}
        # Compute new variables
        for key in data:
            data[key]['rho']    = np.log(np.square(data[key]['m']) / np.square(data[key]['pt']))
            data[key]['rhoDDT'] = np.log(np.square(data[key]['m']) / (data[key]['pt'] * 1.))
            data[key]['tau21DDT'] = data[key]['tau21'] - (fit.GetParameter(1) * (data[key]['rhoDDT'] - 1.5))
            pass
        
        # Profile settings
        slices = {
            'pt': [
                ( 200,  300),
                ( 300,  500),
                ( 500, 1000),
                ],
            }
        
        # Compute profiles
        profiles = {key: {xvar: {yvar: [] for yvar in yvars} for xvar in xvars} for key in data.keys()}
        
        for key, xvar, yvar in itertools.product(data.keys(), xvars, yvars):
            if key == 'sig': # Restrict signal to +/- 20% of pole mass
                msk_sig = np.abs(data[key]['m'] - 160.) / 160. < 0.2
            else:
                msk_sig = np.ones_like(data[key][xvar]).astype(bool)
                pass
            for sl in slices['pt']:
                profile = get_profile(data[key], xvar, yvar, 
                                      bins=axis[xvar], 
                                      mask=(data[key]['pt'] > sl[0]) & (data[key]['pt'] < sl[1]) & msk_sig, 
                                      weights=data[key]['weight'])
                profiles[key][xvar][yvar].append(profile)
                pass
            pass
     
   
        # Comparing groomed and un-groomed
        # ----------------------------------------------------------------------

        # Drawing options
        names = ['[%d, %d] GeV' % (sl[0], sl[1]) for sl in slices['pt']]
        
        # Loop xvars
        for xvar in xvars:
            yvar1, yvar2 = 'tau21', 'tau21_ungroomed'
                        
            c = ap.canvas(num_pads=2, batch=not args.show)

            # -- main pad
            for i, (col, name, prof) in enumerate(zip(colours, names, profiles['incl'][xvar][yvar1])):
                c.plot(prof, linecolor=col, markercolor=col, markerstyle=20, option='PE', label=name, legend_option='LP')
                pass
            for i, (col, name, prof) in enumerate(zip(colours, names, profiles['incl'][xvar][yvar2])):
                c.plot(prof, linecolor=col, markercolor=col, markerstyle=24, option='PE')
                pass

            # -- ratio pad
            for prof1, prof2, col in zip(profiles['incl'][xvar][yvar1], profiles['incl'][xvar][yvar2], colours):
                r = c.ratio_plot((prof1, prof2), option='HIST ][', linecolor=col, default=0)
                pass
            
            # -- decorations
            c.xlabel('Large-#it{R} jet %s' % displayNameUnit(xvar))
            c.ylabel("#LT%s#GT" % displayName(yvar1) + ", #LT%s#GT" % displayName(yvar2))
            c.pads()[1].ylabel("#LT%s#GT" % displayName(yvar1) + " / #LT%s#GT" % displayName(yvar2))

            c.ylim(0, 1.4)# if xvar == 'm' else 1.6)
            c.pads()[1].ylim(0, 1.5)

            c.legend(header='Jet p_{T} in:', categories=[
                    ('Trimmed %s'    % displayName('tau21'), {'markerstyle': 20, 'markercolor': ROOT.kGray + 3, 'option': 'P'}),
                    ('Un-trimmed %s' % displayName('tau21'), {'markerstyle': 24, 'markercolor': ROOT.kGray + 3, 'option': 'P'}),
                    ])
            c.text(lines + (["Jet p_{T} > 2 #times M"] if itreename > 0 else []) + (["Jet %s > 1.5" % displayName('rhoDDT')] if itreename > 1 else []),
                   qualifier=qualifier)

            savename = 'plots/decorrelation_trimmed_untrimmed_%s_vs_%s_tree%d.pdf' % (yvar1, xvar, itreename)
            if args.save: c.save(savename)
            if args.show: c.show()
            pass


        # Only groomed
        # ----------------------------------------------------------------------

        # Loop xvars
        for xvar, yvar in itertools.product(xvars, yvars):
               
            c = ap.canvas(batch=not args.show)

            for i, (col, name, prof) in enumerate(zip(colours, names, profiles['incl'][xvar][yvar])):
                c.plot(prof, linecolor=col, markercolor=col, markerstyle=20, label=name, option='PE', legend_option='LP')
                pass

            c.xlabel('Large-#it{R} jet %s' % displayNameUnit(xvar))
            c.ylabel("#LT%s#GT" % displayName(yvar))
            c.ylim(0, 1.4)# if xvar == 'm' else 1.6)

            c.legend(header='Jet p_{T} in:')
            c.text(lines + (["Jet p_{T} > 2 #times M"] if itreename > 0 else []) + (["Jet %s > 1.5" % displayName('rhoDDT')] if itreename > 1 else []),
                   qualifier=qualifier)

            savename = 'plots/decorrelation_%s_vs_%s_tree%d.pdf' % (yvar, xvar, itreename)
            if args.save: c.save(savename)
            if args.show: c.show()
            pass # end: loop xvars, yvars


        # Comparing QCD and W/Z
        # ----------------------------------------------------------------------

        # Drawing options
        names = ['[%d, %d] GeV' % (sl[0], sl[1]) for sl in slices['pt']]

        # Loop xvars
        for xvar, yvar in itertools.product(xvars, yvars):

            # Save outputs for harmonised plotting
            if itreename == 2 and yvar.startswith('tau21') and xvar in  ['m', 'rhoDDT']:
                outfile = ROOT.TFile('hists_isrgamma_decorrelation_{}_vs_{}.root'.format(yvar, xvar), 'RECREATE')
                profname = 'profile_{{name:s}}__{yvar:s}_vs_{xvar:s}__pT_{{ptmin:d}}_{{ptmax:d}}_GeV'
                profname = profname.format(yvar=yvar, xvar=xvar)
                print ">" * 80
                for prof, sl in zip(profiles['incl'][xvar][yvar], slices['pt']):
                    prof.SetName(profname.format(name='inclphoton', ptmin=sl[0], ptmax=sl[1]))
                    print ">>>", prof.GetName()
                    prof.Write()
                    pass
                for prof, sl in zip(profiles['sig'][xvar][yvar], slices['pt']):
                    prof.SetName(profname.format(name='sig160', ptmin=sl[0], ptmax=sl[1]))
                    print ">>>", prof.GetName()
                    prof.Write()
                    pass
                print ">" * 80
                outfile.Close()
                pass


            c = ap.canvas(batch=not args.show)

            # -- main pad
            for i, (col, name, prof) in enumerate(zip(colours, names, profiles['incl'][xvar][yvar])):
                c.plot(prof, linecolor=col, markercolor=col, markerstyle=20, option='PE', label=name, legend_option='LP')
                pass
            for i, (col, name, prof) in enumerate(zip(colours, names, profiles['sig'][xvar][yvar])):
                c.plot(prof, linecolor=col, markercolor=col, markerstyle=24, option='PE')
                pass

            # -- decorations
            c.xlabel('Large-#it{R} jet %s' % displayNameUnit(xvar))
            c.ylabel("#LT%s#GT" % displayName(yvar))

            c.ylim(0, 1.4)# if xvar == 'm' else 1.6)

            c.legend(header='Jet p_{T} in:', categories=[
                    ("Incl. #gamma MC",      {'markerstyle': 20, 'markercolor': ROOT.kGray + 3, 'option': 'P'}),
                    ("Z'(160 GeV) + #gamma", {'markerstyle': 24, 'markercolor': ROOT.kGray + 3, 'option': 'P'}),
                    ])
            c.text(lines[:1] + [lines[1].replace(', Incl. #gamma', '')] + (["Jet p_{T} > 2 #times M"] if itreename > 0 else []) + (["Jet %s > 1.5" % displayName('rhoDDT')] if itreename > 1 else []),
                   qualifier=qualifier)

            savename = 'plots/decorrelation_inclphoton_sig_%s_vs_%s_tree%d.pdf' % (yvar, xvar, itreename)
            if args.save: c.save(savename)
            if args.show: c.show()
            pass

        pass # end: loop treenames
    
    return


# Utility function(s)
def get_profile (data, xvar, yvar, bins, mask=None, weights=None):
    """ ... """

    # Check(s)
    if mask is None:
        mask    = np.ones_like(data[xvar]).astype(bool)
        pass

    if weights is None:
        weights = np.ones_like(data[xvar])
        pass

    # Initialise output
    profile = ROOT.TProfile('prof_%s_vs_%s_slice%d' % (xvar, yvar, np.sum(mask)), "", *bins)
    
    # Profile fill matrix
    M = np.vstack((data[xvar], data[yvar])).T

    # Weight, masked fill profile
    fill_profile(profile, M[mask,:], weights=weights[mask])

    return profile


# Main function call.
if __name__ == '__main__':
   main()
   pass

