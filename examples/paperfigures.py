#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for producing harmonised plots for the di-jet + ISR paper

@file:   paperfigures.py
@author: Andreas Søgaard
@date:   13 June 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import sys, glob, itertools

# Get ROOT to stop hogging the command-line options
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

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
    

# Main function.
def main ():

    path = "/afs/cern.ch/user/a/asogaard/public/dijetisr/paperfigures/"


    colours = {
        'QCD':  ROOT.kAzure  + 7,
        'WHad': ROOT.kRed    - 4,
        'ZHad': ROOT.kOrange + 4,
        'syst': ROOT.kGreen  + 1,
        }
            

    # W/Z plots
    # --------------------------------------------------------------------------
    if False:
        for isrjet in [True, False]:
            if isrjet:
                f = ROOT.TFile(path + "hists_isrjet_WZ.root")
            else:
                f = ROOT.TFile(path + "hists_isrgamma_WZ.root")
                pass
            
            names = ['QCD', 'QCD_up', 'QCD_down', 'data', 'WHad', 'ZHad']
            
            histograms = dict()
            for name in names:
                hist = f.Get('h_mJ_' + name)
                hist.SetDirectory(0)
                histograms[name] = hist
                pass
            
            f.Close()
            
            # Combine  W/Z component
            histograms['WZHad'] = histograms['WHad'].Clone(histograms['WHad'].GetName().replace('W', 'WZ'))
            histograms['WZHad'].Add(histograms['ZHad'])
            
            # -- canvas
            c = ap.canvas(num_pads=2, fraction=0.45)
            pads = c.pads()
            
            # -- main pad
            c.stack(histograms['QCD'],  fillcolor=colours['QCD'],  label='Background est.')
            c.stack(histograms['WHad'], fillcolor=colours['WHad'], label='W + %s' % ('jets' if isrjet else '#gamma'))
            c.stack(histograms['ZHad'], fillcolor=colours['ZHad'], label='Z + %s' % ('jets' if isrjet else '#gamma'))
            h_sum = c.hist(histograms['QCD'], fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                           label='Stat. uncert.')
            
            h_bkg_up   = c.hist(histograms['QCD_up'],   linestyle=2, linecolor=colours['syst'],
                                label='Syst. uncert.')
            h_bkg_down = c.hist(histograms['QCD_down'], linestyle=2, linecolor=colours['syst'])
            
            c.plot (histograms['data'], label='Data')
            
            # -- diff pad
            pads[1].hist(histograms['WZHad'], fillcolor=colours['WHad'], option='HIST')
            pads[1].hist(histograms['ZHad'],  fillcolor=colours['ZHad'], option='HIST')
            c.diff_plot((h_bkg_up,   histograms['QCD']), option='HIST')
            c.diff_plot((h_bkg_down, histograms['QCD']), option='HIST')
            c.diff_plot((histograms['data'], histograms['QCD']))
            
            # -- decorations
            c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
                    "%s channel" % ('Jet' if isrjet else 'Photon')],
                   qualifier="Internal")
            c.xlabel('Large-#it{R} jet mass [GeV]')
            c.ylabel('Events / 2 GeV')
            c.padding(0.37)
            if isrjet: pads[1].ylim(-1000, 5000)
            else:      pads[1].ylim( -140,  700)
            pads[1].ylabel('Data - background est.')
            pads[1].yline(0)
            c.region("SR", 0.8 * 8.5, 1.2 * 85.)
            c.legend(ymax=0.875, reverse=True)
            c.save('plots/paperplot_wz_%s.pdf' % ('jet' if isrjet else 'gamma'))
            c.show()
            pass
        pass


    # Search plots
    # --------------------------------------------------------------------------
    if True:
        for isrjet, mass, log in itertools.product([True, False], [150, 220], [True, False]):
            if isrjet:
                f = ROOT.TFile(path + "hists_isrjet_%d.root" % mass)
            else:
                f = ROOT.TFile(path + "hists_isrgamma_datadistributions_%dGeV.root" % mass)
                pass
            
            names = ['QCD', 'QCD_up', 'QCD_down', 'data', 'WHad', 'ZHad']#, 'sig']
            
            histograms = dict()
            for name in names:
                hist = f.Get('h_mJ_' + name)
                hist.SetDirectory(0)
                histograms[name] = hist
                pass
            
            f.Close()
            
            # Combine W/Z components
            histograms['WZHad'] = histograms['WHad'].Clone(histograms['WHad'].GetName().replace('W', 'WZ'))
            histograms['WZHad'].Add(histograms['ZHad'])
      
            # Add stat. and syst. uncert. in quadrature
            for bin in range(1, histograms['QCD'].GetXaxis().GetNbins() + 1):
                stat = histograms['QCD']     .GetBinError(bin)
                nom  = histograms['QCD']     .GetBinContent(bin)
                up   = histograms['QCD_up']  .GetBinContent(bin)
                down = histograms['QCD_down'].GetBinContent(bin)
                syst = max(abs(up - nom), abs(down - nom))
                histograms['QCD_up']  .SetBinContent(bin, nom + np.sqrt(np.square(stat) + np.square(syst)))
                histograms['QCD_down'].SetBinContent(bin, nom - np.sqrt(np.square(stat) + np.square(syst)))
                pass
            # ...

            # Add W/Z component to background variations
            histograms['QCD_up']  .Add(histograms['WZHad'])
            histograms['QCD_down'].Add(histograms['WZHad'])
      
            
            # -- canvas
            c = ap.canvas(batch=True, num_pads=2)
            pads = c.pads()
            
            # -- main pad
            if isrjet and False:
                c.stack(histograms['QCD'],   fillcolor=colours['QCD'],  label='Background est.')
                c.stack(histograms['WZHad'], fillcolor=colours['WHad'], label='W/Z + %s' % ('jets' if isrjet else '#gamma'))
            else:
                c.stack(histograms['WZHad'], fillcolor=colours['WHad'], label='W/Z + %s' % ('jets' if isrjet else '#gamma'))
                c.stack(histograms['QCD'],   fillcolor=colours['QCD'],  label='Background est.')
                pass

            h_sum = c.getStackSum()
            h_sum = c.hist(h_sum, fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                           label='Stat. uncert.')
            
            h_bkg_up   = c.hist(histograms['QCD_up'],   linestyle=2, linecolor=colours['syst'],
                                label='Stat. #oplus syst.')
            h_bkg_down = c.hist(histograms['QCD_down'], linestyle=2, linecolor=colours['syst'])
            
            #c.hist(histograms['sig'], label='Data')

            c.plot (histograms['data'], label='Data')
            
            # -- diff pad
            if isrjet: pads[1].ylim(0.95, 1.05)
            else:      pads[1].ylim(0.8, 1.2)

            c.ratio_plot((h_bkg_up,   h_sum), option='HIST')
            c.ratio_plot((h_bkg_down, h_sum), option='HIST')
            c.ratio_plot((h_sum,      h_sum), option='E2')
            c.ratio_plot((histograms['data'], h_sum), oob=True)
            
            # -- decorations
            c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
                    "%s channel" % ('Jet' if isrjet else 'Photon')],
                   qualifier="Internal")
            c.xlabel('Large-#it{R} jet mass [GeV]')
            c.ylabel('Events / 5 GeV')
            c.padding(0.35)
            pads[1].ylabel('Data / est.')
            pads[1].yline(1)
            c.region("SR", 0.8 * mass, 1.2 * mass, offset=(0.14 if (mass == 150 and not isrjet) else 0.07))
            c.legend(reverse=True) # ymax=0.87
            c.log(log)
            c.save('plots/paperplot_%s_%dGeV%s.pdf' % ('jet' if isrjet else 'gamma', mass, '_log' if log else ''))
            c.show()
            pass
        pass

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
