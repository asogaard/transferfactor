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
    from rootplotting.style import *
    from snippets.functions import dict_product
except ImportError:
    print "WARNING: This script uses the 'transferfactor' and 'rootplotting' packages. Clone them as e.g.:"
    print "         $ git clone git@github.com:asogaard/transferfactor.git"
    print "         $ git clone git@github.com:asogaard/rootplotting.git"
    sys.exit()
    pass
    
# Command-line arguments parser
import argparse

parser = argparse.ArgumentParser(description='Produce final plots for paper.')

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


    path = "/afs/cern.ch/user/a/asogaard/public/dijetisr/paperfigures/"

    compcolours = {
        'QCD':  ROOT.kAzure  + 7,
        'WHad': ROOT.kRed    - 4,
        'ZHad': ROOT.kOrange + 4,
        'syst': ROOT.kGreen  + 1,
        }

    #qualifier = "Work In Progress"
    qualifier = "Internal"
   

    # W/Z plots
    # --------------------------------------------------------------------------
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
        c = ap.canvas(num_pads=2, fraction=0.45, batch=not args.show)
        pads = c.pads()
        
        # -- main pad
        h_mc   = c.stack(histograms['QCD'],  fillcolor=compcolours['QCD'],  label='Background est.')
        h_zhad = c.stack(histograms['ZHad'], fillcolor=compcolours['ZHad'], label='Z + %s' % ('jets' if isrjet else '#gamma'))
        h_whad = c.stack(histograms['WHad'], fillcolor=compcolours['WHad'], label='W + %s' % ('jets' if isrjet else '#gamma'))
        h_sum = c.hist(histograms['QCD'], fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                       label='Stat. uncert.')
        
        h_bkg_up   = c.hist(histograms['QCD_up'],   linestyle=2, linecolor=compcolours['syst'], label='Syst. uncert.')
        h_bkg_down = c.hist(histograms['QCD_down'], linestyle=2, linecolor=compcolours['syst'])
        
        h_data = c.plot (histograms['data'], label='Data')
        
        # -- diff pad
        pads[1].hist(histograms['WZHad'], fillcolor=compcolours['WHad'], option='HIST')
        pads[1].hist(histograms['ZHad'],  fillcolor=compcolours['ZHad'], option='HIST')
        c.diff_plot((h_bkg_up,   histograms['QCD']), option='HIST')
        c.diff_plot((h_bkg_down, histograms['QCD']), option='HIST')
        c.diff_plot((histograms['data'], histograms['QCD']))
        
        # -- statistical test
        h_mc.Add(h_zhad)
        h_mc.Add(h_whad)
        for bin in range(1, h_mc.GetXaxis().GetNbins()):
            stat = h_mc.GetBinError(bin)

            nom  = h_sum     .GetBinContent(bin)
            up   = h_bkg_up  .GetBinContent(bin)
            down = h_bkg_down.GetBinContent(bin)
            syst = max(abs(up-nom), abs(down-nom))

            h_mc.SetBinError(bin, np.sqrt( np.square(stat) + np.square(syst) ))
            pass
        chi2 = h_data.Chi2Test      (h_mc, "UW CHI2 P")
        KS   = h_data.KolmogorovTest(h_mc)
        ndf = h_mc.GetXaxis().GetNbins() - 1

        # -- decorations
        c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
                "%s channel" % ('Jet' if isrjet else 'Photon')],
               qualifier=qualifier)
        c.xlabel('Large-#it{R} jet mass [GeV]')
        c.ylabel('Events / 2 GeV')
        c.padding(0.47) # (0.45)
        if isrjet: pads[1].ylim(-1000, 5000)
        else:      pads[1].ylim( -140,  700)
        pads[1].ylabel('Data - background est.')
        pads[1].yline(0, linestyle=1, linecolor=ROOT.kBlack)
        c.region("SR", 0.8 * 85., 1.2 * 85.)
        c.legend(sort=True) # ymax=0.875,
        
        stats_string = "KS prob. = %.3f  |  #chi^{2}/N_{dof} (prob.) = %.1f/%d (%.3f)" % (KS, chi2, ndf, ROOT.TMath.Prob(chi2, ndf))
        #c.latex(stats_string, 0.95 , 0.96, NDC=True, align=31, textsize=16)

        if args.save: c.save('plots/paperplot_wz_%s.pdf' % ('jet' if isrjet else 'gamma'))
        if args.show: c.show()
        pass


    # Search plots
    # --------------------------------------------------------------------------
    for isrjet, mass, log in itertools.product([True, False], [160, 220], [True]):
        
        if isrjet:
            f = ROOT.TFile(path + "hists_isrjet_%d.root" % mass)
        else:
            try:
                # Try to get file with signal...
                f = ROOT.TFile(path + "hists_isrgamma_datadistributions_%dGeV_mu100.root" % mass)
            except:
                # ... otherwise, default to file without
                f = ROOT.TFile(path + "hists_isrgamma_datadistributions_%dGeV.root" % mass)
                pass
            pass
        
        names = ['QCD', 'QCD_up', 'QCD_down', 'data', 'WHad', 'ZHad', 'sig']
        
        histograms = dict()
        for name in names:
            hist = f.Get('h_mJ_' + name)
            try:
                hist.SetDirectory(0)
                histograms[name] = hist
            except: 
                # This mass point doesn't have a signal shape
                warning("Mass point %d doesn't have a signal shape." % mass)
                pass
            pass
        
        f.Close()

        # @TEMP: Blind search region for Laser's plots
        '''
        for bin in np.arange(histograms['data'].GetXaxis().GetNbins(), dtype=int) + 1:
            if histograms['data'].GetXaxis().GetBinUpEdge(bin) > 100.:
                histograms['data'].SetBinContent(bin, 0)
                histograms['data'].SetBinError  (bin, 0)
                pass
            pass
            '''
        
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

        
        # Add W/Z component to background variations
        histograms['QCD_up']  .Add(histograms['WZHad'])
        histograms['QCD_down'].Add(histograms['WZHad'])
        
        # -- canvas
        c = ap.canvas(num_pads=2, batch=not args.show)
        pads = c.pads()
        
        # -- main pad
        h_mc = c.stack(histograms['WZHad'], fillcolor=compcolours['WHad'], label='W/Z + %s' % ('jets' if isrjet else '#gamma'))
        h_qcd = c.stack(histograms['QCD'],   fillcolor=compcolours['QCD'],  label='Background est.')
          
        h_sum = c.getStackSum()
        h_sum = c.hist(h_sum, fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                       label='Stat. uncert.')
        
        h_bkg_up   = c.hist(histograms['QCD_up'],   linestyle=2, linecolor=compcolours['syst'],
                            label='Stat. #oplus syst.')
        h_bkg_down = c.hist(histograms['QCD_down'], linestyle=2, linecolor=compcolours['syst'])
        
        if 'sig' in histograms:
            c.hist(histograms['sig'], linestyle=1, linecolor=ROOT.kViolet+2, label="Z' (%d GeV)" % mass)
            pass
        
        h_data = c.plot (histograms['data'], label='Data')
        
        # -- diff pad
        if isrjet: pads[1].ylim(0.95, 1.05)
        else:      pads[1].ylim(0.8, 1.2)
        
        c.ratio_plot((h_bkg_up,   h_sum), option='HIST')
        c.ratio_plot((h_bkg_down, h_sum), option='HIST')
        c.ratio_plot((h_sum,      h_sum), option='E2')
        c.ratio_plot((histograms['data'], h_sum), oob=True) # @TEMP:, oob=True) taken out for Laser's plots
        
        # -- statistical test
        h_mc.Add(h_qcd)
        for bin in range(1, h_mc.GetXaxis().GetNbins()):
            nom  = h_mc      .GetBinContent(bin)
            up   = h_bkg_up  .GetBinContent(bin)
            down = h_bkg_down.GetBinContent(bin)
            statPlusSyst = max(abs(up-nom), abs(down-nom))

            # h_bkg_{up,down} are stats. oplus syst.
            h_mc.SetBinError(bin, statPlusSyst)
            pass
        chi2 = h_data.Chi2Test      (h_mc, "UW CHI2 P")
        KS   = h_data.KolmogorovTest(h_mc)
        ndf  = h_mc.GetXaxis().GetNbins() - 1

        # -- decorations
        c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
                "%s channel" % ('Jet' if isrjet else 'Photon')],
               qualifier=qualifier)
        c.xlabel('Large-#it{R} jet mass [GeV]')
        c.ylabel('Events / 5 GeV')

        if isrjet:
            c.ymin(150.)
        else:
            c.ymin(5.)
            pass

        if log:
            c.padding(0.45) # 0.50
        else:
            c.padding(0.35)
            pass
        pads[1].ylabel('Data / est.')
        pads[1].yline(1, linestyle=1, linecolor=ROOT.kBlack)
        c.region("SR", 0.8 * mass, 1.2 * mass, offset=0.07) #, offset=(0.14 if (mass == 150 and not isrjet) else 0.07))
        c.legend(sort=True) # ymax=0.87

        # @TEMP: For Laser's blinded plots
        '''
        print "===>", c.pads()[0].ylim()
        c.pads()[0].xline(100., linestyle=1, linewidth=1, linecolor=ROOT.kGray + 2, text='Blinded', text_align='TR')
        c.pads()[1].xline(100., linestyle=1, linewidth=1, linecolor=ROOT.kGray + 2,)
        '''

        c.log(log)

        stats_string = "KS prob. = %.3f  |  #chi^{2}/N_{dof} (prob.) = %.1f/%d (%.3f)" % (KS, chi2, ndf, ROOT.TMath.Prob(chi2, ndf))
        #c.latex(stats_string, 0.95 , 0.955, NDC=True, align=31, textsize=16)

        if args.save: c.save('plots/paperplot_%s_%dGeV%s.pdf' % ('jet' if isrjet else 'gamma', mass, '_log' if log else ''))
        if args.show: c.show()
        pass


    # tau21DDT cut optimisation
    # --------------------------------------------------------------------------
    mass = 160
    for isrjet in [True, False]: # [True, False]:
        
        if isrjet:
            f = ROOT.TFile(path + "hists_isrjet_cutoptimisation.root")
        else:
            f = ROOT.TFile(path + "hists_isrgamma_cutoptimisation.root")
            pass
              
        if not isrjet:
            names = ['h_tau21DDT_bkg', 'h_tau21DDT_sig%d' % mass, 'gr_improvements_sig%d' % mass]
        else:
            names = ['h_bkg_plot', 'h_sig_plot', 'Graph']
            pass
        
        histograms = dict()
        for name, title in zip(names, ['bkg', 'sig', 'impr']):
            hist = f.Get(name)
            try:
                hist.SetDirectory(0)
            except:
                # TGraph
                pass
            histograms[title] = hist
            pass
        
        f.Close()
 
        # Normalise graph
        gr = histograms['impr']
        N = gr.GetN()
        base, x, y = 0, ROOT.Double(0), ROOT.Double(-inf)
        gr.GetPoint(N - 1, x, y)
        base = float(y)
        for i in range(N):
            gr.GetPoint(i, x, y)
            gr.SetPoint(i, float(x), float(y)/base)
            pass
       
        # -- canvas
        c = ap.canvas(num_pads=1, batch=not args.show)
        pads = c.pads()
        
        # -- main pad
        c.hist(histograms['bkg'], fillcolor=compcolours['QCD'],  label="Incl. #gamma MC" if not isrjet else "Multijet MC")
        c.hist(histograms['sig'], linecolor=compcolours['WHad'], label="Z' (%d GeV)" % mass)
        
        # -- overlay
        o = ap.overlay(c, color=ROOT.kViolet)
        o.graph(histograms['impr'], linecolor=ROOT.kViolet, linewidth=2, option='L')
        
        # Get point of maximum improvement
        gr = histograms['impr']
        N = gr.GetN()
        ymax, xmax, x, y = -inf, 0, ROOT.Double(0), ROOT.Double(-inf)
        for i in range(N):
            gr.GetPoint(i, x, y)
            if float(y) > ymax:
                ymax = float(y)
                xmax = float(x)
                pass
            pass
        o.line(xmax, 0, xmax, ymax)
        
        # -- decorations
        c.text(["#sqrt{s} = 13 TeV,  L = 36.1 fb^{-1}",
                "%s channel" % ('Jet' if isrjet else 'Photon')],
               qualifier="Simulation " + qualifier)
        c.xlabel("Large-#it{R} jet #tau_{21}^{DDT}")
        c.ylabel("Events (normalised to unit)")
        o.label("Relative improvement in significance")
        c.legend(width=0.26)
        
        if args.show: c.show()
        if args.save: c.save('plots/paperplot_%s_%dGeV_cutoptimisation.pdf' % ('jet' if isrjet else 'gamma', mass))
        pass


    # Signal acceptance
    # --------------------------------------------------------------------------

    # Create canvas
    c = ap.canvas(batch=not args.show)

    for idx, isrjet in enumerate([False, True]):

        if isrjet:
            f = ROOT.TFile(path + "hists_isrjet_sigacc.root")
        else:
            f = ROOT.TFile(path + "hists_isrgamma_acceptance.root")
            pass

        if isrjet:
            acc    = f.Get('g_ptj')
            acceff = f.Get('g_tau')
        else:
            acc    = f.Get('acc')
            acceff = f.Get('acceff')
            pass

        # -- Draw graphs
        opts = {'linecolor': colours[idx], 'fillcolor': colours_light[idx], 'markercolor': colours[idx], 'linewidth': 2}
        c.graph(acc,     option=('A' if idx == 0 else '') + '3', label="%s channel" % ('Photon' if not isrjet else 'Jet'), legend_option='L', **opts)
        c.graph(acc,     option='LX', **opts)
        opts['fillcolor'] = colours[idx]
        c.graph(acceff, option='3L', fillstyle=3245, **opts)
        pass

    c.padding(0.5)
    c.xlabel("Signal mass point [GeV]")
    c.ylabel("Acceptance, acc. #times efficiency (Z' #rightarrow large-#it{R} jet + #gamma/jet)")
    c.text(["#sqrt{s} = 13 TeV"], qualifier="Simulation " + qualifier)
    c.legend(categories=[
            ("Acceptance",       {'fillcolor': ROOT.kGray,                        'option': 'LF'}),
            ("Acc. #times eff.", {'fillcolor': ROOT.kGray + 2, 'fillstyle': 3245, 'option': 'LF'}),
            ], reverse=True)
    if args.save: c.save('plots/paperplot_acceptance.pdf')
    if args.show: c.show()


    # Tau21 distributions
    # --------------------------------------------------------------------------

    for isrjet in [True, False]:

        # Create canvas
        c = ap.canvas(batch=not args.show)
        
        if isrjet:
            f = ROOT.TFile(path + "hists_isrjet_tau21DDTslices.root")
        else:
            f = ROOT.TFile(path + "hists_isrgamma_tau21distributions.root")
            pass
        
        if isrjet:
            slices = {
                'pt': [
                    ( 500,  700),
                    ( 700,  900),
                    ( 900, 1100),
                    ],
                'm': [
                    (100, 150),
                    (150, 200),
                    (200, 250),
                    ],
                }
        else:
            slices = {
                'pt': [
                    ( 200,  300),
                    ( 300,  500),
                    ( 500, 1000),
                    ],
                'm': [
                    (100, 150),
                    (150, 200),
                    (200, 250),
                    ],
                }
            pass

        if isrjet:
            histograms = [f.Get('h_%d_%d_%d_%d'      % (sl['pt'][0], sl['pt'][1], sl['m'][0], sl['m'][1])) for sl in dict_product(slices)]
        else:
            histograms = [f.Get('h_m_%d_%d_pt_%d_%d' % (sl['m'][0], sl['m'][1], sl['pt'][0], sl['pt'][1])) for sl in dict_product(slices)]
            pass

        category_names = ["[%.0f, %.0f] GeV" % (slices['m'][i][0], slices['m'][i][1]) for i in range(len(slices['m']))]
        
        # Draw
        for i, hist in enumerate(histograms):
            name = hist.GetName()
            if isrjet:
                pt1, pt2 = int(name.split('_')[1]), int(name.split('_')[2])
            else:
                pt1, pt2 = int(name.split('_')[5]), int(name.split('_')[6])
                pass
            #print "--", hist.GetName()
            label = "[%.0f, %.0f] GeV" % (pt1, pt2)
            c.hist(hist, fillstyle=0, linecolor=ROOT.kRed + (i%3) * 2, linestyle=1 + (i//3), label=label if i < 3 else None, normalise=True, option='HIST')
            pass

        # Decorations
        c.xlabel('Signal jet #tau_{21}^{DDT}')
        c.ylabel('Jets (normalised to unit)')
        c.text(["#sqrt{s} = 13 TeV", "%s selection" % ('Jet' if isrjet else 'Photon')], qualifier="Simulation " + qualifier)
        c.padding(0.50)

        ymax = 0.735
        c.legend(header='Jet p_{T} in:',
                 ymax=ymax)
        c.legend(header='Jet mass in:', categories=[
                        (category_names[idx], {'linestyle': 1 + idx}) for idx in range(len(category_names))
                        ],
                 ymax=ymax, xmin=0.19)

        if args.show: c.show()
        if args.save: c.save('plots/paperplot_%s_tau21distribution.pdf' % ('jet' if isrjet else 'gamma'))
        pass

    return
    
    
# Main function call.
if __name__ == '__main__':
    main()
    pass
