#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for producing harmonised plots for the di-jet + ISR paper

@file:   paperfigures.py
@author: Andreas Søgaard
@date:   13 June 2017
@email:  andreas.sogaard@cern.ch
"""

# Basic import(s)
import re, sys, glob, itertools

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
    import root_numpy
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

    qualifier = "" # "Internal" # "Internal"

    # Setup
    # -- Limit the number of axis digits to 4 from 5, to reduce clutter
    ROOT.TGaxis.SetMaxDigits(4)


    # Tau21 decorrelation plots
    # --------------------------------------------------------------------------

    basedir = {
        'isrjet':   '/eos/atlas/unpledged/group-wisc/users/lkaplan/dijetISR_all/dijetISR/plotsForNote/outputs/',
        'isrgamma': '/afs/cern.ch/user/a/asogaard/public/dijetisr/',
        }

    plottypes = ['correlation', 'decorrelation']
    for plottype in plottypes:

        if plottype == 'decorrelation':
            filenames = {
                'isrjet':   'hists_isrjet_decorrelation.root',
                'isrgamma': 'hists_isrgamma_decorrelation_tau21DDT_vs_m.root',
            }
        elif plottype == 'correlation':
            filenames = {
                'isrjet':   'hists_isrjet_correlation.root',
                'isrgamma': 'hists_isrgamma_decorrelation_tau21_vs_rhoDDT.root'
            }
        else:
            assert False
            pass

        for key, filename in filenames.iteritems():
            # Open file
            f = ROOT.TFile.Open(basedir[key] + filename)

            # Get list of profile names
            names = sorted(list(map(lambda tkey: tkey.GetName(), f.GetListOfKeys())))

            # Load profiles
            profiles = list()
            for name in names:
                
                # Load profie
                h = f.Get(name)
                h.SetDirectory(0)

                # Get pT range
                m = re.search('.*\_([\d]+)\_([\d]+)\_.*', name)
                pt = map(int, (m.group(1), m.group(2)))

                # Get signal status
                sig = ('sig' in name)
                profiles.append((h, sig, pt))
                pass

            # Create canvas
            c = ap.canvas(batch=not args.show)

            # Plot profiles
            colors = [ROOT.kBlack, compcolours['WHad'], compcolours['QCD']]
            for plot_signal in [False, True]:
                idx = 0
                for profile, sig, pt in profiles:
                    # Select desired class
                    if sig != plot_signal: continue

                    # Define style variables
                    color = colors[idx]
                    markerstyle = (24 if plot_signal else 20)

                    # Plot profile
                    c.plot(profile, markerstyle=markerstyle, markersize=1, markercolor=color, linecolor=color, option='PE', legend_option='PE', label="p_{T} #in  [%d, %d] GeV" % (pt[0], pt[1]))

                    # Increment counter
                    idx += 1
                    pass

                # Draw class-specific header
                header = "Signal, Z'(160 GeV):" if plot_signal else "Background:"
                if   plottype == 'decorrelation':
                    c.legend(header=header,
                             width=0.25, ymax = 0.38,
                             xmin=0.17 + 0.38 * plot_signal)
                elif plottype == 'correlation':
                    c.legend(header=header,
                             width=0.25, ymax = 0.38 + 0.35 * plot_signal,
                             xmin=0.17 + 0.38 * plot_signal)
                else:
                    assert False
                    pass
                pass

            # Re-draw axes on top
            c.pads()[0]._primitives[0].Draw('AXIS SAME')

            # Decorations
            c.ylim(0, 1)
            if   plottype == 'decorrelation':
                c.xlim(50, 300)
                c.xlabel("Large-#it{R} jet mass [GeV]")
                c.ylabel("#LT#tau_{21}^{DDT}#GT")
            elif plottype == 'correlation':
                c.xlim(1, 6)
                c.xlabel("Large-#it{R} jet mass [GeV]")
                c.ylabel("#LT#tau_{21}^{DDT}#GT")
            else:
                assert False
                pass

            c.text(["#sqrt{s} = 13 TeV",
                    "{} channel".format("Jet" if key == 'isrjet' else "Photon")],
                    qualifier="Simulation")

            # Show/save
            if args.show:
                c.show()
                pass

            if args.save:
                c.save('plots/fig2_{}_{}.pdf'.format(plottype, key.replace('isr', '')))
                pass
            pass
        pass


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
            # Store histogram
            histograms[name] = hist
            pass

        f.Close()

        # Combine  W/Z component
        histograms['WZHad'] = histograms['WHad'].Clone(histograms['WHad'].GetName().replace('W', 'WZ'))
        histograms['WZHad'].Add(histograms['ZHad'])

        # Compute systematics band
        histograms['syst'] = histograms['QCD'].Clone("h_mJ_syst")
        for bin in range(1, histograms['syst'].GetXaxis().GetNbins() + 1):
            nom  = histograms['QCD']     .GetBinContent(bin)
            up   = histograms['QCD_up']  .GetBinContent(bin)
            down = histograms['QCD_down'].GetBinContent(bin)
            histograms['syst'].SetBinError(bin, max(abs(up - nom), abs(nom - down)))
            pass

        # -- canvas
        c = ap.canvas(num_pads=2, fraction=0.45, batch=not args.show)
        pads = c.pads()

        # -- main pad
        h_mc   = c.stack(histograms['QCD'],  fillcolor=compcolours['QCD'],  label='Background est.')
        h_zhad = c.stack(histograms['ZHad'], fillcolor=compcolours['ZHad'], label='Z + %s' % ('jets' if isrjet else '#gamma'))
        h_whad = c.stack(histograms['WHad'], fillcolor=compcolours['WHad'], label='W + %s' % ('jets' if isrjet else '#gamma'))
        h_sum = c.hist(histograms['QCD'], fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                       label='Bkg. stat. uncert.')

        """ Using green dashed line to indicate syst error"
        h_bkg_up   = c.hist(histograms['QCD_up'],   linestyle=2, linecolor=compcolours['syst'], label='Syst. uncert.')
        h_bkg_down = c.hist(histograms['QCD_down'], linestyle=2, linecolor=compcolours['syst'])
        """
        h_syst = c.hist(histograms['syst'], fillstyle=1001, alpha=0.30, fillcolor=ROOT.kGray + 3, option='E2',
                        label='Bkg. syst. uncert.')


        h_data = c.plot (histograms['data'], label='Data', option='PE0X0', legend_option='PE')

        # -- diff pad
        pads[1].hist(histograms['WZHad'], fillcolor=compcolours['WHad'], option='HIST')
        pads[1].hist(histograms['ZHad'],  fillcolor=compcolours['ZHad'], option='HIST')
        #c.diff_plot((h_bkg_up,   histograms['QCD']), option='HIST')
        #c.diff_plot((h_bkg_down, histograms['QCD']), option='HIST')
        c.diff_plot((histograms['syst'],  histograms['QCD']), fillcolor=ROOT.kGray + 3, alpha=0.30, option='E2', uncertainties=False)
        c.diff_plot((histograms['data'],  histograms['QCD']), option='PE0X0')

        # -- statistical test
        h_mc.Add(h_zhad)
        h_mc.Add(h_whad)
        for bin in range(1, h_mc.GetXaxis().GetNbins()):
            stat = h_mc.GetBinError(bin)

            #nom  = h_sum     .GetBinContent(bin)
            #up   = h_bkg_up  .GetBinContent(bin)
            #down = h_bkg_down.GetBinContent(bin)
            #syst = max(abs(up-nom), abs(down-nom))
            syst = histograms['syst'].GetBinError(bin)

            h_mc.SetBinError(bin, np.sqrt( np.square(stat) + np.square(syst) ))
            pass
        chi2 = h_data.Chi2Test      (h_mc, "UW CHI2 P")
        KS   = h_data.KolmogorovTest(h_mc)
        ndf = h_mc.GetXaxis().GetNbins() - 1

        # -- decorations
        c.text(["#sqrt{s} = 13 TeV,  36.1 fb^{-1}",
                "W/Z validation",
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
        c.legend(ymax=0.891) # ..., sort=True)

        stats_string = "KS prob. = %.3f  |  #chi^{2}/N_{dof} (prob.) = %.1f/%d (%.3f)" % (KS, chi2, ndf, ROOT.TMath.Prob(chi2, ndf))
        #c.latex(stats_string, 0.95 , 0.96, NDC=True, align=31, textsize=16)

        if args.save: c.save('plots/paperplot_wz_%s.pdf' % ('jet' if isrjet else 'gamma'))
        if args.show: c.show()
        pass


    # Search plots
    # --------------------------------------------------------------------------
    for isrjet, mass, log in itertools.product([True, False], [140, 150, 160, 220], [True]):

        if isrjet and mass == 140: continue # Sample not available

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

        def normalise (hist):
            """Divide histogram by binwidth."""
            for bin in range(1, hist.GetXaxis().GetNbins() + 1):
                hist.SetBinContent(bin, hist.GetBinContent(bin) / float(hist.GetBinWidth(bin)))
                hist.SetBinError  (bin, hist.GetBinError  (bin) / float(hist.GetBinWidth(bin)))
                pass
            return hist

        histograms = dict()
        for name in names:
            hist = f.Get('h_mJ_' + name)
            try:
                hist.SetDirectory(0)

                if not isrjet:
                    hist = normalise(hist)
                    pass
                histograms[name] = hist

            except:
                # This mass point doesn't have a signal shape
                warning("Mass point %d doesn't have a signal shape." % mass)

                # Try to access the interpolated signal
                interp_path="/eos/atlas/user/a/asogaard/Analysis/2016/BoostedJetISR/StatsInputs/2017-08-06/"
                interp_filename="{}_{:d}.root".format("sig" if isrjet else "ISRgamma_signal", mass)
                interp_f = ROOT.TFile(interp_path + interp_filename, "READ")
                interp_treename = ("Signal_ISRjet_{:d}" if isrjet else "SIGNAL_{:d}_Nominal").format(mass)
                interp_t = interp_f.Get(interp_treename)
                array = root_numpy.tree2array(interp_t)
                fields = list(array.dtype.names)

                # -- Fill new histogram for interpolated signal
                bins = tf.config['massbins']
                c = ap.canvas(batch=True)
                warning("Fields : [{}, {}]".format(*fields))
                h_interp = c.hist(array[fields[0]], bins=bins, weights=array[fields[1]], display=False)
                warning("  Setting interpolated histogram as '{}'".format(name))
                histograms[name] = normalise(h_interp)
                pass
            pass

        f.Close()

        # Divide by bins width
        """
        if not isrjet:
            for hist in histograms.itervalues():
                for bin in range(1, hist.GetXaxis().GetNbins() + 1):
                    hist.SetBinContent(bin, hist.GetBinContent(bin) / float(hist.GetBinWidth(bin)))
                    hist.SetBinError  (bin, hist.GetBinError  (bin) / float(hist.GetBinWidth(bin)))
                    pass
                pass
            pass
        #"""

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

        if 'sig' in histograms:
            c.hist(histograms['sig'], linestyle=1, linecolor=ROOT.kViolet+2, label="Z' (%d GeV)" % mass)
        else:
            warning("Key 'sig' not in `histograms`")
            pass

        h_sum = c.getStackSum()
        h_sum = c.hist(h_sum, fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                       label='Bkg. stat. uncert.')

        """ Using green dashed line to indicate stat + syst error"
        h_bkg_up   = c.hist(histograms['QCD_up'],   linestyle=2, linecolor=compcolours['syst'],
                            label='Stat. #oplus syst.')
        h_bkg_down = c.hist(histograms['QCD_down'], linestyle=2, linecolor=compcolours['syst'])
        """
        h_bkg_var = h_sum.Clone("h_bkg_var")
        for bin in range(1, h_bkg_var.GetXaxis().GetNbins() + 1):
            nom  = h_sum                 .GetBinContent(bin)
            up   = histograms['QCD_up']  .GetBinContent(bin)
            down = histograms['QCD_down'].GetBinContent(bin)
            h_bkg_var.SetBinError(bin, max(abs(up - nom), abs(nom - down)))
            pass
        h_bkg_var = c.hist(h_bkg_var, fillstyle=1001, alpha=0.30, fillcolor=ROOT.kGray + 3, option='E2',
                           label='Bkg. stat. #oplus syst.')

        h_data = c.plot (histograms['data'], label='Data', option='PE0X0', legend_option='PE')

        # -- diff pad
        if isrjet: pads[1].ylim(0.98, 1.02)
        else:      pads[1].ylim(0.85, 1.15)

        pulls_stat_syst = c.ratio_plot((h_bkg_var,  h_sum), option='E2')
        pulls_bkg_stat  = c.ratio_plot((h_sum,              h_sum), option='E2')
        pulls_data      = c.ratio_plot((histograms['data'], h_sum), option='PE0X0', oob=True)

        # -- statistical test
        h_mc.Add(h_qcd)
        for bin in range(1, h_mc.GetXaxis().GetNbins()):
            nom  = h_mc      .GetBinContent(bin)
            #up   = h_bkg_up  .GetBinContent(bin)
            #down = h_bkg_down.GetBinContent(bin)
            #statPlusSyst = max(abs(up-nom), abs(down-nom))
            statPlusSyst = h_bkg_var.GetBinError(bin)

            # h_bkg_{up,down} are stats. oplus syst.
            h_mc.SetBinError(bin, statPlusSyst)
            pass
        chi2 = h_data.Chi2Test      (h_mc, "UW CHI2 P")
        KS   = h_data.KolmogorovTest(h_mc)
        ndf  = h_mc.GetXaxis().GetNbins() - 1

        # -- decorations
        c.text(["#sqrt{s} = 13 TeV,  36.1 fb^{-1}",
                 "%s channel" % ('Jet' if isrjet else 'Photon')
                ],
               qualifier=qualifier)
        c.xlabel('Large-#it{R} jet mass [GeV]')
        c.ylabel('Events / GeV') # ... / 5 GeV

        if isrjet:
            c.ymin(150. / 5.)
        else:
            c.ymin(10. / 5.)
            pass

        if log:
            c.padding(0.42) # 0.45 | 0.50
        else:
            c.padding(0.35)
            pass
        pads[1].ylabel('Data / est.')
        pads[1].yline(1, linestyle=1, linecolor=ROOT.kBlack)
        c.region("SR", 0.8 * mass, 1.2 * mass, offset=0.07) #, offset=(0.14 if (mass == 150 and not isrjet) else 0.07))
        c.legend(ymax=0.89) # ..., sort=True)

        c.log(log)


        stats_string = "KS prob. = %.3f  |  #chi^{2}/N_{dof} (prob.) = %.1f/%d (%.3f)" % (KS, chi2, ndf, ROOT.TMath.Prob(chi2, ndf))

        if args.save: c.save('plots/paperplot_%s_%dGeV%s.pdf' % ('jet' if isrjet else 'gamma', mass, '_log' if log else ''))
        if args.show: c.show()


        # Temp: Pre-fit pull plots for journal reviewer
        if mass not in [160, 220]: continue
        print "=" * 40
        print "isrjet: {}, mass: {}".format(isrjet, mass)
        pulls = list()
        diffs = list()
        errs  = list()
        CRs   = list()
        for ibin in range(1, pulls_stat_syst.GetXaxis().GetNbins() + 1):
            # @TODO: Only in signal region?
            CRs.append(abs(pulls_stat_syst.GetXaxis().GetBinCenter(ibin) - mass) / mass > 0.2)
            diff     = pulls_data.GetBinContent(ibin) - 1.
            err_bkg  = pulls_stat_syst.GetBinError(ibin)
            err_data = pulls_data.GetBinError(ibin)
            err      = np.sqrt( err_bkg**2 + err_data**2 )
            pull = diff / np.sqrt( err_bkg**2 + err_data**2 )

            diffs.append(diff)
            errs .append(err)
            pulls.append(pull)
            pass
        pulls_cr = [pull for cr, pull in zip(CRs, pulls) if cr]

        print "CR: mean +/- rms (eom) of pulls: {} +/- {} ({})".format(np.mean(pulls_cr), np.std(pulls_cr), np.std(pulls_cr) / np.sqrt(len(pulls_cr)))
        c1 = ap.canvas(batch=not args.show)
        bins = np.linspace(-3, 3, 6 * 2 + 1, endpoint=True)
        h1    = c1.hist(pulls, bins=bins, linecolor=ROOT.kRed, label='All bins')
        h1_cr = c1.hist(pulls_cr, bins=bins, linecolor=ROOT.kBlue, linestyle=2, label='Only CR bins')
        
        print "==== FULL:"
        h1   .Fit('gaus', 'V+')
        print "\n==== CR:"
        h1_cr.Fit('gaus', 'V+')

        c1.xlabel('Pre-fit pulls')
        c1.ylabel('Number of large-#it{R} jet mass bins')
        c1.text(["{} channel".format('Jet' if isrjet else 'Photon'), 
                 "Mass point: {} GeV".format(mass)], qualifier='Internal')
        c1.legend()
        c1.legend()
        c1.save('plots/temp_pulls_{}_{}.pdf'.format('jet' if isrjet else 'photon', mass))
        #c1.show()
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
        c.text(["#sqrt{s} = 13 TeV,  36.1 fb^{-1}",
                "%s channel" % ('Jet' if isrjet else 'Photon')],
               qualifier="Simulation " + qualifier)
        c.xlabel("Large-#it{R} jet #tau_{21}^{DDT}")
        c.ylabel("Fraction of events")
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


        print ""
        print "-" * 40
        print "ISR jet:", isrjet
        x, y = ROOT.Double(), ROOT.Double()
        print "Acceptance:"
        for i in range(acc.GetN()):
            acc.GetPoint(i,x,y)
            print "  (x,y) = ({},{})".format(x,y)
            pass
        print "Acceptance x efficiency:"
        for i in range(acceff.GetN()):
            acceff.GetPoint(i,x,y)
            print "  (x,y) = ({},{})".format(x,y)
            pass

        # -- Draw graphs
        opts = {'linecolor': colours[idx], 'fillcolor': colours_light[idx], 'markercolor': colours[idx], 'linewidth': 2}
        c.graph(acc,     option=('A' if idx == 0 else '') + '3', label="%s channel" % ('Photon' if not isrjet else 'Jet'), legend_option='L', **opts)
        c.graph(acc,     option='LX', **opts)
        opts['fillcolor'] = colours[idx]
        c.graph(acceff, option='3L', fillstyle=3245, **opts)
        pass

    c.log()
    #c.padding(0.4)
    #c.ymin(1.0E-04)
    c.ylim(1E-04, 1E+00)
    c.xlabel("m_{Z'} [GeV]")
    c.ylabel("Acceptance, acc. #times efficiency (Z' #rightarrow large-#it{R} jet + #gamma/jet)")
    c.text(["#sqrt{s} = 13 TeV"], qualifier="Simulation " + qualifier)
    c.legend(categories=[
            ("Acceptance",       {'fillcolor': ROOT.kGray,                        'option': 'LF'}),
            ("Acc. #times eff.", {'fillcolor': ROOT.kGray + 2, 'fillstyle': 3245, 'option': 'LF'}),
            ], reverse=True)#, ymax=0.84)
    if args.save: c.save('plots/paperplot_acceptance.pdf')
    if args.show: c.show()


    # Tau21 distributions
    # --------------------------------------------------------------------------

    for isrjet, ddt in itertools.product([True, False], ['', 'DDT']):

        # Create canvas
        c = ap.canvas(batch=not args.show)

        if isrjet:
            f = ROOT.TFile(path + "hists_isrjet_tau21{ddt}slices.root".format(ddt=ddt))
        else:
            f = ROOT.TFile(path + "hists_isrgamma_tau21{ddt}distributions.root".format(ddt=ddt))
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
        c.xlabel('Signal jet #tau_{21}^{%s}' % (ddt.upper()))
        c.ylabel('Fraction of events')
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
        if args.save: c.save('plots/paperplot_%s_tau21%sdistribution.pdf' % ('jet' if isrjet else 'gamma', ddt))
        pass

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
