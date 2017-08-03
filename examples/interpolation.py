#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for ...

@file:   interpolation.py
@author: Andreas Søgaard
@date:   14 June 2017
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
    
    import matplotlib.pyplot as plt
except ImportError:
    print "WARNING: One or more scientific python packages were not found. If you're in lxplus, try running:"
    print "         $ source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh"
    sys.exit()
    pass

# Local import(s)
try:
    import transferfactor as tf
    from transferfactor.utils import make_directories, check_make_dir
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

parser = argparse.ArgumentParser(description='Perform closure test.')

parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--isrjet', dest='isrjet', action='store_const',
                    const=True, default=False,
                    help='Run on ISR jet channel (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')

# Main function.
def main ():

    # Parse command-line arguments
    args = parser.parse_args()

    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Get files
    if not args.isrjet:
        files = glob.glob(tf.config['base_path'] + 'objdef_MC_30836*.root')
    else:
        #files = glob.glob('/afs/cern.ch/user/l/lkaplan/public/forAndreas/signals/sig_*.root')
        files = glob.glob('/afs/cern.ch/user/l/lkaplan/public/forAndreas/signals_fatjetunc/sig_*.root')
        pass

    if len(files) == 0:
        warning("No files found.")
        return

    # Get names of all available systematic variations
    f = ROOT.TFile(files[0], 'READ')
    if not args.isrjet:
        f.cd("BoostedJet+ISRgamma")
        pass
    variations = [key.GetName() for key in ROOT.gDirectory.GetListOfKeys()]
    f.Close()

    colours = [ROOT.kViolet + 7, ROOT.kAzure + 7, ROOT.kTeal, ROOT.kSpring - 2, ROOT.kOrange - 3, ROOT.kPink]

    # Load data
    if not args.isrjet:
        data = {var: loadData(files, tf.config['finaltree'].replace('Nominal', var), prefix=tf.config['prefix']) for var in variations}
        info = {var: loadData(files, tf.config['outputtree'].replace('Nominal', var), stop=1) for var in variations}
        
        # Scaling by cross section
        xsec = loadXsec(tf.config['xsec_file'])
        
        # Append new DSID field # @TODO: Make more elegant?
        for var in variations:
            data[var] = append_fields(data[var], 'DSID', np.zeros((data[var].size,)), dtypes=int)
            for idx in info[var]['id']:    
                msk = (data[var]['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
                DSID = info[var]['DSID'][idx]  # Get DSID for this file
                data[var]['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
                data[var]['DSID']  [msk] = DSID        # Store DSID
                pass
            data[var]['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity
            pass

    else:
        data = {var: loadData(files, var) for var in variations}
        for var in variations:
            data[var] = append_fields(data[var], 'DSID', np.zeros((data[var].size,)), dtypes=int)
            for idx in list(set(data[var]['id'])):
                msk = (data[var]['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
                data[var]['DSID'][msk] = idx
                pass
            pass
        pass


    # Check output.
    if not args.isrjet and data['Nominal'].size == 0:
        warning("No data was loaded.")
        return 

    nominal = 'Nominal' if not args.isrjet else 'Signal_ISRjet'
    m       = 'm' if not args.isrjet else 'mJ'

    # Perform interpolation of signal mass peaks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    closure = False

    # Store mass points
    massPoints = np.array([100, 130, 160, 190, 220])
    massVec = ROOT.TVectorD(len(massPoints))
    for i in range(len(massPoints)):
        massVec[i] = massPoints[i]
        pass

    # Prepare workspaces
    workspace = ROOT.RooWorkspace()

    # Create variables
    mZ = workspace.factory('mZ[50,300]') # this is our continuous interpolation parameter

    mJ = workspace.factory('mJ[50,300]')
    mJ.setBins(50)
    frame = mJ.frame()


    # Create fit model
    for i in range(len(massPoints)):
        workspace.factory('Gaussian::signal{i}(mJ, mu_sig_{i}[{mean},0,300], sigma_sig_{i}[{width},0,50])'.format(i=i, mean=massPoints[i], width=massPoints[i]*0.1))
        if not args.isrjet:
            workspace.factory('Gaussian::background{i}(mJ, mu_bkg_{i}[0,0,0], sigma_bkg_{i}[100,50,500])'.format(i=i))
        else:
            workspace.factory('Exponential::background{i}(mJ, tau_bkg_{i}[-0.01,-100,0])'.format(i=i))
            pass

        workspace.factory('SUM::model{i}(norm_sig_{i}[2,0,10]*signal{i}, norm_bkg_{i}[1,1,1]*background{i})'.format(i=i))
        pass

    # Create p.d.f.s
    c = ap.canvas(batch=not args.show)
    bins = np.linspace(50, 300, 50 + 1, endpoint=True)

    pdfs = ROOT.RooArgList(*[workspace.pdf('model{i}'.format(i=i)) for i in range(len(massPoints))])
    integrals = np.zeros((len(massPoints),), dtype=float)
    events    = np.zeros((len(massPoints),), dtype=float)
    for idx in range(len(massPoints)):
        DSID = idx + (308363 if not args.isrjet else 0)
        print "==> %d" % DSID

        # Create data histogram for fitting
        msk = (data[nominal]['DSID'] == DSID)        
        hist = c.hist(data[nominal][m][msk], bins=bins, weights=data[nominal]['weight'][msk], normalise=True, display=False)
        integrals[idx] = np.sum(data[nominal]['weight'][msk])
        events   [idx] = np.sum(msk)

        # Add minimal error to empty bins
        emax = np.finfo(float).max
        for bin in range(1, hist.GetXaxis().GetNbins() + 1):
            if hist.GetBinError(bin) == 0: hist.SetBinError(bin, emax)
            pass

        # Create RooFit histogram and p.d.f.
        rdh = ROOT.RooDataHist('rdh', 'rdh', ROOT.RooArgList(mJ), hist)
        rhp = ROOT.RooHistPdf ('rhp', 'rhp', ROOT.RooArgSet (mJ), rdh)

        # Fit data with model
        pdfs[idx].chi2FitTo(rdh, ROOT.RooLinkedList())

        # Plot data histogram and fit
        rhp      .plotOn(frame, ROOT.RooFit.LineColor(colours[idx]), ROOT.RooFit.LineStyle(1), ROOT.RooFit.LineWidth(3))
        pdfs[idx].plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlack), ROOT.RooFit.LineStyle(2))
        pass

    setting = ROOT.RooMomentMorph.Linear
    morph = ROOT.RooMomentMorph('morph', 'morph', mZ, ROOT.RooArgList(mJ), pdfs, massVec, setting)
    getattr(workspace,'import')(morph) # work around for morph = w.import(morph)
    morph.Print('v')

    # Make plots of interpolated p.d.f.s
    interpolationMassPoints = np.linspace(100, 220, (220 - 100) / 10 + 1, endpoint=True)
    interpolationMassPoints = sorted(list(set(interpolationMassPoints) - set(massPoints)))
    """ @TEMP: BEGIN """
    # -- Interpolate logarithmically between integrals
    interpolationIntegrals = np.exp(np.interp(interpolationMassPoints, massPoints, np.log(integrals)))
        
    # -- Interpolate linearly between event counts
    interpolationEvents    =        np.interp(interpolationMassPoints, massPoints, events).astype(int)   
    """ @TEMP: END """

    #for i, (n, integral, mass) in enumerate(zip(interpolationEvents, interpolationIntegrals, interpolationMassPoints)):
    for i, mass in enumerate(interpolationMassPoints):
        print "=" * 80
        mZ.setVal(mass)
        mZ.Print()
        morph.Print()
        morph.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed - 4), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineWidth(1))
        pass

    # Interpolation closure
    if closure:
        for i in [2]:
            DSID = 308363 + i
            print "==> %d" % DSID
            
            # Create data histogram for fitting
            msk = (data[nominal]['DSID'] == DSID)        
            hist = c.hist(data[nominal][m][msk], bins=bins, weights=data[nominal]['weight'][msk], normalise=True, display=False)
            
            # Add minimal error to empty bins
            emax = np.finfo(float).max
            for bin in range(1, hist.GetXaxis().GetNbins() + 1):
                if hist.GetBinError(bin) == 0: hist.SetBinError(bin, emax)
                pass
            
            # Create RooFit histogram and p.d.f.
            rdh = ROOT.RooDataHist('rdh', 'rdh', ROOT.RooArgList(mJ), hist)
            rhp = ROOT.RooHistPdf ('rhp', 'rhp', ROOT.RooArgSet (mJ), rdh)
            
            # Plot data histogram and fit
            rhp    .plotOn(frame, ROOT.RooFit.LineColor(1), ROOT.RooFit.LineStyle(1))
            pass
        
        for massval in np.linspace(135, 185, 9, endpoint=True):
            mZ.setVal(massval)
            morph.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineWidth(3 if massval == 160 else 1))    
            pass
        pass

    # Draw frame
    if args.show or args.save:
        c = ap.canvas(batch=not args.show)
        c.pads()[0]._bare().cd()
        frame.Draw()
        frame.GetXaxis().SetTitle("Large-#it{R} jet mass [GeV]")
        frame.GetYaxis().SetTitle("Signal p.d.f.")
        frame.GetYaxis().SetRangeUser(0, 0.22)
        c.text(["#sqrt{s} = 13 TeV", "ISR %s channel" % ('#gamma' if not args.isrjet else 'jet')], qualifier="Simulation Internal")
        c.legend(header="Signal MC shape:",
                 categories= [ ("Z' (%d GeV)" % mass, {'linecolor': colours[idx], 'linewidth': 3, 'option': 'L'}) for idx, mass in enumerate(massPoints)]
                + [
                ("Fitted  shape",      {'linecolor': ROOT.kBlack,   'linestyle': 2, 'linewidth': 2, 'option': 'L'}),
                ("Interpolated shape", {'linecolor': ROOT.kRed - 4, 'linestyle': 2, 'linewidth': 1, 'option': 'L'}),
                ])
        c.ylim(0, 0.22)
        if args.save: c.save('plots/interpolation_frame_isr%s.pdf' % ('gamma' if not args.isrjet else 'jet'))
        if args.show: c.show()
        pass

    # Write outputs
    if args.save:
        check_make_dir('output')
        
        # -- Interpolate logarithmically between integrals
        interpolationIntegrals = np.exp(np.interp(interpolationMassPoints, massPoints, np.log(integrals)))
        
        # -- Interpolate linearly between event counts
        interpolationEvents    =        np.interp(interpolationMassPoints, massPoints, events).astype(int)
        
        mult = 100
        for n, integral, mass in zip(mult * interpolationEvents, interpolationIntegrals, interpolationMassPoints):
            
            # Generate mass- and weight dataset for interpolated mass point
            weight = integral / float(n)
            mZ.setVal(mass)
            dataset = morph.generate(ROOT.RooArgSet(mJ), n)
            vector_m = np.array([dataset.get(idx).getRealValue('mJ') for idx in range(n)])
            
            # Get closest mass points on either side
            idx1 = np.where(massPoints - mass < 0)[0][-1]
            idx2 = np.where(massPoints - mass > 0)[0][0]
            DSID1 = idx1 + (308363 if not args.isrjet else 0)
            DSID2 = idx2 + (308363 if not args.isrjet else 0)
            
            # Open output file
            DSID = int("200%03d" % mass)
            if not args.isrjet:
                filename = 'objdef_MC_{DSID:6d}.root'.format(DSID=DSID)
            else:
                filename = 'sig_%d.root' % mass
                pass
            outputdir = '/eos/atlas/user/a/asogaard/Analysis/2016/BoostedJetISR/StatsInputs/2017-07-19/'
            #output = ROOT.TFile('output/' + filename, 'RECREATE')
            output = ROOT.TFile(outputdir + filename, 'RECREATE')

            print "Writing %d events to file '%s'" % (n, output.GetName())
            ROOT.TH1.AddDirectory(0)
            
            # Get nominal histograms for each of the two closes mass points
            
            msk1 = data[nominal]['DSID'] == DSID1
            msk2 = data[nominal]['DSID'] == DSID2
            
            hist_nom1 = c.hist(data[nominal][m][msk1], bins=bins, weights=data[nominal]['weight'][msk1], display=False)
            hist_nom2 = c.hist(data[nominal][m][msk2], bins=bins, weights=data[nominal]['weight'][msk2], display=False)
            
            # Loop systematic variations
            for var in variations:

                # Reset for each variation!
                vector_w = np.ones((n,)) * weight
                
                print "  Writing arrays to file: %s" % var
                
                msk1 = data[var]['DSID'] == DSID1
                msk2 = data[var]['DSID'] == DSID2
                
                hist_var1 = c.hist(data[var][m][msk1], bins=bins, weights=data[var]['weight'][msk1], display=False)
                hist_var2 = c.hist(data[var][m][msk2], bins=bins, weights=data[var]['weight'][msk2], display=False)
                
                # Plot variations
                c2 = ap.canvas(num_pads=2, batch=not args.show)
                p0, p1 = c2.pads()
                c2.hist(vector_m, bins=bins, weights=vector_w, label="%d GeV (interp.)" % mass)
                c2.hist(hist_nom1, linecolor=ROOT.kRed,  label="%d GeV" % massPoints[idx1])
                c2.hist(hist_nom2, linecolor=ROOT.kBlue, label="%d GeV" % massPoints[idx2])
                c2.hist(hist_var1, linecolor=ROOT.kRed,  linestyle=2)
                c2.hist(hist_var2, linecolor=ROOT.kBlue, linestyle=2)
                
                # -- Get the variation/nominal ratio for each closest masspoint
                hist_rel1 = c2.ratio_plot((hist_var1, hist_nom1), linecolor=ROOT.kRed,  linestyle=1, option='HIST')
                hist_rel2 = c2.ratio_plot((hist_var2, hist_nom2), linecolor=ROOT.kBlue, linestyle=1, option='HIST')
                
                # -- Get the number of bins to shift the systematic variation of each closest masspoint
                binwidth = bins[1] - bins[0]
                shift1 = int(np.round((mass - massPoints[idx1]) / binwidth))
                shift2 = int(np.round((mass - massPoints[idx2]) / binwidth))
                arr_rel1 = hist2array(hist_rel1)
                arr_rel2 = hist2array(hist_rel2)
                
                # -- Shift the variations by the appropriate number of bins by _rolling_
                arr_relshift1 = np.roll(arr_rel1, shift1)
                arr_relshift2 = np.roll(arr_rel2, shift2)
                arr_relshift1[:shift1] = 1
                arr_relshift2[shift2:] = 1
                
                # -- Get weights for the averaging of the two shifted systematics
                w1 = 1. / float(abs(shift1))
                w2 = 1. / float(abs(shift2))
                
                # -- Take the weighted average
                arr_relshift = (w1 * arr_relshift1 + w2 * arr_relshift2) / (w1 + w2)
                
                hist_relshift1 = array2hist(arr_relshift1, hist_rel1.Clone("hist_relshift1"))
                hist_relshift2 = array2hist(arr_relshift2, hist_rel2.Clone("hist_relshift2"))
                hist_relshift  = array2hist(arr_relshift,  hist_rel1.Clone("hist_relshift"))
                
                p1.hist(hist_relshift1, linecolor=ROOT.kRed,   linestyle=2)
                p1.hist(hist_relshift2, linecolor=ROOT.kBlue,  linestyle=2)
                p1.hist(hist_relshift,  linecolor=ROOT.kBlack, linestyle=1)
                
                # -- Compute the jet-by-jet weights for the current variation
                vector_w_var = vector_w
                for bin in range(1, hist_relshift.GetXaxis().GetNbins() + 1):
                    # Get mask for all jet in 'bin'
                    msk = (vector_m >= hist_relshift.GetXaxis().GetBinLowEdge(bin)) & (vector_m < hist_relshift.GetXaxis().GetBinUpEdge(bin))
                    vector_w_var[msk] *= hist_relshift.GetBinContent(bin)
                    pass
                
                c2.hist(vector_m, bins=bins, weights=vector_w_var, linecolor=ROOT.kBlack, linestyle=2)
                
                p1.ylim(0.5, 1.5)
                p1.yline(1)
                c2.xlabel('Large-#it{R} jet mass [GeV]')
                c2.ylabel('Events / 5 GeV')
                p1.ylabel('Variation / nominal')
                c2.legend(header="Signal mass point:", categories=[
                        ('Nominal', {}),
                        (var,       {'linestyle': 2}),
                        ])
                c2.text(["#sqrt{s} = 13 TeV",
                         "ISR %s channel" % ('gamma' if not args.isrjet else 'jet'),
                         ], qualifier='Simulation Internal')
                if args.save: c2.save('plots/interpolation_shift_isr%s_%s_%dGeV.pdf' % ('gamma' if not args.isrjet else 'jet', var, mass))
                if args.show: c2.show()
                
                
                # Prepare DISD and isMC vectors
                vector_DSID = np.ones_like(vector_w_var) * DSID
                vector_isMC = np.ones_like(vector_w_var).astype(bool)
                
                jetmassname = (tf.config['prefix'] + 'm' if not args.isrjet else 'mJ')
                
                array1 = np.array(zip(vector_m, vector_w_var),
                                  dtype = [(jetmassname, np.float64),
                                           ('weight',    np.float64)])
                
                array2 = np.array(zip(vector_DSID, vector_isMC),
                                  dtype = [('DSID', np.uint32),
                                           ('isMC', np.bool_)])
                
                # Mass and weight branch
                if not args.isrjet:
                    treename1 = tf.config['finaltree'].replace('Nominal', var)
                else:
                    treename1 = var.replace('ISRjet', 'ISRjet_%d' % mass)
                    pass
                make_directories('/'.join(treename1.split('/')[:-1]), fromDir=output)
                tree1 = ROOT.TTree(treename1.split('/')[-1], "")
                array2tree(array1, tree=tree1)
                
                # outputTree
                if not args.isrjet:
                    treename2 = tf.config['outputtree'].replace('Nominal', var)
                    make_directories('/'.join(treename2.split('/')[:-1]), fromDir=output)
                    tree2 = ROOT.TTree(treename2.split('/')[-1], "")
                    array2tree(array2, tree=tree2)
                    pass
                
                output.Write()
                pass
            
            output.Close()
            pass

        pass

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
