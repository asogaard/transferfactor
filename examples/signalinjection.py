#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Script for performing signal injection tests.

@file:   signalinjection.py
@author: Andreas Søgaard
@date:   24 April 2017
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
    from transferfactor.utils import get_signal_DSID, get_histogram
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

parser = argparse.ArgumentParser(description='Perform signal injection test.')

parser.add_argument('--mass', dest='mass', type=int,
                    required=True,
                    help='Center of excluded mass window')
parser.add_argument('--show', dest='show', action='store_const',
                    const=True, default=False,
                    help='Show plots (default: False)')
parser.add_argument('--save', dest='save', action='store_const',
                    const=True, default=False,
                    help='Save plots (default: False)')
parser.add_argument('--inject', dest='inject', action='store_const',
                    const=True, default=False,
                    help='Inject signal (default: False)')
parser.add_argument('--toys', dest='toys', action='store_const',
                    const=True, default=False,
                    help='Use toys (default: False)')


# Main function.
def main ():

    # Parse command-line arguments
    args = parser.parse_args()


    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Get signal file
    sig_DSID = get_signal_DSID(args.mass, tolerance=10)
    if sig_DSID is None:
        return
    sig_file = 'objdef_MC_{DSID:6d}.root'.format(DSID=sig_DSID)

    # Load data
    files = glob.glob(tf.config['base_path'] + 'objdef_MC_3610*.root') + [tf.config['base_path'] + sig_file]

    if len(files) == 0:
        warning("No files found.")
        return

    data = loadData(files, tf.config['tree'], prefix=tf.config['prefix']) 
    info = loadData(files, tf.config['outputtree'], stop=1)

    # Scaling by cross section
    xsec = loadXsec(tf.config['xsec_file'])
    
    # Append new DSID field # @TODO: Make more elegant?
    data = append_fields(data, 'DSID', np.zeros((data.size,)), dtypes=int)
    for idx in info['id']:    
        msk = (data['id'] == idx) # Get mask of all 'data' entries with same id, i.e. from same file
        DSID = info['DSID'][idx]  # Get DSID for this file
        data['weight'][msk] *= xsec[DSID] # Scale by cross section x filter eff. for this DSID
        data['DSID']  [msk] = DSID        # Store DSID
        pass
    data['weight'] *= tf.config['lumi'] # Scale all events (MC) by luminosity

    # Check output.
    if data.size == 0:
        warning("No data was loaded.")
        return 

    # Compute new variables
    data = append_fields(data, 'logpt', np.log(data['pt']))

    # Separate out signal MC
    msk_sig  = (data['DSID'] == sig_DSID)
    msk_data = ~msk_sig

    print "DATA STATISTICS:", np.sum(data[msk_data]['weight'])

    signal = data[msk_sig]
    if not args.inject:
        # If we're not injecting signal, explicitly remove it from the 'data' array
        data = data[~msk_sig]
        pass


    # Toys
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    if args.toys:

        # Get masks
        msk_pass = tf.config['pass'](data)
        msk_fail = ~msk_pass
        
        # Create histograms
        if args.inject:
            pdf_pass = get_histogram(data, tf.config['params'], tf.config['axes'], mask=msk_pass & ~msk_sig)
            pdf_fail = get_histogram(data, tf.config['params'], tf.config['axes'], mask=msk_fail & ~msk_sig)
        else:
            pdf_pass = get_histogram(data, tf.config['params'], tf.config['axes'], mask=msk_pass)
            pdf_fail = get_histogram(data, tf.config['params'], tf.config['axes'], mask=msk_fail)
            pass

        # Smooth (only leading background)
        for _ in range(2):
            pdf_pass.Smooth()
            pdf_fail.Smooth()
            pass

        # Inject afterwards
        if args.inject:
            pdf_pass.Add(get_histogram(data, tf.config['params'], tf.config['axes'], mask=msk_pass & msk_sig))
            pdf_fail.Add(get_histogram(data, tf.config['params'], tf.config['axes'], mask=msk_fail & msk_sig))
        

        # Create p.d.f.s
        # -- Define variables
        rhoDDT = ROOT.RooRealVar('rhoDDT', 'rhoDDT', tf.config['axes'][0][0], tf.config['axes'][0][-1])
        logpt  = ROOT.RooRealVar('logpt',  'logpt',  tf.config['axes'][1][0], tf.config['axes'][1][-1])

        rhoDDT.setBins(len(tf.config['axes'][0]) - 1)
        logpt .setBins(len(tf.config['axes'][1]) - 1)
                    
        # -- Define histograms
        rdh_pass = ROOT.RooDataHist('rdh_pass', 'rdh_pass', ROOT.RooArgList(rhoDDT, logpt), pdf_pass)
        rdh_fail = ROOT.RooDataHist('rdh_fail', 'rdh_fail', ROOT.RooArgList(rhoDDT, logpt), pdf_fail)
        
        # -- Turn histograms into pdf's
        rhp_pass = ROOT.RooHistPdf('rhp_pass', 'rhp_pass', ROOT.RooArgSet(rhoDDT, logpt), rdh_pass)
        rhp_fail = ROOT.RooHistPdf('rhp_fail', 'rhp_fail', ROOT.RooArgSet(rhoDDT, logpt), rdh_fail)

        # Generate toys
        mult = 1.
        N_pass = int(np.sum(data['weight'][msk_pass]) * mult)
        N_fail = int(np.sum(data['weight'][msk_fail]) * mult)

        dtype = ['rhoDDT','logpt', 'tau21DDT', 'pt', 'm', 'weight']
        dtype = [(var, 'f8') for var in dtype]
        toys_pass = np.zeros(N_pass, dtype=dtype)
        toys_fail = np.zeros(N_fail, dtype=dtype)

        print "Generating toys (pass: %d, fail: %d)" % (N_pass, N_fail)
        rds_pass = rhp_pass.generate(ROOT.RooArgSet(rhoDDT, logpt), N_pass, True, False)
        rds_fail = rhp_fail.generate(ROOT.RooArgSet(rhoDDT, logpt), N_fail, True, False)

        for idx in range(N_pass):
            toys_pass['rhoDDT']  [idx] = rds_pass.get(idx).getRealValue('rhoDDT')
            toys_pass['logpt']   [idx] = rds_pass.get(idx).getRealValue('logpt')
            toys_pass['pt']      [idx] = np.exp(toys_pass['logpt'][idx])
            toys_pass['m']       [idx] = np.sqrt(np.exp(toys_pass['rhoDDT'][idx]) * toys_pass['pt'][idx] * 1.)
            toys_pass['weight']  [idx] = 1. / float(mult)
            toys_pass['tau21DDT'][idx] = 0.
            pass

        for idx in range(N_fail):
            toys_fail['rhoDDT']  [idx] = rds_fail.get(idx).getRealValue('rhoDDT')
            toys_fail['logpt']   [idx] = rds_fail.get(idx).getRealValue('logpt')
            toys_fail['pt']      [idx] = np.exp(toys_fail['logpt'][idx])
            toys_fail['m']       [idx] = np.sqrt(np.exp(toys_fail['rhoDDT'][idx]) * toys_fail['pt'][idx] * 1.)
            toys_fail['weight']  [idx] = 1. / float(mult)
            toys_fail['tau21DDT'][idx] = 1.
            pass

        data = np.concatenate((toys_pass, toys_fail)) # ???
        pass


    # Transfer factor
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    calc = tf.calculator(data=data, config=tf.config) # Using default configuration
    calc.mass = args.mass
    calc.fullfit()

    # Pass/fail masks
    msk_data_pass = tf.config['pass'](data)
    msk_data_fail = ~msk_data_pass
    msk_sig_pass  = tf.config['pass'](signal)
    msk_sig_fail  = ~msk_sig_pass
    
    print "  -- Computing data weights"
    w_nom, w_up, w_down = calc.fullweights(data[msk_data_fail])
    print "  -- Computing signal weights"
    w_sig, _, _ = calc.fullweights(signal[msk_sig_fail])
    print "  -- Final fit done"
    if args.show or args.save: calc.plot(show=args.show, save=args.save, prefix='plots/signalinjection_%s%s_' % ("toys_" if args.toys else "", "injected" if args.inject else "notinjected"))
    

    # Performing signal injection test
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    
    if True or args.show or args.save:
        
        bestfit_mu = None
        
        for mu, fit, prefit, subtract in zip([0, 1, 1, None],
                                             [False, False, True,  False],
                                             [True,  True,  True,  False],
                                             [True,  True,  False, True]):
            
            if not prefit:
                mu = bestfit_mu[0]
                pass


            # Plotting
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            
            c = ap.canvas(num_pads=2, batch=not args.show)
            p0, p1 = c.pads()

            # -- Histograms: Main pad
            bins = tf.config['massbins']

            h_bkg      = c.hist(data['m'][msk_data_fail], bins=bins, weights=data['weight'][msk_data_fail] * w_nom,  display=False)
            h_bkg_up   = c.hist(data['m'][msk_data_fail], bins=bins, weights=data['weight'][msk_data_fail] * w_up,   display=False)
            h_bkg_down = c.hist(data['m'][msk_data_fail], bins=bins, weights=data['weight'][msk_data_fail] * w_down, display=False)
            
            h_sig = c.hist(signal['m'][msk_sig_pass], bins=bins, weights=signal['weight'][msk_sig_pass],         scale=mu, display=False)
            h_sfl = c.hist(signal['m'][msk_sig_fail], bins=bins, weights=signal['weight'][msk_sig_fail] * w_sig, scale=mu, display=False)
            
            if not fit:
                h_bkg     .Add(h_sfl, -1) # Subtracting signal
                h_bkg_up  .Add(h_sfl, -1) # --
                h_bkg_down.Add(h_sfl, -1) # --
                pass

            h_bkg = c.stack(h_bkg,
                            fillcolor=ROOT.kAzure + 7, 
                            label='Background pred.')
            h_sig = c.stack(h_sig,
                            fillcolor=ROOT.kOrange + 7,
                            label="Z' (#mu = %s)" % ("%.0f" % mu if prefit else "%.2f #pm %.2f" % (mu, bestfit_mu[1])))
            
            h_sum = h_bkg
            h_sum = c.hist(h_sum, 
                           fillstyle=3245, fillcolor=ROOT.kGray + 3, option='E2',
                           label='Stat. uncert.')

            h_bkg_up   = c.hist(h_bkg_up,
                                linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST',
                                label='Syst. uncert.')
            h_bkg_down = c.hist(h_bkg_down,
                                linecolor=ROOT.kGreen + 1, linestyle=2, option='HIST')
        
            h_data = c.plot(data['m'][msk_data_pass], bins=bins, weights=data['weight'][msk_data_pass], display=False)
            h_data = c.plot(h_data,
                            label='Pseudo-data')

            # -- Histograms: Ratio pad
            c.ratio_plot((h_sig,      h_sum), option='HIST', offset=1, fillcolor=ROOT.kOrange + 7)
            c.ratio_plot((h_sum,      h_sum), option='E2')
            c.ratio_plot((h_bkg_up,   h_sum), option='HIST')
            c.ratio_plot((h_bkg_down, h_sum), option='HIST')
            c.ratio_plot((h_data,     h_sum))

            # -- Text
            c.text(["#sqrt{s} = 13 TeV,  L = %s fb^{-1}" % tf.config['lumi'],
                    "Sherpa incl. #gamma MC",
                    "Trimmed anti-k_{t}^{R=1.0} jets",
                    "ISR #gamma selection",
                    "Window: %d GeV #pm %d %%" % (args.mass, 20.),
                    ("Signal" if args.inject else "No signal") + " injected",
                   ] + (["Using toys"] if args.toys else []), qualifier='Simulation Internal')

            # -- Axis labels
            c.xlabel('Signal jet mass [GeV]')
            c.ylabel('Events')
            p1.ylabel('Data / Est.')

            # -- Log
            c.log()            

            # -- Axis limits
            p0.padding(0.55)
            p1.ylim(0.7, 1.3)

            # -- Line(s)
            p1.yline(1.0)

            # -- Region(s)        
            c.region("SR", 0.8 * args.mass, 1.2 * args.mass)
            
            c.legend()
            if args.show and not fit: c.show()
            if args.save and not fit: c.save('plots/signalinjection_%s%dGeV_pm%d_%s_%s.pdf' % (
                    "toys_" if args.toys else "",
                    args.mass, 
                    20.,
                    ('prefit_mu%d' % mu if prefit else 'postfit'),
                    ('injected' if args.inject else 'notinjected')
                    ))

            
            
            # Fitting
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

            if fit:

                bestfit_mu = list()

                hs_save = [
                    h_bkg_down.Clone('h_save_down'),
                    h_bkg     .Clone('h_save_nom'),
                    h_bkg_up  .Clone('h_save_up'),
                    ]
                
                for variation in range(3):
                    
                    print "Variation: " + ("Nominal" if variation == 1 else ("Up" if variation == 0 else "Down"))

                    # Get correct histogram fore this variation
                    h_bkg_use = hs_save[variation]
                    
                    # -- Define jet mass variable
                    mJ = ROOT.RooRealVar('mJ', 'mJ', 50, 300)
                    mJ.setBins(50)
                    
                    # -- Define histograms
                    rdh_bkg = ROOT.RooDataHist('rdh_bkg', 'rdh_bkg', ROOT.RooArgList(mJ), h_bkg_use)
                    rdh_sig = ROOT.RooDataHist('rdh_sig', 'rdh_sig', ROOT.RooArgList(mJ), h_sig)
                    rdh_sfl = ROOT.RooDataHist('rdh_sfl', 'rdh_sfl', ROOT.RooArgList(mJ), h_sfl)
                    
                    # -- Turn histograms into pdf's
                    rhp_bkg = ROOT.RooHistPdf('rhp_bkg', 'rhp_bkg', ROOT.RooArgSet(mJ), rdh_bkg)
                    rhp_sig = ROOT.RooHistPdf('rhp_sig', 'rhp_sig', ROOT.RooArgSet(mJ), rdh_sig)
                    rhp_sfl = ROOT.RooHistPdf('rhp_sfl', 'rhp_sfl', ROOT.RooArgSet(mJ), rdh_sfl)
                    
                    # -- Define integrals as constants
                    n_bkg = ROOT.RooRealVar('n_bkg', 'n_bkg', h_bkg_use.Integral())
                    n_sig = ROOT.RooRealVar('n_sig', 'n_sig', h_sig.Integral())
                    n_sfl = ROOT.RooRealVar('n_sfl', 'n_sfl', h_sfl.Integral())
                    
                    # -- Define signal strength and constant(s)
                    mu   = ROOT.RooRealVar('mu',   'mu', 1, 0, 5)
                    neg1 = ROOT.RooRealVar('neg1', 'neg1', -1)
                    
                    # -- Define fittable normalisation factors
                    c_bkg = ROOT.RooFormulaVar('c_bkg', 'c_bkg', '@0',           ROOT.RooArgList(n_bkg))
                    c_sig = ROOT.RooFormulaVar('c_sig', 'c_sig', '@0 * @1',      ROOT.RooArgList(mu, n_sig))
                    c_sfl = ROOT.RooFormulaVar('c_sfl', 'c_sfl', '@0 * @1 * @2', ROOT.RooArgList(neg1, mu, n_sfl))
                    
                    # -- Construct combined pdf
                    pdf = ROOT.RooAddPdf('pdf', 'pdf', ROOT.RooArgList(rhp_bkg, rhp_sig, rhp_sfl), ROOT.RooArgList(c_bkg, c_sig, c_sfl))
                    
                    # -- Construct data histogram
                    rdh_data = ROOT.RooDataHist('rdh_data', 'rdh_data', ROOT.RooArgList(mJ), h_data)
                    
                    # -- Fit pdf to data histogram
                    pdf.chi2FitTo(rdh_data, ROOT.RooLinkedList(), ROOT.Strategy(0))
                    
                    print "Best fit mu: %.3f +/- %.3f" % (mu.getValV(), mu.getError())
                    bestfit_mu.append( (mu.getValV(), mu.getError()) )
                    pass

                bestfit_mu = bestfit_mu[1][0], np.sqrt(np.power(abs(bestfit_mu[0][0] - bestfit_mu[2][0]) / 2., 2.) + np.power(bestfit_mu[1][1], 2.))
                pass

            pass
        
        pass
    
    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
