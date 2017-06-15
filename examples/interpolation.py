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
except ImportError:
    print "WARNING: One or more scientific python packages were not found. If you're in lxplus, try running:"
    print "         $ source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh"
    sys.exit()
    pass

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

    # Parse command-line arguments
    #args = parser.parse_args()


    # Setup.
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Load data
    files = glob.glob(tf.config['base_path'] + 'objdef_MC_30836*.root')

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


    # Perform interpolation of signal mass peaks
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    # Store mass points
    massPoints = np.array([100, 130, 160, 190, 220])
    massVec = ROOT.TVectorD(len(massPoints))
    for i in range(5):
        massVec[i] = massPoints[i]
        pass

    # Prepare workspaces
    workspace = ROOT.RooWorkspace()

    # Create variables
    mJ = workspace.factory('mJ[50,300]')
    mZ = workspace.factory('mZ[50,300]') # this is our continuous interpolation parameter
    mJ.setBins(50)
    frame = mJ.frame()

    # Create fit model
    for i in range(len(massPoints)):
        workspace.factory('Gaussian::g{i}(mJ, mu{i}[{mean},0,300], sigma{i}[{width},0,50])'.format(i=i, mean=massPoints[i], width=massPoints[i]*0.1))
        workspace.factory('Gaussian::gb{i}(mJ, mu_b{i}[0,0,50], sigma_b{i}[100,50,150])'.format(i=i))
        workspace.factory('SUM::model{i}(s{i}[50,0,100]*g{i}, b{i}[100,0,1000]*gb{i})'.format(i=i))
        pass

    # Create p.d.f.s
    c = ap.canvas()
    bins = np.linspace(50, 300, 50 + 1, endpoint=True)

    pdfs = ROOT.RooArgList(*[workspace.pdf('model{i}'.format(i=i)) for i in range(len(massPoints))])
    for i in range(len(massPoints)):
        DSID = 308363 + i
        print "==> %d" % DSID

        # Create data histogram for fitting
        msk = (data['DSID'] == DSID)        
        hist = c.hist(data['m'][msk], bins=bins, weights=data['weight'][msk], normalise=True, display=False)

        # Add minimal error to empty bins
        emax = np.finfo(float).max
        for bin in range(1, hist.GetXaxis().GetNbins() + 1):
            if hist.GetBinError(bin) == 0: hist.SetBinError(bin, emax)
            pass

        # Create RooFit histogram and p.d.f.
        rdh = ROOT.RooDataHist('rdh', 'rdh', ROOT.RooArgList(mJ), hist)
        rhp = ROOT.RooHistPdf ('rhp', 'rhp', ROOT.RooArgSet (mJ), rdh)

        # Fit data with model
        pdfs[i].chi2FitTo(rdh, ROOT.RooLinkedList())

        # Plot data histogram and fit
        rhp    .plotOn(frame, ROOT.RooFit.LineColor(1), ROOT.RooFit.LineStyle(1))
        pdfs[i].plotOn(frame, ROOT.RooFit.LineColor(2), ROOT.RooFit.LineStyle(2))
        pass

    setting = ROOT.RooMomentMorph.Linear
    morph = ROOT.RooMomentMorph('morph', 'morph', mZ, ROOT.RooArgList(mJ), pdfs, massVec, setting)
    getattr(workspace,'import')(morph) # work around for morph = w.import(morph)
    morph.Print('v')

    # Make plots of interpolated pdf
    for i, mass in enumerate(massPoints[:-1] + np.diff(massPoints) / 2.):
        print i, mass
        mZ.setVal(mass)
        mZ.Print()
        morph.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(2))
        pass
    
    frame.Draw()
    ROOT.gPad.Update()        
    wait()

    return


# Main function call.
if __name__ == '__main__':
    main()
    pass
