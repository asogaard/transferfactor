# -*- coding: utf-8 -*-

""" Common utility functions related to the transfer factor method. 

@file:   utils.py
@author: Andreas Søgaard
@date:   29 March 2017
@email:  andreas.sogaard@gmail.com
"""

# Basic
import sys, itertools, os, inspect

# Scientific
try:
    # ROOT
    from ROOT import *
    from root_numpy import fill_hist
    
    # Numpy
    import numpy as np
except ImportError:
    print "WARNING: One or more scientific python packages were not found. If you're in lxplus, try running:"
    print "         $ source /cvmfs/sft.cern.ch/lcg/views/LCG_88/x86_64-slc6-gcc49-opt/setup.sh"
    sys.exit()
    pass


eps = np.finfo(float).eps


def warning (string):
    print '\033[91m\033[1mWARNING\033[0m ' + string
    return


def error (string, e=Exception):
    print '\033[91m\033[1mERROR\033[0m   ' + string
    return


def in_CR (data, mass, window):
    """ Return array indicating whether each jet is in the control region (CR). """
    if mass is None or window is None:
        return np.ones_like(data['m']).astype(bool)

    if window > 0 and window < 1:
        window *= mass 
        pass
    return (data['m'] < mass - window) | (data['m'] > mass + window)


def in_SR (data, mass, window):
    """ Return array indicating whether each jet is in the signal region (SR). """
    return ~in_CR(data, mass, window)


def get_ratio_numpy (num, denom):
    """ Return numpy array containing the mean ratio, and error, of two ROOT TH2 hitsograms. """

    # Divide ROOT histograms
    ratio = num.Clone('ratio_%d' % np.random.randint(100000))
    ratio.Divide(denom)

    # Get bin numbers
    nx, ny = num.GetXaxis().GetNbins(), num.GetYaxis().GetNbins()

    # Initialise mean and error arrays
    mean = np.zeros((ny,nx))
    err  = np.zeros((ny,nx))
        
    # Fill arrays from ROOT histogram
    for (i,j) in itertools.product(range(nx), range(ny)):
        mean[j,i] = ratio.GetBinContent(i + 1, j + 1) # Transposing due to different ordering
        err [j,i] = ratio.GetBinError  (i + 1, j + 1)
        pass

    return mean, err


def get_histogram (data, params, axes, mask=None):
    """ Utility function returning the mean and error of the Npass/Nfail profiles in the signal- and control regions. """

    # Check(s)
    if mask is None:
        mask = np.ones_like(data['weight']).astype(bool)
        pass

    # Initialise ROOT histogram
    hist = TH2F('h2_%d' % np.random.randint(100000), "", 
        len(axes[0]), axes[0][0], axes[0][-1], 
        len(axes[1]), axes[1][0], axes[1][-1])
    hist.Sumw2()

    # Prepare arrays for filling 
    x, y, w = data[params[0]], data[params[1]], data['weight']

    # Fill ROOT histogram
    matrix = np.vstack((x[mask], y[mask])).T
    fill_hist(hist, matrix, weights=w[mask])

    return hist


# Return 'iterable' in batches of (at most) 'n'.
def batch (iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]
        pass
    pass

# Use classifier tp produce predictions (wrapper function)
def predict (clf, data, eval_MSE=None):
    if eval_MSE is None:
        return clf.predict(data)
    return clf.predict(data, eval_MSE)

# Produce prediction using asynchronous processes.
from multiprocessing import Pool
import time

def asyncPredict (clf, data, num_processes=10, batch_size=10000, quiet=False, eval_MSE=None):

    # Create multiprocessing pool.
    pool = Pool()
    timeout = 9999999999

    # Get number of examples for which to produce predictions.
    num_examples = data.shape[0]

    # Compute suitable batch size.
    num_rounds = int(num_examples/float(num_processes * batch_size)) + 1

    # Loop batches.
    if not quiet: print "Total number of examples: %d (%d is maximal index)" % (num_examples, num_examples - 1)
    predictions = list()
    for iround, round_indices in enumerate(batch(np.arange(num_examples), num_processes * batch_size), start = 1):
        results = list()
        if not quiet: print "asyncPredict: Round %d of %d." % (iround, num_rounds)
        for indices in batch(round_indices, batch_size):

            # Submit prediction as asynchronous process.
            args = [clf, data[indices,:], eval_MSE]
            if num_processes > 1:
                results.append( pool.apply_async(predict, args) )
            else:
                results.append( predict(*args) )
                pass
            pass

        # Collect predictions.
        if num_processes > 1:
            predictions += [result.get(timeout = timeout) for result in results]
        else:
            predictions += results
            pass
        pass

    # Return predictions.
    return np.hstack(predictions)


def make_directories (path, fromDir=None):
    """ Create nested ROOT TDirectories (if necessary). """
    if fromDir is not None:
        fromDir.cd()
        pass

    dirs = path.split('/')

    if len(dirs) == 0:
        return

    tdir = gDirectory
    for d in dirs:
        if tdir.GetDirectory(d):
            # Exists; move there
            tdir.cd(d)
        else:
            # Doesn't exist; create it
            tdir = tdir.mkdir(d)
            tdir.cd()
            pass
        pass

    return tdir


def check_make_dir (path):
    """ Check if a system directory exists, and create it if not. """
    # Check(s)
    if path.startswith('/'):
        warning("Not accepting relative paths. Recieved: '%s'. Exiting" % path)
        return

    if not os.path.exists(path):
        os.makedirs(path)
        pass
    return


def make_serializable (iterable):
    """ Turn numpy arrays into lists, lambda functions into text, in order to make them JSON serializable. """
    result = iterable
    if   type(result) == dict:
        for key, element in result.iteritems():
            result[key] = make_serializable(element)
            pass
    elif type(result) == list:
        for idx, element in enumerate(result):
            result[idx] = make_serializable(element)
            pass
    elif type(result).__module__.startswith(np.__name__):
        result = result.tolist()
    elif type(result).__name__ == 'function':
        result = inspect.getsource(result)
        pass
    return result


def get_signal_DSID (mass, tolerance=0):
    """ From signal mass window midpoint, get DSID of corresponding signal sample. """

    # Check(s)
    #...

    # Output
    sig_DSID = None

    # Get appropriate signal file 
    theory_masses = [100, 130, 160, 190, 220]
    theory_mass = None

    for tm in theory_masses:
        # Compare input mass to theory masses within 10 GeV
        if abs(mass - tm) <= tolerance:
            theory_mass = tm
            break
        pass

    if theory_mass is None:
        print "Requested mass does not have a signal sample. Exiting."
        return sig_DSID

    #sig_file = 'objdef_MC_{DSID:6d}.root'
    if   theory_mass == 100: sig_DSID = 308363
    elif theory_mass == 130: sig_DSID = 308364
    elif theory_mass == 160: sig_DSID = 308365
    elif theory_mass == 190: sig_DSID = 308366
    elif theory_mass == 220: sig_DSID = 308367
    else: 
        print "This shouldn't happen... Exiting."
        pass
    return sig_DSID


def fixHist (h, m, w):
    """ Remove bins that overlap with the SR/VR mass window.
    @author: Laser Kaplan (@lkaplan)
    """
    if m <= 0: return h
    for j in xrange(h.GetNbinsY()):
        minrho = np.log(pow(m*(1-w), 2) / h.GetYaxis().GetBinUpEdge(j+1))
        maxrho = np.log(pow(m*(1+w), 2) / h.GetYaxis().GetBinLowEdge(j+1))
        for i in xrange(h.GetNbinsX()):
            if not (h.GetXaxis().GetBinUpEdge(i+1) < minrho or h.GetXaxis().GetBinLowEdge(i+1) > maxrho):
                h.SetBinContent(i+1, j+1, 0)
                pass
            pass
        pass
    return h
